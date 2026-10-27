"""Fit CNN on exact existing BERT batches, then apply the same retrospective scoring."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
import torch

from ..pr01_02_embedding_v2.budget_validation_20261026.pipeline import candidate_cache
from ..pr01_02_embedding_v2.cnn_boundary_20261027.pipeline import LocalCNN, cnn_pools, setup, train
from ..pr01_02_embedding_v2.discovery import ALPHABET, read_meme, similarity
from ..pr01_02_embedding_v2.refinement_20261026.pipeline import sequence_windows
from ..practice_motif_20261027.build import evaluate_partial

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE=ROOT/'experiments/pr01_02_embedding_v2/budget_validation_20261026'
ART=HERE/'artifacts'
OUT=ROOT/'reports/same_batch_motif_20261027'

def read(p):
    return json.loads(p.read_text(encoding='utf8'))

def sha(p):
    with p.open('rb') as f: return hashlib.file_digest(f,'sha256').hexdigest()

def save(p,obj):
    p.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n',encoding='utf8')

def config():
    p=read(HERE/'protocol.json'); c=read(BASE/'protocol.json'); c['cnn']=p['cnn']; c['threads']=p['threads']
    return p,c

def checked(path,lock):
    rel=path.relative_to(BASE).as_posix()
    if sha(path)!=lock['outputs'].get(rel): raise ValueError(f'Prior locked source differs: {rel}')
    return path

def fit():
    p,c=config(); setup(c); lock=read(BASE/'summary/fit_lock.json')
    for name,digest in lock['source_sha256'].items():
        if sha(ROOT/name)!=digest: raise ValueError(f'Prior source changed: {name}')
    ART.mkdir(exist_ok=False)
    source={}; done=0
    for scenario in p['scenario_ids']:
        for seed in p['data_seeds']:
            name=f'{scenario}__s{seed}'; old=BASE/'artifacts'/name; folder=ART/name;folder.mkdir()
            data_path=checked(old/'data.json',lock); data=read(data_path)
            source[data_path.relative_to(ROOT).as_posix()]=sha(data_path)
            for alg in p['algorithm_seeds']:
                parent=folder/f'a{alg}';parent.mkdir();torch.manual_seed(alg)
                model=LocalCNN(c['cnn']['filters'],c['width'])
                curve,label_hash=train(model,data['discovery']+data['discovery_control'],
                    [1]*len(data['discovery'])+[0]*len(data['discovery_control']),c['cnn'],alg)
                torch.save(model.state_dict(),parent/'model.pt')
                ranks=cnn_pools(model,data['discovery'],data['discovery_control'],c,'TIED',parent)
                save(parent/'fit_metadata.json',{'algorithm_seed':alg,'label_sha256':label_hash,'curve':curve})
                for method, ranking in ranks.items():
                    save(parent/method/'cap5.json',[r['candidate'] for r in ranking[:p['candidate_budget']]])
                for method in ['EMB128','OH24','OH48']:
                    for filename in ['ranking.json','motifs.meme']:
                        q=checked(old/f'{method}__a{alg}'/'BG_ONLY'/filename,lock)
                        source[q.relative_to(ROOT).as_posix()]=sha(q)
                done+=1
            print(f'fitted {name}: {done}/45 CNN models',flush=True)
    output={x.relative_to(HERE).as_posix():sha(x) for x in ART.rglob('*') if x.is_file()}
    save(HERE/'fit_lock.json',{'protocol_sha256':sha(HERE/'protocol.json'),'pipeline_sha256':sha(HERE/'pipeline.py'),
        'prior_lock_sha256':sha(BASE/'summary/fit_lock.json'),'source_sha256':source,'output_sha256':output,
        'cnn_fits':done,'test_note':'existing batch was seen in earlier experiments'})
    print('FIT_LOCKED',flush=True)

def score():
    p,c=config(); lock=read(HERE/'fit_lock.json')
    if sha(HERE/'protocol.json')!=lock['protocol_sha256'] or sha(HERE/'pipeline.py')!=lock['pipeline_sha256']:
        raise ValueError('Protocol or pipeline changed after fitting')
    for rel,digest in lock['source_sha256'].items():
        if sha(ROOT/rel)!=digest: raise ValueError(f'Source changed: {rel}')
    for rel,digest in lock['output_sha256'].items():
        if sha(HERE/rel)!=digest: raise ValueError(f'CNN artifact changed: {rel}')
    rows=[]
    for scenario in p['scenario_ids']:
        for seed in p['data_seeds']:
            name=f'{scenario}__s{seed}'; old=BASE/'artifacts'/name; folder=ART/name;data=read(old/'data.json')
            prepared={k:sequence_windows(data[k],c['width']) for k in ['calibration_control','development','development_control','test','test_control']}
            reference=np.full((c['width'],4),data['scenario']['mutation']/3)
            for i,char in enumerate(c['synthetic_pattern']):reference[i,ALPHABET.index(char)]=1-data['scenario']['mutation']
            for alg in p['algorithm_seeds']:
                for method in p['methods']:
                    if method.startswith('TIED'):
                        directory=folder/f'a{alg}'/method; names=read(directory/'cap5.json')
                    else:
                        directory=old/f'{method}__a{alg}'/'BG_ONLY'
                        names=[r['candidate'] for r in read(directory/'ranking.json')[:p['candidate_budget']]]
                    motifs=read_meme(directory/'motifs.meme'); cache,_=candidate_cache(motifs,data,c,prepared)
                    matches={n:similarity(pwm,reference,minimum_overlap=p['partial_min_overlap_bp']) for n,pwm in motifs.items()}
                    for split in ['development','test']:
                        rows.append({'scenario':scenario,'data_seed':seed,'algorithm_seed':alg,'method':method,
                            **evaluate_partial(cache,matches,names,data,split,p['position_tolerance_bp'])})
            print(f'scored {name}',flush=True)
    OUT.mkdir(parents=True,exist_ok=False)
    frame=pd.DataFrame(rows);assert len(frame)==450
    frame.to_csv(OUT/'runs.tsv',sep='\t',index=False,lineterminator='\n')
    per=frame.groupby(['scenario','data_seed','method','split'],as_index=False).agg(f1=('f1','mean'),
        blank_fpr=('family_blank_fpr','mean'),candidate_count=('candidate_count','mean'),overlap_bp=('overlap_bp','mean'))
    per.to_csv(OUT/'per_data_seed.tsv',sep='\t',index=False,lineterminator='\n')
    summary=per.groupby(['scenario','method','split'],as_index=False).agg(f1_mean=('f1','mean'),
        f1_min=('f1','min'),f1_max=('f1','max'),blank_fpr=('blank_fpr','mean'),n_data=('data_seed','nunique'))
    summary.to_csv(OUT/'summary.tsv',sep='\t',index=False,lineterminator='\n')
    save(OUT/'manifest.json',{'protocol':p,'fit_lock_sha256':sha(HERE/'fit_lock.json'),
        'run_rows':len(frame),'data_seed_rows':len(per),'summary_rows':len(summary),
        'outputs_sha256':{x.name:sha(x) for x in OUT.iterdir() if x.is_file()}})
    print('SCORE_COMPLETED',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['fit','score']);args=parser.parse_args()
    {'fit':fit,'score':score}[args.action]()
