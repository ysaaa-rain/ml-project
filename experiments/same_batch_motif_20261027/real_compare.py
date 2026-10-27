"""Uniform sequence scanning of existing real-data candidate PWMs."""
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from ..pr01_02_embedding_v2.audit import load
from ..pr01_02_embedding_v2.discovery import background, read_meme
from ..pr01_02_embedding_v2.refinement_20261026.pipeline import sequence_windows, scan_windows

ROOT=Path(__file__).resolve().parents[2]
BRANCH=ROOT/'experiments/pr01_02_embedding_v2'
STREME=ROOT/'results/motif/tju_streme_v2_20261008'
OUT=ROOT/'reports/same_batch_real_equal_20261027'
SEEDS=[20261005,20261006,20261007]

def sha(path):
    with path.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def main():
    sources=sorted(p.parent.name for p in STREME.glob('tjupan_*/selected_main.meme'))
    assert len(sources)==6
    provenance={};candidate_rows=[];run_rows=[]
    for source in sources:
        data_path=ROOT/f'data/processed/pr01_02_data_v2/{source}/main/sequence_master.tsv'
        positive=load(source,'discovery',1).sequence.tolist();blank=load(source,'discovery',0).sequence.tolist()
        devpos=load(source,'development',1).sequence.tolist();devblank=load(source,'development',0).sequence.tolist()
        assert min(map(len,[positive,blank,devpos,devblank]))>0
        provenance[data_path.relative_to(ROOT).as_posix()]=sha(data_path)
        bg=background(blank);windows={}
        pools=[('STREME','original',STREME/source/'selected_main.meme',None)]
        for method in ['onehot','prokbert']:
            for seed in SEEDS:
                folder=BRANCH/'runs'/f'{source}__{method}__L6W10__s{seed}'
                meta=json.loads((folder/'manifest.json').read_text(encoding='utf8'))
                assert meta['status']=='COMPLETED' and not meta['holdout_used']
                path=folder/'motifs.meme'; assert sha(path)==meta['output_sha256']['motifs.meme']
                key=str(data_path.relative_to(ROOT)).replace('/','\\')
                assert sha(data_path)==meta['input_hashes'][key]
                provenance[(folder/'manifest.json').relative_to(ROOT).as_posix()]=sha(folder/'manifest.json')
                pools.append((method,str(seed),path,meta))
        common_budget=min(5,*(len(read_meme(path)) for _,_,path,_ in pools))
        if common_budget<1:raise ValueError(f'No common candidates for {source}')
        for method,seed,path,meta in pools:
            provenance[path.relative_to(ROOT).as_posix()]=sha(path)
            pwms=read_meme(path); assert pwms
            rankings=[]
            for name,pwm in pwms.items():
                width=len(pwm)
                if width not in windows:
                    windows[width]=[sequence_windows(x,width) for x in [positive,blank,devpos,devblank]]
                scores=[scan_windows(w,pwm,bg)[0] for w in windows[width]]
                threshold=float(np.quantile(scores[1],.95,method='higher'))
                hit=[v>threshold for v in scores]
                score=float(hit[0].mean()-hit[1].mean())
                information=float((pwm*np.log2(np.maximum(pwm,1e-12)/bg)).sum())
                rankings.append({'name':name,'width':width,'score':score,'information':information,'threshold':threshold,'hit':hit})
            rankings.sort(key=lambda r:(-r['score'],-r['information'],r['name']))
            chosen=rankings[:common_budget]
            for rank,r in enumerate(chosen,1):
                candidate_rows.append({'source':source,'method':method,'seed':seed,'candidate':r['name'],'rank':rank,
                    'width':r['width'],'discovery_presence_difference':r['score'],'threshold':r['threshold'],
                    'dev_positive_rate':float(r['hit'][2].mean()),'dev_blank_rate':float(r['hit'][3].mean())})
            pos=np.any([r['hit'][2] for r in chosen],axis=0);neg=np.any([r['hit'][3] for r in chosen],axis=0)
            run_rows.append({'source':source,'method':method,'seed':seed,'candidate_count':len(chosen),
                'positive_n':len(devpos),'control_n':len(devblank),'positive_hits':int(pos.sum()),
                'control_hits':int(neg.sum()),'common_budget':common_budget,'positive_rate':float(pos.mean()),'control_rate':float(neg.mean()),
                'rate_difference':float(pos.mean()-neg.mean())})
        print('completed',source,flush=True)
    OUT.mkdir(parents=True,exist_ok=False)
    candidates=pd.DataFrame(candidate_rows);runs=pd.DataFrame(run_rows)
    assert len(runs)==42 and not runs.duplicated(['source','method','seed']).any()
    candidates.to_csv(OUT/'candidates.tsv',sep='\t',index=False,lineterminator='\n')
    runs.to_csv(OUT/'runs.tsv',sep='\t',index=False,lineterminator='\n')
    summary=runs.groupby(['source','method'],as_index=False).agg(n_runs=('seed','size'),candidate_count=('candidate_count','mean'),
        positive_rate=('positive_rate','mean'),control_rate=('control_rate','mean'),rate_difference=('rate_difference','mean'))
    summary.to_csv(OUT/'summary.tsv',sep='\t',index=False,lineterminator='\n')
    files={x.name:sha(x) for x in OUT.iterdir() if x.is_file()}
    manifest={'status':'COMPLETED','source_sha256':provenance,'script_sha256':sha(Path(__file__)),
        'output_sha256':files,'count_runs':len(runs),'split':'Development','source_method':'Discovery-selected existing PWMs',
        'ranking':'Discovery positive-control rate difference; ties information then name; same min(5, available across all runs) for each species',
        'scan':'same background per species, two strand log odds, each PWM threshold 95th percentile Discovery control, strict greater',
        'limits':'STREME original selected main subset; Embedding runs L6W10; sequence-level family OR; previously observed Development; CNN unavailable'}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(summary.to_string(index=False),flush=True)

if __name__=='__main__':main()
