"""Bounded supervised CNN candidate discovery on frozen real promoter data."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
import torch

from ..audit import load, fasta
from ..cnn_boundary_20261027.pipeline import LocalCNN, choose_shift, shifted_fragment, direction_maxima, setup, train
from ..discovery import ALPHABET, background, make_pwm, read_meme, similarity, write_meme
from ..refinement_20261026.pipeline import sequence_windows, scan_windows

HERE=Path(__file__).resolve().parent; ROOT=HERE.parents[2]; ART=HERE/'artifacts'; SUMMARY=HERE/'summary'
DATA=ROOT/'data/processed/pr01_02_data_v2/tjupan_bacillus_subtilis/main'
STREME=ROOT/'results/motif/tju_streme_v2_20261008/tjupan_bacillus_subtilis/selected_main.meme'
PROTOCOL=HERE/'protocol.json'

def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()

def read(p):return json.loads(p.read_text(encoding='utf8'))

def save(p,x):p.write_text(json.dumps(x,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')

def groups(frame,flags):
    return pd.DataFrame({'group':frame.leakage_group.to_numpy(),'hit':flags}).groupby('group').hit.max()

def load_split(split):
    pos=load('tjupan_bacillus_subtilis',split,1);neg=load('tjupan_bacillus_subtilis',split,0)
    assert len(pos) and len(neg) and pos.sequence.str.len().eq(81).all() and neg.sequence.str.len().eq(81).all()
    overlap=set(pos.leakage_group)&set(neg.leakage_group)
    assert not pos.sequence_id.duplicated().any() and not neg.sequence_id.duplicated().any()
    records=fasta(DATA/f'{split}.dinucleotide_null.fasta'); mapping={k.removesuffix('__dinucleotide_null'):s for k,s in records}
    assert len(mapping)==len(records)==len(pos) and set(mapping)==set(pos.sequence_id)
    return pos,neg,[mapping[k] for k in pos.sequence_id],overlap

def rates(pos,neg,positive_hit,negative_hit,overlap):
    p=groups(pos,positive_hit).drop(index=list(overlap),errors='ignore')
    n=groups(neg,negative_hit).drop(index=list(overlap),errors='ignore')
    assert len(p) and len(n)
    return float(p.mean()),float(n.mean()),len(p),len(n)

def rank(motifs,pos,neg,bg,q,overlap):
    caches={};r=[]
    for name,pwm in motifs.items():
        w=len(pwm)
        if w not in caches:caches[w]=(sequence_windows(pos.sequence.tolist(),w),sequence_windows(neg.sequence.tolist(),w))
        a,b=(scan_windows(x,pwm,bg)[0] for x in caches[w]);threshold=float(np.quantile(b,q,method='higher'))
        pr,nr,_,_=rates(pos,neg,a>threshold,b>threshold,overlap)
        info=float(np.sum(pwm*np.log2(np.maximum(pwm,1e-12)/bg)))
        r.append({'candidate':name,'discovery_score':pr-nr,'information':info,'threshold':threshold,'width':w})
    return sorted(r,key=lambda x:(-x['discovery_score'],-x['information'],x['candidate']))

def fit():
    cfg=read(PROTOCOL);setup(cfg);pos,neg,_,overlap=load_split('discovery')
    baseline_manifest=ROOT/'experiments/pr01_02_embedding_v2/runs/tjupan_bacillus_subtilis__prokbert__L6W10__s20261005/manifest.json'
    previous=read(baseline_manifest)
    for p in [DATA/'sequence_master.tsv',STREME]:
        key=str(p.relative_to(ROOT)).replace('/','\\')
        assert sha(p)==previous['input_hashes'][key],f'Frozen input hash mismatch: {p}'
    ART.mkdir(exist_ok=False);SUMMARY.mkdir(exist_ok=False)
    bg=background(neg.sequence.tolist()); sources=[PROTOCOL,Path(__file__),baseline_manifest,DATA/'sequence_master.tsv',DATA/'discovery.dinucleotide_null.fasta',DATA/'development.dinucleotide_null.fasta',STREME,
      ROOT/'experiments/pr01_02_embedding_v2/cnn_features_20261026/pipeline.py',ROOT/'experiments/pr01_02_embedding_v2/cnn_boundary_20261027/pipeline.py',ROOT/'experiments/pr01_02_embedding_v2/discovery.py',ROOT/'experiments/pr01_02_embedding_v2/refinement_20261026/pipeline.py']
    source_hashes={p.relative_to(ROOT).as_posix():sha(p) for p in sources}
    baseline=read_meme(STREME);assert len(baseline)==3
    save(ART/'streme_ranking.json',rank(baseline,pos,neg,bg,cfg['scan_control_quantile'],overlap))
    for seed in cfg['algorithm_seeds']:
        out=ART/f'seed_{seed}';out.mkdir();torch.manual_seed(seed);model=LocalCNN(cfg['cnn']['filters'],cfg['cnn']['width'])
        curve,label_hash=train(model,pos.sequence.tolist()+neg.sequence.tolist(),[1]*len(pos)+[0]*len(neg),cfg['cnn'],seed)
        torch.save(model.state_dict(),out/'model.pt');save(out/'training.json',{'curve':curve,'label_sha256':label_hash})
        score,start,forward=direction_maxima(model,pos.sequence.tolist(),'TIED');null,_,_=direction_maxima(model,neg.sequence.tolist(),'TIED')
        thresholds=np.quantile(null,cfg['activation_control_quantile'],axis=0,method='higher'); pwm={'RAW':{},'SHIFT':{}};instances={'RAW':[],'SHIFT':[]};decisions=[]
        for j in range(score.shape[1]):
            ids=np.flatnonzero(score[:,j]>thresholds[j]);name=f'F{j:02d}'
            selected=[{'sequence_index':int(i),'start0':int(start[i,j]),'forward':bool(forward[i,j]),'activation':float(score[i,j])} for i in ids]
            eligible,shift,values=choose_shift(pos.sequence.tolist(),selected,bg,cfg['cnn']['width'])
            kept=len(eligible)>=cfg['minimum_support'];decisions.append({'candidate':name,'selected':len(ids),'eligible':len(eligible),'retained':kept,'shift':shift,'information_by_shift':values})
            if not kept:continue
            for variant,s in [('RAW',0),('SHIFT',shift)]:
                for item in eligible:
                    begin,fragment=shifted_fragment(pos.sequence.iloc[item['sequence_index']],item['start0'],item['forward'],s,cfg['cnn']['width'])
                    instances[variant].append({'candidate':name,'sequence_id':pos.sequence_id.iloc[item['sequence_index']],'start0':begin,'strand':'+' if item['forward'] else '-','fragment':fragment})
                _,pwm[variant][name]=make_pwm([r['fragment'] for r in instances[variant] if r['candidate']==name])
        save(out/'decisions.json',decisions)
        for variant in ['RAW','SHIFT']:
            d=out/variant;d.mkdir();write_meme(d/'motifs.meme',pwm[variant],bg)
            pd.DataFrame(instances[variant],columns=['candidate','sequence_id','start0','strand','fragment']).to_csv(d/'instances.tsv',sep='\t',index=False)
            save(d/'ranking.json',rank(read_meme(d/'motifs.meme'),pos,neg,bg,cfg['scan_control_quantile'],overlap))
        print(f'fit seed {seed}: RAW {len(pwm["RAW"])} SHIFT {len(pwm["SHIFT"])}',flush=True)
    save(SUMMARY/'fit_lock.json',{'source_sha256':source_hashes,'artifact_sha256':{x.relative_to(HERE).as_posix():sha(x) for x in ART.rglob('*') if x.is_file()},
      'positive_n':len(pos),'natural_control_n':len(neg),'mixed_groups_excluded':len(overlap),'holdout_used':False,'fit_complete':True})
    print('FIT_LOCKED',flush=True)

def validate_lock():
    lock=read(SUMMARY/'fit_lock.json')
    for p,d in lock['source_sha256'].items():assert sha(ROOT/p)==d,p
    for p,d in lock['artifact_sha256'].items():assert sha(HERE/p)==d,p
    return lock

def score():
    validate_lock();cfg=read(PROTOCOL);disc,negative,_,_=load_split('discovery');pos,neg,shuffle,overlap=load_split('development');bg=background(negative.sequence.tolist())
    refs=read_meme(STREME);sets=[('STREME','existing',refs,read(ART/'streme_ranking.json'))]
    for seed in cfg['algorithm_seeds']:
        for v in ['RAW','SHIFT']:
            d=ART/f'seed_{seed}'/v;sets.append(('CNN_'+v,str(seed),read_meme(d/'motifs.meme'),read(d/'ranking.json')))
    results=[];candidate=[]
    for method,seed,motifs,ranking in sets:
        hits={};windows={}
        for item in ranking:
            name=item['candidate'];pwm=motifs[name];w=len(pwm)
            if w not in windows:windows[w]=[sequence_windows(x,w) for x in [pos.sequence.tolist(),neg.sequence.tolist(),shuffle]]
            hits[name]=[scan_windows(x,pwm,bg)[0]>item['threshold'] for x in windows[w]]
            pr,nr,pn,nn=rates(pos,neg,hits[name][0],hits[name][1],overlap)
            pg=groups(pos,hits[name][0]);sg=groups(pos,hits[name][2])
            best=max(((similarity(pwm,r,minimum_overlap=min(6,len(r),len(pwm))),n) for n,r in refs.items()),key=lambda x:x[0][0])
            candidate.append({'method':method,'seed':seed,**item,'positive_groups':pn,'natural_groups':nn,'dev_positive_rate':pr,'dev_natural_rate':nr,'dev_shuffle_rate':float(sg.mean()),
                'nearest_streme_id':best[1],'streme_cosine':float(best[0][0]),'streme_offset':int(best[0][1]),'streme_orientation':best[0][2]})
        for view,names in [('top3',[x['candidate'] for x in ranking[:3]]),('all',[x['candidate'] for x in ranking])]:
            if names:
                hit=[np.logical_or.reduce([hits[n][i] for n in names]) for i in range(3)]
            else:hit=[np.zeros(len(pos),bool),np.zeros(len(neg),bool),np.zeros(len(pos),bool)]
            pr,nr,pn,nn=rates(pos,neg,hit[0],hit[1],overlap);pg=groups(pos,hit[0]);sg=groups(pos,hit[2])
            results.append({'method':method,'seed':seed,'view':view,'candidates':len(names),'positive_groups':pn,'natural_groups':nn,
                'positive_rate':pr,'natural_rate':nr,'natural_difference':pr-nr,'paired_positive_rate':float(pg.mean()),'shuffle_rate':float(sg.mean()),
                'paired_positive_only':int((pg&~sg).sum()),'paired_shuffle_only':int((sg&~pg).sum()),'paired_groups':len(pg)})
    pd.DataFrame(candidate).to_csv(SUMMARY/'candidates.tsv',sep='\t',index=False,lineterminator='\n')
    pd.DataFrame(results).to_csv(SUMMARY/'family.tsv',sep='\t',index=False,lineterminator='\n')
    save(SUMMARY/'score_manifest.json',{'status':'COMPLETED','fit_lock_sha256':sha(SUMMARY/'fit_lock.json'),
        'outputs':{x.name:sha(x) for x in [SUMMARY/'candidates.tsv',SUMMARY/'family.tsv']},'development_positive':len(pos),'development_natural':len(neg),
        'development_mixed_groups':len(overlap),'holdout_used':False,'evaluation':'retrospective Development; same-data STREME reference'})
    print(pd.DataFrame(results).to_string(index=False),flush=True)

def verify():
    lock=validate_lock();m=read(SUMMARY/'score_manifest.json');assert m['fit_lock_sha256']==sha(SUMMARY/'fit_lock.json')
    for n,h in m['outputs'].items():assert sha(SUMMARY/n)==h
    c=pd.read_csv(SUMMARY/'candidates.tsv',sep='\t');f=pd.read_csv(SUMMARY/'family.tsv',sep='\t')
    assert len(f)==10 and set(f.method)=={'STREME','CNN_RAW','CNN_SHIFT'}
    assert (f[['positive_rate','natural_rate','shuffle_rate']]>=0).all().all() and (f[['positive_rate','natural_rate','shuffle_rate']]<=1).all().all()
    assert len(list(ART.glob('seed_*/model.pt')))==2
    for seed in read(PROTOCOL)['algorithm_seeds']:
        for variant in ['RAW','SHIFT']:
            d=ART/f'seed_{seed}'/variant;motifs=read_meme(d/'motifs.meme');instances=pd.read_csv(d/'instances.tsv',sep='\t')
            assert set(motifs)==set(instances.candidate)
            for name,pwm in motifs.items():
                _,rebuilt=make_pwm(instances[instances.candidate==name].fragment.tolist());assert np.allclose(rebuilt,pwm,atol=1e-8)
    save(SUMMARY/'verification.json',{'status':'PASSED','model_count':2,'candidate_rows':len(c),'family_rows':len(f),
        'source_hashes':len(lock['source_sha256']),'locked_files':len(lock['artifact_sha256']),'pwm_instance_reconstruction':'passed'})
    print('VERIFIED',len(c),'candidate rows',flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['fit','score','verify']);args=parser.parse_args();{'fit':fit,'score':score,'verify':verify}[args.action]()
