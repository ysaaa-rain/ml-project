"""End-to-end toy discovery recovery. Ground truth is evaluation-only."""
import json
import time
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans
from threadpoolctl import threadpool_limits
from .audit import BRANCH, save, sha
from .encoder import load_model, extract
from .discovery import make_pwm, background, max_scan, similarity, test_enrichment, bh, write_meme


def main(tag='synthetic'):
    out = BRANCH/tag
    out.mkdir(exist_ok=False)
    cfg = json.loads((BRANCH/'config.json').read_text(encoding='utf-8'))
    rng = np.random.default_rng(cfg['seed'])
    motif = 'TTGACATATAAT'
    def generate(n,plant):
        sequences = [''.join(rng.choice(list('ACGT'),81)) for _ in range(n)]
        if plant:
            sequences = [s[:25]+motif+s[37:] for s in sequences]
        return sequences
    disc,dev,calibration,control = generate(256,True),generate(128,True),generate(256,False),generate(128,False)
    frame = pd.DataFrame({'sequence_id':[f'synthetic_{i}' for i in range(256)],'sequence':disc,
                          'leakage_group':[f'synthetic_{i}' for i in range(256)],'split':'discovery'})
    frame.to_csv(out/'discovery.tsv',sep='\t',index=False)
    save(out/'input.json',{'motif':motif,'start0':25,'discovery':disc,'development':dev,
         'calibration_control':calibration,'development_control':control,'seed':cfg['seed'],
         'scope':'strong planted 12bp toy; cannot establish biological functionality or weak-motif ability'})
    model,tokenizer = load_model(cfg)
    extract(frame,out/'cache',cfg,model,tokenizer)
    bg = background(calibration)
    _,truth = make_pwm([motif]*256)
    rows, instances = [], []
    started = time.time()
    for method in ['prokbert','tokenwindow','onehot']:
        width=10; n_windows=72
        if method=='prokbert':
            features=np.load(out/'cache/window_L6_W10.npy').reshape(-1,384)
        elif method=='tokenwindow':
            token=np.load(out/'cache/token_L6.npy')
            features=np.stack([token[:,i:i+5].mean(axis=1) for i in range(72)],axis=1).reshape(-1,384)
        else:
            features=np.stack([np.eye(4,dtype=np.float32)[['ACGT'.index(c) for c in s[i:i+width]]].ravel()
                              for s in disc for i in range(n_windows)])
        with threadpool_limits(limits=cfg['threads']):
            z=PCA(n_components=24,svd_solver='randomized',random_state=cfg['seed']).fit_transform(features)
            k=MiniBatchKMeans(n_clusters=24,n_init=5,batch_size=2048,random_state=cfg['seed'])
            labels=k.fit_predict(z)
            distance=np.linalg.norm(z-k.cluster_centers_[labels],axis=1)
        motifs={}
        for cluster in range(24):
            ids=np.flatnonzero(labels==cluster)
            selected={}
            for idx in ids[np.argsort(distance[ids],kind='stable')]:
                selected.setdefault(int(idx//n_windows),int(idx%n_windows))
            if len(selected)<20:
                continue
            fragments=[disc[i][s:s+width] for i,s in selected.items()]
            _,pwm=make_pwm(fragments)
            name=f'{method}_C{cluster:02d}'
            motifs[name]=pwm
            purity=float(np.mean([max(0,min(s+width,37)-max(s,25))>=8 for s in selected.values()]))
            for i,s in selected.items():
                instances.append({'motif':name,'sequence_id':f'synthetic_{i}','start0':s,'end0':s+width,'sequence':disc[i][s:s+width]})
            sim,offset,strand=similarity(pwm,truth,minimum_overlap=10)
            threshold=float(np.quantile([a[0] for a in max_scan(calibration,pwm,bg)],.95,method='higher'))
            pos=np.array([a[0]>threshold for a in max_scan(dev,pwm,bg)])
            neg=np.array([a[0]>threshold for a in max_scan(control,pwm,bg)])
            rows.append({'method':method,'motif':name,'localization_purity_overlap_ge8':purity,
                         'truth_fullwidth_similarity':sim,'offset':offset,'strand':strand,**test_enrichment(pos,neg)})
        write_meme(out/f'{method}.meme',motifs,bg)
    results=pd.DataFrame(rows)
    results['q']=bh(results.p)
    results['recovered']=(results.localization_purity_overlap_ge8>=.5)&(results.truth_fullwidth_similarity>=.8)&(results.q<=.05)
    results.to_csv(out/'recovery.tsv',sep='\t',index=False)
    pd.DataFrame(instances).to_csv(out/'instances.tsv',sep='\t',index=False)
    save(out/'manifest.json',{'status':'COMPLETED','config':cfg,'elapsed_post_encoding_seconds':time.time()-started,
         'recovery_criteria':'PWM cosine >=0.8 over >=10bp, localization overlap >=8 purity >=0.5, dev q<=0.05',
         'recovered':results.groupby('method').recovered.sum().astype(int).to_dict(),
         'code_sha256':{p.name:sha(p) for p in BRANCH.glob('*.py')},
         'output_sha256':{p.name:sha(p) for p in out.iterdir() if p.is_file()}})
    print(results.groupby('method').recovered.sum().to_string(),flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--tag',default='synthetic')
    main(parser.parse_args().tag)
