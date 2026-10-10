"""A2 real-data replication without changing first-round artifacts or split."""
import json
import time
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans
from threadpoolctl import threadpool_limits
from scipy.stats import binomtest
from ..audit import BRANCH,DATA,load,fasta,save,sha
from ..discovery import make_pwm,background,max_scan,write_meme,group_presence,test_enrichment,bh,similarity
from ..run import logo
from .common import NEXT,SOURCES,start,finish,pool,pool_weights


def algebra():
    target,manifest=start('pooling_algebra',{'widths':[6,10,16],'tolerance':1e-7},[NEXT/'protocol.json'])
    rows=[]
    for width in [6,10,16]:
        a0,a2=pool_weights(width,'A0'),pool_weights(width,'A2')
        for position,(x,y) in enumerate(zip(a0,a2)):
            rows.append({'width':width,'start0':position,'max_weight_difference':float(abs(x-y).max()),
                         'equivalent':bool(np.allclose(x,y,atol=1e-7)),'edge_window':position<5 or position+width>76})
    path=NEXT/'pooling_ablation/weight_equivalence.tsv'
    pd.DataFrame(rows).to_csv(path,sep='\t',index=False)
    finish(target,manifest,[path],{'finding':'A0 and A2 algebraically coincide for windows away from the sequence ends because all covered bases have six covering tokens.'})


def run(source,seed,method='A2'):
    cfg=json.loads((BRANCH/'config.json').read_text(encoding='utf-8'))
    cfg['seed']=seed
    width=10
    cache=BRANCH/'cache'/source
    cachemanifest=json.loads((cache/'cache_manifest.json').read_text(encoding='utf-8'))
    tokenpath=cache/'token_L6.npy'
    assert sha(tokenpath)==cachemanifest['files'][tokenpath.name]
    discovery=load(source,'discovery',1)
    assert discovery.sequence_id.tolist()==cachemanifest['ordered_ids']
    out=NEXT/'real_data_comparison'/f'{source}__{method}__s{seed}'
    out.mkdir(exist_ok=False)
    target,manifest=start(out.name,cfg,[DATA/source/'main/sequence_master.tsv',tokenpath,NEXT/'protocol.json'])
    token=np.load(tokenpath,mmap_mode='r')
    features=pool(token,width,'A2' if method=='A2' else 'A0').reshape(-1,384)
    old=np.load(cache/'window_L6_W10.npy',mmap_mode='r')
    difference=np.linalg.norm(features.reshape(old.shape)-old,axis=-1)
    pd.DataFrame({'start0':np.arange(72),'mean_vector_l2_difference_A0_A2':difference.mean(axis=0),
                  'maximum_vector_l2_difference_A0_A2':difference.max(axis=0)}).to_csv(out/'pooling_difference.tsv',sep='\t',index=False)
    with threadpool_limits(limits=cfg['threads']):
        pca=PCA(n_components=24,svd_solver='randomized',random_state=seed)
        z=pca.fit_transform(features)
        km=MiniBatchKMeans(n_clusters=24,n_init=5,batch_size=2048,random_state=seed)
        labels=km.fit_predict(z)
        distances=np.linalg.norm(z-km.cluster_centers_[labels],axis=1)
    motifs,instances,merges={},[],[]
    for cluster in range(24):
        idxs=np.flatnonzero(labels==cluster)
        selected={}
        for idx in idxs[np.argsort(distances[idxs],kind='stable')]:
            selected.setdefault(int(idx//72),int(idx))
        if len(selected)<20:continue
        name=f'{method}_C{cluster:02d}'
        fragments=[]
        for row,idx in selected.items():
            pos=idx%72;r=discovery.iloc[row];seq=r.sequence[pos:pos+width];fragments.append(seq)
            instances.append({'motif':name,'sequence_id':r.sequence_id,'leakage_group':r.leakage_group,
                  'start0':pos,'end0_exclusive':pos+width,'start1':pos+1,'stop1':pos+width,
                  'center_tss':pos+(width-1)/2-60,'strand':'+','sequence':seq,'split':'discovery'})
        pfm,pwm=make_pwm(fragments)
        other=next((n for n,m in motifs.items() if similarity(pwm,m,width)[0]>=.95),None)
        np.savetxt(out/f'{name}.pfm.tsv',pfm,delimiter='\t')
        if other:merges.append({'motif':name,'representative':other})
        else:motifs[name]=pwm;logo(out/f'{name}.png',pwm,name)
    natural_disc=load(source,'discovery',0);dev=load(source,'development',1);natural_dev=load(source,'development',0)
    parents=dict(zip(dev.sequence_id,dev.leakage_group))
    shuffled=pd.DataFrame([{'sequence_id':i,'sequence':s,'leakage_group':parents[i.removesuffix('__dinucleotide_null')]}
                   for i,s in fasta(DATA/source/'main/development.dinucleotide_null.fasta')])
    bg=background(natural_disc.sequence)
    write_meme(out/'motifs.meme',motifs,bg)
    save(out/'merge_map.json',merges)
    pd.DataFrame(instances).to_csv(out/'instances.tsv',sep='\t',index=False)
    hits,evaluation,thresholds=[],[],{}
    for name,pwm in motifs.items():
        threshold=float(np.quantile([x[0] for x in max_scan(natural_disc.sequence,pwm,bg)],.95,method='higher'))
        thresholds[name]=threshold
        presence={}
        for control,frame in [('positive',dev),('natural',natural_dev),('shuffle',shuffled)]:
            scanned=max_scan(frame.sequence,pwm,bg)
            presence[control]=group_presence(frame,[x[0]>threshold for x in scanned])
            for r,(score,pos,strand) in zip(frame.itertuples(),scanned):
                hits.append({'motif':name,'control':control,'sequence_id':r.sequence_id,'leakage_group':r.leakage_group,
                   'score':score,'present':score>threshold,'start1':pos+1,'stop1':pos+width,'strand':strand,'center_tss':pos+(width-1)/2-60})
        for control in ['natural','shuffle']:
            common=set(presence['positive'].index)&set(presence[control].index)
            if control=='natural':
                stats=test_enrichment(presence['positive'].drop(list(common)),presence[control].drop(list(common)))
                stats['shared_groups_excluded']=len(common)
            else:
                a=presence['positive'].loc[sorted(common)].to_numpy();b=presence['shuffle'].loc[sorted(common)].to_numpy()
                discord_pos=int((a&~b).sum());discord_null=int((~a&b).sum())
                stats={'positive_groups':len(a),'control_groups':len(b),'positive_hits':int(a.sum()),'control_hits':int(b.sum()),
                       'discordant_pos_only':discord_pos,'discordant_null_only':discord_null,
                       'p':float(binomtest(discord_pos,discord_pos+discord_null,.5,alternative='greater').pvalue) if discord_pos+discord_null else 1.}
            evaluation.append({'method':method,'motif':name,'control':control,'threshold':threshold,**stats})
    ev=pd.DataFrame(evaluation)
    for control in ['natural','shuffle']:
        mask=ev.control==control;ev.loc[mask,'q_run']=bh(ev.loc[mask,'p'])
    ev.to_csv(out/'development_enrichment.tsv',sep='\t',index=False)
    pd.DataFrame(hits).to_csv(out/'development_scans.tsv',sep='\t',index=False)
    save(out/'thresholds.json',thresholds)
    # No changes to old run manifests; the second-round manifest captures all new artifacts.
    finish(target,manifest,[p for p in out.iterdir() if p.is_file()],{'source':source,'seed':seed,'method':method,'motifs':len(motifs),
              'source_cache_model_revision':cachemanifest['config']['revision'],'discovery_n':len(discovery),'development_n':len(dev)})
    print(f'Completed {out.name}, motifs={len(motifs)}',flush=True)


def main():
    algebra()
    for source in SOURCES:
        for seed in [20261005,20261006,20261007]:run(source,seed)


if __name__=='__main__':main()
