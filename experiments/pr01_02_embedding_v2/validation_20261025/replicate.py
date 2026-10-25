"""Unmodified prior fit recipe, independently locked and reported."""
import argparse
import time

import numpy as np
import pandas as pd

from ..bottleneck_20261025.pilot import fit as original_fit
from .pipeline import (ART,HERE,OLD,PROTOCOL,SUMMARY,assert_locked,background,evaluate,
                       onehot,read,read_meme,save,sha,write_meme)

RECIPE=HERE/'replication_protocol.json'
DEST=HERE/'replication_artifacts'


def fit():
    assert_locked(); cfg=read(PROTOCOL); DEST.mkdir(exist_ok=False)
    manifest={'status':'RUNNING','started_unix':time.time(),'protocol_sha256':sha(PROTOCOL),
              'recipe_sha256':sha(RECIPE),'source_sha256':sha(HERE/'replicate.py'),
              'original_fit_source_sha256':sha(OLD/'pilot.py'),'runs':[],'real_holdout_used':False}
    save(SUMMARY/'replication_manifest.json',manifest)
    try:
        for scenario in cfg['scenarios']:
            for seed in cfg['data_seeds']:
                condition=f'{scenario["id"]}__s{seed}'; data=read(ART/condition/'data.json')
                features=np.load(ART/condition/'local_discovery.npy',mmap_mode='r')
                oh=onehot(data['discovery'],cfg['width']); folder=DEST/condition; folder.mkdir()
                for method in read(RECIPE)['methods']:
                    for algorithm_seed in cfg['algorithm_seeds']:
                        out=folder/f'{method}__a{algorithm_seed}'; out.mkdir()
                        fitted_cfg={**cfg,'algorithm_seed':algorithm_seed,'pca_components':24}
                        motifs,variance=original_fit(oh if method.startswith('OH') else features,
                            data['discovery'],fitted_cfg,out)
                        write_meme(out/'motifs.meme',motifs,background(data['discovery_control']))
                        item={'condition':condition,'scenario':scenario['id'],'data_seed':seed,
                              'algorithm_seed':algorithm_seed,'method':method,'candidate_count':len(motifs),
                              'pca_variance':variance,'status':'COMPLETED'}
                        manifest['runs'].append(item)
                save(SUMMARY/'replication_manifest.json',manifest)
                print(f'recipe replication {condition}: {len(manifest["runs"])}/90',flush=True)
                del features,oh
        lock={'protocol_sha256':sha(PROTOCOL),'recipe_sha256':sha(RECIPE),'source_sha256':sha(HERE/'replicate.py'),
              'original_fit_source_sha256':sha(OLD/'pilot.py'),
              'outputs':{p.relative_to(HERE).as_posix():sha(p) for p in DEST.rglob('*') if p.is_file()},
              'locked_before_test':True}
        save(SUMMARY/'replication_lock.json',lock)
        manifest.update(status='FITS_COMPLETED',elapsed_seconds=time.time()-manifest['started_unix'])
        save(SUMMARY/'replication_manifest.json',manifest)
    except Exception as error:
        manifest.update(status='FAILED',error=str(error)); save(SUMMARY/'replication_manifest.json',manifest); raise


def score():
    assert_locked(); cfg=read(PROTOCOL); lock=read(SUMMARY/'replication_lock.json')
    assert lock['protocol_sha256']==sha(PROTOCOL) and lock['recipe_sha256']==sha(RECIPE)
    assert lock['source_sha256']==sha(HERE/'replicate.py') and lock['original_fit_source_sha256']==sha(OLD/'pilot.py')
    assert all(sha(HERE/name)==digest for name,digest in lock['outputs'].items())
    manifest=read(SUMMARY/'replication_manifest.json'); rows=[]
    assert len(manifest['runs'])==read(RECIPE)['expected_fits']
    for item in manifest['runs']:
        data=read(ART/item['condition']/'data.json')
        motifs=read_meme(DEST/item['condition']/f'{item["method"]}__a{item["algorithm_seed"]}'/'motifs.meme')
        for split in ['development','test']:
            metrics,_=evaluate(motifs,data,cfg,split)
            rows.extend({**item,**r} for r in metrics)
    frame=pd.DataFrame(rows); frame.to_csv(SUMMARY/'replication_metrics.tsv',sep='\t',index=False)
    primary=frame[(frame.split=='test')&(frame.tolerance_bp==2)]
    per_seed=primary.groupby(['scenario','data_seed','method']).agg(f1=('f1','mean'),
        f1_algorithm_min=('f1','min'),f1_algorithm_max=('f1','max'),success_fraction=('reference_eligible','mean')).reset_index()
    per_seed.to_csv(SUMMARY/'replication_per_data_seed.tsv',sep='\t',index=False)
    per_seed.groupby(['scenario','method']).agg(f1_mean=('f1','mean'),f1_data_min=('f1','min'),
        f1_data_max=('f1','max'),success_fraction=('success_fraction','mean')).reset_index().to_csv(SUMMARY/'replication_summary.tsv',sep='\t',index=False)
    manifest.update(status='COMPLETED'); save(SUMMARY/'replication_manifest.json',manifest)
    print('Original recipe replication scored, 90 actual fits',flush=True)


def verify():
    assert_locked(); cfg=read(PROTOCOL); lock=read(SUMMARY/'replication_lock.json')
    assert lock['protocol_sha256']==sha(PROTOCOL) and lock['recipe_sha256']==sha(RECIPE)
    assert lock['source_sha256']==sha(HERE/'replicate.py') and lock['original_fit_source_sha256']==sha(OLD/'pilot.py')
    assert all(sha(HERE/name)==digest for name,digest in lock['outputs'].items())
    manifest=read(SUMMARY/'replication_manifest.json')
    assert manifest['status']=='COMPLETED' and len(manifest['runs'])==read(RECIPE)['expected_fits']
    instances=0
    for item in manifest['runs']:
        data=read(ART/item['condition']/'data.json')
        path=DEST/item['condition']/f'{item["method"]}__a{item["algorithm_seed"]}'/'instances.tsv'
        for row in pd.read_csv(path,sep='\t').itertuples():
            assert data['discovery'][row.sequence_index][row.start0:row.end0]==row.fragment
            instances+=1
    result={'status':'PASSED','fits':len(manifest['runs']),'hashes_checked':len(lock['outputs']),
            'instances_checked':instances,'real_holdout_used':False}
    save(SUMMARY/'replication_verification.json',result);print(result,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['fit','score','verify'])
    {'fit':fit,'score':score,'verify':verify}[parser.parse_args().action]()
