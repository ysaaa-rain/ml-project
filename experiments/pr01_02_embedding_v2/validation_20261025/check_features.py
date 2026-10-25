"""Numerical QA against prior crop encoder; no discovery or test scoring."""
import numpy as np

from ..bottleneck_20261025.pilot import local_features as original_local_features
from .pipeline import ART,PROTOCOL,SUMMARY,assert_locked,load_encoder,read,save,sha


def main():
    assert_locked();cfg=read(PROTOCOL);before=sha(SUMMARY/'encoder.json')
    model,tokenizer=load_encoder(cfg);rows=[]
    for scenario,tier in [(s,'main') for s in cfg['scenarios']]+[(s,'stress') for s in cfg['stress_scenarios']]:
        seed=cfg['data_seeds'][0] if tier=='main' else cfg['stress_seeds'][0]
        folder=ART/f'{scenario["id"]}__s{seed}';data=read(folder/'data.json')
        _,reference=original_local_features(data['discovery'][:1],model,tokenizer,{**cfg,'batch_size':64})
        cached=np.load(folder/'local_discovery.npy',mmap_mode='r')[:1]
        max_error=float(np.abs(reference-cached).max())
        matches=bool(np.allclose(reference,cached,rtol=1e-4,atol=1e-5))
        rows.append({'condition':folder.name,'checked_sequences':1,'checked_core_windows':70,
                     'reference_batch':64,'current_batch':cfg['batch_size'],'max_absolute_error':max_error,
                     'allclose_rtol1e_4_atol1e_5':matches})
        if not matches:raise ValueError(f'Crop cache differs from original encoder: {folder.name}')
    assert sha(SUMMARY/'encoder.json')==before
    result={'status':'PASSED','checks':rows,'limits':'six Discovery sequences only; numerical equivalence, not bitwise or exhaustive replay'}
    save(SUMMARY/'feature_equivalence.json',result);print(result,flush=True)


if __name__=='__main__':main()
