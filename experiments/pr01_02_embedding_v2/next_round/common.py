import json
from pathlib import Path
import time
import importlib.metadata
import subprocess
import numpy as np
from ..audit import BRANCH,ROOT,DATA,sha,save

NEXT=BRANCH/'next_round'
PRESENTATION=BRANCH/'presentation'
SOURCES=['tjupan_bacillus_subtilis','tjupan_baumannii','tjupan_bradyrhizobium','tjupan_diphtheria','tjupan_escherichia_coli','tjupan_staphylococcus']
for directory in ['configs','stability_diagnostics','gc_diagnostics','pooling_ablation','synthetic_benchmark','real_data_comparison','figures','run_manifests']:
    (NEXT/directory).mkdir(parents=True,exist_ok=True)


def start(name,config,inputs):
    target=NEXT/'run_manifests'/f'{name}.json'
    if target.exists():
        raise FileExistsError(target)
    manifest={'status':'RUNNING','config':config,'started_unix':time.time(),'holdout_used':False,
      'device':'cpu','git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
      'input_sha256':{str(Path(p).relative_to(ROOT)):sha(p) for p in inputs},
      'code_sha256':{p.name:sha(p) for p in NEXT.glob('*.py')},
      'dependencies':{p:importlib.metadata.version(p) for p in ['torch','numpy','pandas','scipy','scikit-learn','transformers']}}
    save(target,manifest)
    return target,manifest


def finish(target,manifest,outputs,extra=None):
    import psutil
    memory=psutil.Process().memory_info()
    manifest.update(status='COMPLETED',elapsed_seconds=time.time()-manifest['started_unix'],
        process_peak_memory_bytes=getattr(memory,'peak_wset',memory.rss),
        memory_scope='process high-water mark; shared sequential experiments are not isolated peaks',
        output_sha256={str(Path(p).relative_to(ROOT)):sha(p) for p in outputs})
    if extra: manifest.update(extra)
    save(target,manifest)


def pool_weights(width,method,length=81,k=6):
    """Rows are normalized window weights over real token intervals."""
    starts=np.arange(length-k+1)
    positions=np.arange(length)
    covers=(positions[:,None]>=starts)&(positions[:,None]<starts+k)
    base_weights=covers/covers.sum(axis=1,keepdims=True)
    rows=[]
    for start in range(length-width+1):
        end=start+width
        overlap=np.maximum(0,np.minimum(starts+k,end)-np.maximum(starts,start))
        if method=='A0':
            weights=base_weights[start:end].mean(axis=0)
        elif method=='A1':
            weights=((starts>=start)&(starts+k<=end)).astype(float)
        elif method=='A2':
            weights=overlap/k
        elif method=='A3':
            weights=(overlap/k)*np.exp(-.5*((starts+(k-1)/2-(start+(width-1)/2))/(width/4))**2)
        else:
            raise ValueError(method)
        if weights.sum()==0: raise ValueError('No supported full token inside target window')
        rows.append(weights/weights.sum())
    return np.array(rows,dtype=np.float32)


def pool(token,width,method):
    return np.einsum('wt,ntd->nwd',pool_weights(width,method),token,optimize=True)
