"""Paired DNA fragment refinement; parent sequences and historical outputs are read-only."""
import argparse
import gc
import importlib.metadata
import importlib.util
import platform
import subprocess
import sys
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from ..validation_20261025.pipeline import (ALPHABET, background, cluster, create_data,
    download, local_order, make_pwm, nested_pca, onehot, read, read_meme,
    reverse_complement, save, select_instances, sha, similarity, write_meme)

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PREVIOUS=HERE.parent/'validation_20261025'
PILOT=HERE.parent/'bottleneck_20261025'
PROTOCOL=HERE/'protocol.json'
ART=HERE/'artifacts'
SUMMARY=HERE/'summary'
LOOKUP=np.zeros(256,dtype=np.int8)
LOOKUP[np.frombuffer(b'ACGT',dtype=np.uint8)]=np.arange(4)
INSTANCE_COLUMNS=['motif','sequence_index','original_start0','start0','end0','strand','shift','fragment','distance']


def sources():
    paths=[PROTOCOL,HERE/'pipeline.py',HERE/'__init__.py',PREVIOUS/'pipeline.py',PILOT/'pilot.py',
           HERE.parent/'encoder.py',HERE.parent/'discovery.py']
    return {p.relative_to(ROOT).as_posix():sha(p) for p in paths}


def load_encoder(cfg):
    import torch
    torch.set_num_threads(cfg['threads'])
    info=read(PILOT/'model_info.json')['models']['mini']
    folder=PILOT/'model_cache/mini'
    if not all((folder/name).is_file() for name in info['files']):
        folder=HERE/'model_cache'
        for name,digest in info['files'].items():
            repository='neuralbioinfo/nbrg-transformers' if name.endswith('.py') else info['repository']
            revision=info['code_revision'] if name.endswith('.py') else info['revision']
            download(f'https://huggingface.co/{repository}/resolve/{revision}/{name}',folder/name,digest)
    for name,digest in info['files'].items():
        if sha(folder/name)!=digest: raise ValueError('Encoder cache differs from pinned source')
    modules={}
    for name in ['models','tokenizer']:
        spec=importlib.util.spec_from_file_location('emb_refinement_'+name,folder/(name+'.py'))
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module
        spec.loader.exec_module(module);modules[name]=module
    tokenizer=modules['tokenizer'].LCATokenizer(vocab_file=str(folder/'vocab.txt'),kmer=6,shift=1)
    config=modules['models'].ProkBertConfig.from_pretrained(str(folder),local_files_only=True)
    model,loading=modules['models'].ProkBertModel.from_pretrained(str(folder),config=config,
        local_files_only=True,output_loading_info=True)
    if any(loading.get(k) for k in ['missing_keys','unexpected_keys','mismatched_keys','error_msgs']):
        raise ValueError(f'Incomplete model: {loading}')
    model.eval().requires_grad_(False)
    save(SUMMARY/'encoder.json',{**info,'loading':loading,'device':'cpu'})
    return model,tokenizer


def encode_dna(strings):
    return LOOKUP[np.frombuffer(''.join(strings).encode('ascii'),dtype=np.uint8)].reshape(len(strings),-1)


def sequence_windows(sequences,width):
    return np.lib.stride_tricks.sliding_window_view(encode_dna(sequences),width,axis=1)


def scan_windows(windows,pwm,bg):
    width=len(pwm)
    forward=np.log2(pwm/bg)[np.arange(width),windows].sum(axis=-1)
    reverse=np.log2(pwm[::-1,::-1]/bg)[np.arange(width),windows].sum(axis=-1)
    scores=np.concatenate([forward,reverse],axis=1);indices=scores.argmax(axis=1)
    return scores[np.arange(len(windows)),indices],indices%forward.shape[1],indices<forward.shape[1]


def raw_instances(sequences,selected,distance,width):
    return [{'sequence_index':row,'original_start0':start,'start0':start,'end0':start+width,
             'strand':'+','shift':0,'fragment':sequences[row][start:start+width],
             'distance':float(distance[row*(82-width)+start])} for row,start in selected.items()]


def align_instances(sequences,instances,cfg,bg,allow_shift):
    """Reference-free three-round refinement of actual DNA; offsets stay relative to original."""
    if not instances: return []
    width=cfg['width']; options=[]; option_meta=[]; mask=[]
    shifts=cfg['alignment']['shift_order'] if allow_shift else [0]
    for instance in instances:
        row=[];meta=[];valid=[]
        sequence=sequences[instance['sequence_index']]
        for shift in shifts:
            start=instance['original_start0']+shift
            legal=0<=start<=len(sequence)-width
            fragment=sequence[start:start+width] if legal else 'A'*width
            for strand in cfg['alignment']['strand_order']:
                row.append(fragment if strand=='+' else reverse_complement(fragment))
                meta.append((start,strand,shift));valid.append(legal)
        options.append(row);option_meta.append(meta);mask.append(valid)
    shape=(len(instances),len(options[0]),width)
    encoded=encode_dna([s for row in options for s in row]).reshape(shape)
    mask=np.array(mask,dtype=bool)
    seed=min(instances,key=lambda r:(r['distance'],r['sequence_index'],r['original_start0']))
    _,pwm=make_pwm([seed['fragment']]*cfg['alignment']['initial_seed_copies'])
    chosen=None
    for _ in range(cfg['alignment']['iterations']):
        scores=np.log2(pwm/bg)[np.arange(width),encoded].sum(axis=-1)
        scores[~mask]=-np.inf
        chosen=scores.argmax(axis=1)
        fragments=[options[i][j] for i,j in enumerate(chosen)]
        _,pwm=make_pwm(fragments)
    result=[]
    for i,j in enumerate(chosen):
        start,strand,shift=option_meta[i][j]
        result.append({**instances[i],'start0':start,'end0':start+width,'strand':strand,
                       'shift':shift,'fragment':options[i][j]})
    return result


def reject_background(instances,pwm,blank_windows,cfg,bg):
    scores,_,_=scan_windows(blank_windows,pwm,bg)
    threshold=float(np.quantile(scores,cfg['background_rejection']['quantile'],method='higher'))
    if not instances:return [],threshold
    encoded=encode_dna([r['fragment'] for r in instances])
    values=np.log2(pwm/bg)[np.arange(cfg['width']),encoded].sum(axis=-1)
    return [r for r,value in zip(instances,values) if value>threshold],threshold


def refine_clusters(sequences,blank_sequences,labels,distance,cfg,destination):
    """No ground truth input. BG is post-cluster rejection, not a new KMeans objective."""
    selections,_=select_instances(labels,distance,82-cfg['width'],1.0)
    bg=background(blank_sequences);blank_windows=sequence_windows(blank_sequences,cfg['width'])
    motifs={variant:{} for variant in cfg['variants']}
    instances={variant:[] for variant in cfg['variants']};statistics=[]
    for label,selected in selections.items():
        name=f'C{label:02d}';raw=raw_instances(sequences,selected,distance,cfg['width'])
        if len(raw)<cfg['minimum_promoters']:continue
        rc=align_instances(sequences,raw,cfg,bg,False)
        aligned=align_instances(sequences,raw,cfg,bg,True)
        _,raw_pwm=make_pwm([r['fragment'] for r in raw])
        _,aligned_pwm=make_pwm([r['fragment'] for r in aligned])
        raw_kept,raw_threshold=reject_background(raw,raw_pwm,blank_windows,cfg,bg)
        aligned_kept,aligned_threshold=reject_background(aligned,aligned_pwm,blank_windows,cfg,bg)
        variants=dict(zip(cfg['variants'],[raw,rc,aligned,raw_kept,aligned_kept]))
        for variant,rows in variants.items():
            retained=len(rows)>=cfg['minimum_promoters']
            threshold=raw_threshold if variant=='BG_ONLY' else aligned_threshold if variant=='RC_SHIFT_BG' else None
            statistics.append({'candidate':name,'variant':variant,'original_promoters':len(raw),
                'retained_promoters':len(rows),'candidate_retained':retained,'background_threshold':threshold,
                'reversed_instances':sum(r['strand']=='-' for r in rows),
                'shifted_instances':sum(r['shift']!=0 for r in rows)})
            if retained:
                _,motifs[variant][name]=make_pwm([r['fragment'] for r in rows])
                instances[variant].extend({'motif':name,**r} for r in rows)
    for variant in cfg['variants']:
        out=destination/variant;out.mkdir()
        write_meme(out/'motifs.meme',motifs[variant],bg)
        pd.DataFrame(instances[variant],columns=INSTANCE_COLUMNS).to_csv(out/'instances.tsv',sep='\t',index=False)
    pd.DataFrame(statistics).to_csv(destination/'refinement_stats.tsv',sep='\t',index=False)
    return {variant:len(motifs[variant]) for variant in cfg['variants']}


def protected_files():
    protected=read(PREVIOUS/'summary/protected_before.json')
    for name,digest in read(PREVIOUS/'summary/artifact_hashes.json').items():
        protected[(PREVIOUS/name).relative_to(ROOT).as_posix()]=digest
    for name,digest in read(PREVIOUS/'summary/replication_lock.json')['outputs'].items():
        protected[(PREVIOUS/name).relative_to(ROOT).as_posix()]=digest
    for folder in [PREVIOUS,PREVIOUS/'summary']:
        for path in folder.iterdir():
            if path.is_file() and path.suffix in ['.py','.md','.json','.png','.tsv']:
                protected[path.relative_to(ROOT).as_posix()]=sha(path)
    for name,digest in protected.items():
        if sha(ROOT/name)!=digest:raise ValueError(f'Prior evidence differs: {name}')
    return protected


def fit():
    cfg=read(PROTOCOL);ART.mkdir(exist_ok=False);SUMMARY.mkdir(exist_ok=True)
    manifest={'status':'RUNNING','started_unix':time.time(),'protocol_sha256':sha(PROTOCOL),'source_sha256':sources(),
        'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'python':platform.python_version(),'dependencies':{p:importlib.metadata.version(p) for p in ['torch','transformers','numpy','pandas','scipy','scikit-learn']},
        'real_holdout_used':False,'cluster_fits':0,'runs':[]}
    save(SUMMARY/'manifest.json',manifest);save(SUMMARY/'protected_before.json',protected_files())
    try:
        model,tokenizer=load_encoder(cfg)
        for scenario in cfg['scenarios']:
            for seed in cfg['data_seeds']:
                tag=f'{scenario["id"]}__s{seed}';folder=ART/tag;folder.mkdir();began=time.time()
                data=create_data(cfg,scenario,seed);save(folder/'data.json',data)
                features=local_order(data['discovery'],model,tokenizer,cfg,folder/'local_discovery.npy')
                pca,z=nested_pca(features,128,cfg,folder/'embedding_pca.joblib')
                oh=onehot(data['discovery'],cfg['width']);oh_pca,oh_z=nested_pca(oh,24,cfg,folder/'onehot_pca.joblib')
                for representation in cfg['representations']:
                    reduced=oh_z if representation['id']=='OH' else z
                    for algorithm_seed in cfg['algorithm_seeds']:
                        km,labels,distance=cluster(reduced,algorithm_seed,cfg);manifest['cluster_fits']+=1
                        out=folder/f'{representation["id"]}__a{algorithm_seed}';out.mkdir()
                        joblib.dump(km,out/'kmeans.joblib');np.savez_compressed(out/'assignments.npz',labels=labels,distance=distance)
                        counts=refine_clusters(data['discovery'],data['discovery_control'],labels,distance,cfg,out)
                        for variant,count in counts.items():
                            item={'condition':tag,'scenario':scenario['id'],'data_seed':seed,'algorithm_seed':algorithm_seed,
                                  'representation':representation['id'],'variant':variant,'candidate_count':count,'status':'COMPLETED'}
                            manifest['runs'].append(item)
                        save(SUMMARY/'manifest.json',manifest)
                print(f'completed {tag}: {len(manifest["runs"])}/240 outputs ({time.time()-began:.1f}s)',flush=True)
                del features,z,pca,oh,oh_z,oh_pca,reduced;gc.collect()
        save(SUMMARY/'fit_lock.json',{'protocol_sha256':sha(PROTOCOL),'source_sha256':sources(),'before_test_scoring':True,
            'outputs':{p.relative_to(HERE).as_posix():sha(p) for p in ART.rglob('*') if p.is_file()}})
        manifest.update(status='FITS_COMPLETED',fit_elapsed_seconds=time.time()-manifest['started_unix'])
        save(SUMMARY/'manifest.json',manifest)
    except Exception as error:
        manifest.update(status='FAILED',error=f'{type(error).__name__}: {error}');save(SUMMARY/'manifest.json',manifest);raise


def assert_locked():
    lock=read(SUMMARY/'fit_lock.json')
    if lock['protocol_sha256']!=sha(PROTOCOL) or lock['source_sha256']!=sources():raise ValueError('Protocol/source changed after fit')
    if not all(sha(HERE/name)==digest for name,digest in lock['outputs'].items()):raise ValueError('Fitted outputs changed')


def evaluate(motifs,data,cfg,split,prepared):
    width=cfg['width'];mutation=data['scenario']['mutation'];reference=np.full((width,4),mutation/3)
    for i,char in enumerate(cfg['synthetic_pattern']):reference[i,ALPHABET.index(char)]=1-mutation
    bg=background(data['discovery_control'])
    matches=[(similarity(pwm,reference,minimum_overlap=width),name) for name,pwm in motifs.items()]
    match,winner=max(matches,key=lambda x:(x[0][0],x[1])) if matches else ((0.,0,'+'),None)
    rows=[];nulls=[];family_null=np.zeros(len(data[split+'_control']),dtype=bool)
    for name,pwm in motifs.items():
        calibration,_,_=scan_windows(prepared['calibration_control'],pwm,bg)
        threshold=float(np.quantile(calibration,cfg['evaluation']['calibration_quantile'],method='higher'))
        scores,starts,forward=scan_windows(prepared[split],pwm,bg)
        blank,_,_=scan_windows(prepared[split+'_control'],pwm,bg)
        present=scores>threshold;false=blank>threshold;family_null|=false
        nulls.append({'split':split,'candidate':name,'raw_blank_fpr':float(false.mean()),
                      'raw_blank_hits':int(false.sum()),'blank_n':len(false),'threshold':threshold})
        if name!=winner:continue
        eligible=match[0]>=cfg['evaluation']['reference_cosine_gate']
        actual=np.array([r['start0'] for r in data[split+'_truth']])
        actual_forward=np.array([r['strand']=='+' for r in data[split+'_truth']])
        canonical=forward if match[2]=='+' else ~forward
        for tolerance in cfg['evaluation']['position_tolerances']:
            tp=int((present&(canonical==actual_forward)&(np.abs(starts-actual)<=tolerance)).sum()) if eligible else 0
            predictions=int(present.sum()+false.sum()) if eligible else 0
            precision=tp/predictions if predictions else 0.;recall=tp/len(actual)
            rows.append({'split':split,'tolerance_bp':tolerance,'candidate':winner,'candidate_count':len(motifs),
                'pwm_cosine':match[0],'reference_eligible':eligible,'precision':precision,'recall':recall,
                'f1':2*precision*recall/(precision+recall) if precision+recall else 0.,'tp':tp,
                'prediction_count':predictions,'truth_n':len(actual),'raw_blank_fpr':float(false.mean())})
    if winner is None:
        rows=[{'split':split,'tolerance_bp':t,'candidate':'','candidate_count':0,'pwm_cosine':0.,
               'reference_eligible':False,'precision':0.,'recall':0.,'f1':0.,'tp':0,'prediction_count':0,
               'truth_n':len(data[split]),'raw_blank_fpr':0.} for t in cfg['evaluation']['position_tolerances']]
    family={'split':split,'candidate_count':len(motifs),'raw_family_blank_fpr':float(family_null.mean()),
            'raw_family_blank_hits':int(family_null.sum()),'blank_n':len(family_null)}
    return rows,nulls,family


def score():
    assert_locked();cfg=read(PROTOCOL);manifest=read(SUMMARY/'manifest.json')
    metrics=[];nulls=[];families=[];refinement=[];data_cache={};prepared_cache={}
    for item in manifest['runs']:
        tag=item['condition']
        if tag not in data_cache:
            data=read(ART/tag/'data.json');data_cache[tag]=data
            prepared_cache[tag]={s:sequence_windows(data[s],cfg['width'])
                for s in ['development','development_control','test','test_control','calibration_control']}
        data=data_cache[tag];prepared=prepared_cache[tag]
        parent=ART/tag/f'{item["representation"]}__a{item["algorithm_seed"]}'
        motifs=read_meme(parent/item['variant']/'motifs.meme')
        for split in ['development','test']:
            rows,blank,family=evaluate(motifs,data,cfg,split,prepared)
            metrics.extend({**item,**row} for row in rows)
            nulls.extend({**item,**row} for row in blank);families.append({**item,**family})
    for scenario in cfg['scenarios']:
        for seed in cfg['data_seeds']:
            tag=f'{scenario["id"]}__s{seed}'
            for representation in cfg['representations']:
                for algorithm_seed in cfg['algorithm_seeds']:
                    table=pd.read_csv(ART/tag/f'{representation["id"]}__a{algorithm_seed}'/'refinement_stats.tsv',sep='\t')
                    common={'condition':tag,'scenario':scenario['id'],'data_seed':seed,
                            'representation':representation['id'],'algorithm_seed':algorithm_seed}
                    refinement.extend({**common,**row} for row in table.to_dict('records'))
    frame=pd.DataFrame(metrics);frame.to_csv(SUMMARY/'metrics.tsv',sep='\t',index=False)
    pd.DataFrame(nulls).to_csv(SUMMARY/'candidate_null_fpr.tsv',sep='\t',index=False)
    pd.DataFrame(families).to_csv(SUMMARY/'family_null_fpr.tsv',sep='\t',index=False)
    pd.DataFrame(refinement).to_csv(SUMMARY/'refinement_stats.tsv',sep='\t',index=False)
    primary=frame[(frame.split=='test')&(frame.tolerance_bp==2)]
    per_seed=primary.groupby(['scenario','data_seed','representation','variant']).agg(
        f1=('f1','mean'),algorithm_min=('f1','min'),algorithm_max=('f1','max'),
        success_fraction=('reference_eligible','mean'),candidate_count=('candidate_count','mean')).reset_index()
    per_seed.to_csv(SUMMARY/'per_data_seed.tsv',sep='\t',index=False)
    per_seed.groupby(['scenario','representation','variant']).agg(f1_mean=('f1','mean'),
        f1_data_min=('f1','min'),f1_data_max=('f1','max'),success_fraction=('success_fraction','mean'),
        candidate_count=('candidate_count','mean')).reset_index().to_csv(SUMMARY/'test_summary.tsv',sep='\t',index=False)
    pivot=per_seed.pivot(index=['scenario','data_seed','representation'],columns='variant',values='f1')
    pivot.subtract(pivot['RAW'],axis=0).reset_index().to_csv(SUMMARY/'paired_delta_vs_RAW.tsv',sep='\t',index=False)
    pivot=per_seed.pivot(index=['scenario','data_seed','variant'],columns='representation',values='f1')
    pivot['EMB_minus_OH']=pivot.EMB-pivot.OH
    pivot.reset_index().to_csv(SUMMARY/'paired_EMB_vs_OH.tsv',sep='\t',index=False)
    manifest.update(status='COMPLETED',elapsed_seconds=time.time()-manifest['started_unix'])
    save(SUMMARY/'manifest.json',manifest);print('All locked refinements scored',flush=True)


def verify():
    assert_locked();cfg=read(PROTOCOL);manifest=read(SUMMARY/'manifest.json')
    assert manifest['status']=='COMPLETED' and manifest['cluster_fits']==cfg['expected_cluster_fits']
    assert len(manifest['runs'])==cfg['expected_discovery_outputs']
    protected=read(SUMMARY/'protected_before.json')
    assert all(sha(ROOT/name)==digest for name,digest in protected.items())
    instances=0;truths=0;reversed_count=0;shifted_count=0
    for folder in ART.iterdir():
        data=read(folder/'data.json')
        for split in ['discovery','development','test']:
            for row in data[split+'_truth']:
                assert data[split][row['sequence_index']][row['start0']:row['end0']]==row['fragment'];truths+=1
        for path in folder.glob('*/*/instances.tsv'):
            table=pd.read_csv(path,sep='\t')
            for row in table.itertuples():
                fragment=data['discovery'][row.sequence_index][row.start0:row.end0]
                if row.strand=='-':fragment=reverse_complement(fragment);reversed_count+=1
                assert fragment==row.fragment and row.end0-row.start0==cfg['width']
                assert row.start0-row.original_start0==row.shift and abs(row.shift)<=2
                shifted_count+=int(row.shift!=0);instances+=1
    inventory={p.relative_to(HERE).as_posix():sha(p) for p in ART.rglob('*') if p.is_file()}
    inventory.update({p.relative_to(HERE).as_posix():sha(p) for p in SUMMARY.iterdir()
                      if p.is_file() and p.name not in ['artifact_hashes.json','verification.json']})
    save(SUMMARY/'artifact_hashes.json',inventory)
    result={'status':'PASSED','cluster_fits':manifest['cluster_fits'],'discovery_outputs':len(manifest['runs']),
            'hashes_checked':len(inventory),'protected_files_unchanged':len(protected),
            'instances_checked':instances,'reversed_instances_checked':reversed_count,
            'shifted_instances_checked':shifted_count,'truth_sites_checked':truths,'real_holdout_used':False,
            'limits':['three data seeds and two clustering seeds; bounded pilot',
                      'no independent rerun or full near-homology audit','background rejection is post-cluster, not new clustering']}
    save(SUMMARY/'verification.json',result);print(result,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['fit','score','verify'])
    {'fit':fit,'score':score,'verify':verify}[parser.parse_args().action]()
