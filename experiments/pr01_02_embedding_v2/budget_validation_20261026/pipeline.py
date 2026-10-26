"""Locked Discovery candidate ranking and fresh Test confirmation at equal budgets."""
import argparse
import gc
import importlib.metadata
import importlib.util
import itertools
import platform
import subprocess
import sys
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd

from ..refinement_20261026.pipeline import (ALPHABET,background,cluster,create_data,download,local_order,
    make_pwm,nested_pca,onehot,raw_instances,read,read_meme,reject_background,save,scan_windows,
    select_instances,sequence_windows,sha,similarity,write_meme)

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PREVIOUS=HERE.parent/'refinement_20261026'
VALIDATION=HERE.parent/'validation_20261025'
PILOT=HERE.parent/'bottleneck_20261025'
ART=HERE/'artifacts';SUMMARY=HERE/'summary';PROTOCOL=HERE/'protocol.json'
COLUMNS=['motif','sequence_index','original_start0','start0','end0','strand','shift','fragment','distance']


def sources():
    paths=[PROTOCOL,HERE/'pipeline.py',HERE/'__init__.py',PREVIOUS/'pipeline.py',VALIDATION/'pipeline.py',
           PILOT/'pilot.py',HERE.parent/'encoder.py',HERE.parent/'discovery.py']
    return {p.relative_to(ROOT).as_posix():sha(p) for p in paths}


def load_encoder(cfg):
    import torch
    torch.set_num_threads(cfg['threads']);info=read(PILOT/'model_info.json')['models']['mini']
    folder=PILOT/'model_cache/mini'
    if not all((folder/name).is_file() for name in info['files']):
        folder=HERE/'model_cache'
        for name,digest in info['files'].items():
            repository='neuralbioinfo/nbrg-transformers' if name.endswith('.py') else info['repository']
            revision=info['code_revision'] if name.endswith('.py') else info['revision']
            download(f'https://huggingface.co/{repository}/resolve/{revision}/{name}',folder/name,digest)
    for name,digest in info['files'].items():
        if sha(folder/name)!=digest:raise ValueError('Pinned encoder hash mismatch')
    modules={}
    for name in ['models','tokenizer']:
        spec=importlib.util.spec_from_file_location('emb_budget_'+name,folder/(name+'.py'))
        module=importlib.util.module_from_spec(spec);sys.modules[spec.name]=module;spec.loader.exec_module(module);modules[name]=module
    tokenizer=modules['tokenizer'].LCATokenizer(vocab_file=str(folder/'vocab.txt'),kmer=6,shift=1)
    config=modules['models'].ProkBertConfig.from_pretrained(str(folder),local_files_only=True)
    model,loading=modules['models'].ProkBertModel.from_pretrained(str(folder),config=config,local_files_only=True,output_loading_info=True)
    if any(loading.get(k) for k in ['missing_keys','unexpected_keys','mismatched_keys','error_msgs']):raise ValueError(loading)
    model.eval().requires_grad_(False);save(SUMMARY/'encoder.json',{**info,'loading':loading,'device':'cpu'})
    return model,tokenizer


def build_pools(sequences,blank,labels,distance,cfg,parent):
    selections,_=select_instances(labels,distance,82-cfg['width'],1.0)
    bg=background(blank);blank_windows=sequence_windows(blank,cfg['width'])
    motifs={v:{} for v in cfg['variants']};instances={v:[] for v in cfg['variants']};stats=[]
    for label,selected in selections.items():
        rows=raw_instances(sequences,selected,distance,cfg['width'])
        if len(rows)<cfg['minimum_promoters']:continue
        name=f'C{label:02d}';_,pwm=make_pwm([r['fragment'] for r in rows])
        kept,threshold=reject_background(rows,pwm,blank_windows,cfg,bg)
        for variant,chosen in [('RAW',rows),('BG_ONLY',kept)]:
            retained=len(chosen)>=cfg['minimum_promoters']
            stats.append({'candidate':name,'variant':variant,'before':len(rows),'after':len(chosen),
                          'retained':retained,'threshold':threshold if variant=='BG_ONLY' else None})
            if retained:
                _,motifs[variant][name]=make_pwm([r['fragment'] for r in chosen])
                instances[variant].extend({'motif':name,**r} for r in chosen)
    for variant in cfg['variants']:
        out=parent/variant;out.mkdir();write_meme(out/'motifs.meme',motifs[variant],bg)
        pd.DataFrame(instances[variant],columns=COLUMNS).to_csv(out/'instances.tsv',sep='\t',index=False)
    pd.DataFrame(stats).to_csv(parent/'rejection.tsv',sep='\t',index=False)


def rank_candidates(motifs,support,sequences,blank,cfg):
    """Discovery-only ranking. Does not accept planted truth, reference, or Test."""
    bg=background(blank);positive=sequence_windows(sequences,cfg['width']);negative=sequence_windows(blank,cfg['width'])
    rows=[]
    for name,pwm in motifs.items():
        null,_,_=scan_windows(negative,pwm,bg);threshold=float(np.quantile(null,.95,method='higher'))
        scores,_,_=scan_windows(positive,pwm,bg)
        pr=float((scores>threshold).mean());nr=float((null>threshold).mean())
        rows.append({'candidate':name,'score':pr-nr,'positive_presence':pr,'blank_presence':nr,
            'threshold':threshold,'information_bits':float((pwm*np.log2(pwm/bg)).sum()),'support':int(support[name])})
    return sorted(rows,key=lambda r:(-r['score'],-r['information_bits'],-r['support'],r['candidate']))


def memberships(rankings,cap):
    views=[]
    for representation,rows in rankings.items():
        names=[r['candidate'] for r in rows]
        for view,selected in [('FULL',names),('CAP5',names[:cap])]:
            views.append({'view':view,'representation':representation,'candidates':selected,'available':len(names),
                          'budget':len(selected),'cap':cap,'underfilled':len(names)<cap,'zero_budget':not selected})
    for comparator in ['OH24','OH48']:
        count=min(cap,len(rankings['EMB128']),len(rankings[comparator]))
        for representation in ['EMB128',comparator]:
            views.append({'view':'MATCH_'+comparator,'representation':representation,
                'candidates':[r['candidate'] for r in rankings[representation][:count]],
                'available':len(rankings[representation]),'budget':count,'cap':cap,
                'underfilled':count<cap,'zero_budget':count==0})
    return views


def protected_files():
    protected=read(PREVIOUS/'summary/protected_before.json')
    for name,digest in read(PREVIOUS/'summary/artifact_hashes.json').items():
        protected[(PREVIOUS/name).relative_to(ROOT).as_posix()]=digest
    for folder in [PREVIOUS,PREVIOUS/'summary']:
        for p in folder.iterdir():
            if p.is_file() and p.suffix in ['.py','.json','.md','.png','.tsv']:protected[p.relative_to(ROOT).as_posix()]=sha(p)
    if not all(sha(ROOT/name)==digest for name,digest in protected.items()):raise ValueError('Historical files changed')
    return protected


def fit():
    cfg=read(PROTOCOL);ART.mkdir(exist_ok=False);SUMMARY.mkdir(exist_ok=True)
    manifest={'status':'RUNNING','started_unix':time.time(),'protocol_sha256':sha(PROTOCOL),'source_sha256':sources(),
        'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'python':platform.python_version(),
        'dependencies':{k:importlib.metadata.version(k) for k in ['torch','transformers','numpy','pandas','scipy','scikit-learn']},
        'real_holdout_used':False,'cluster_fits':0,'pwm_pools':0,'evaluation_sets':0,'conditions':[]}
    save(SUMMARY/'manifest.json',manifest);save(SUMMARY/'protected_before.json',protected_files())
    try:
        model,tokenizer=load_encoder(cfg)
        for scenario in cfg['scenarios']:
            for seed in cfg['data_seeds']:
                tag=f'{scenario["id"]}__s{seed}';folder=ART/tag;folder.mkdir();began=time.time()
                data=create_data(cfg,scenario,seed);save(folder/'data.json',data)
                features=local_order(data['discovery'],model,tokenizer,cfg,folder/'local_discovery.npy')
                pca,emb=nested_pca(features,128,cfg,folder/'embedding_pca.joblib')
                oh=onehot(data['discovery'],cfg['width']);oh_pca,oh24=nested_pca(oh,24,cfg,folder/'onehot_pca.joblib')
                reduced={'EMB128':emb,'OH24':oh24,'OH48':oh.reshape(-1,48)}
                for algorithm_seed in cfg['algorithm_seeds']:
                    rankings={v:{} for v in cfg['variants']}
                    for representation in cfg['representations']:
                        name=representation['id'];parent=folder/f'{name}__a{algorithm_seed}';parent.mkdir()
                        km,labels,distance=cluster(reduced[name],algorithm_seed,cfg);manifest['cluster_fits']+=1
                        joblib.dump(km,parent/'kmeans.joblib');np.savez_compressed(parent/'assignments.npz',labels=labels,distance=distance)
                        build_pools(data['discovery'],data['discovery_control'],labels,distance,cfg,parent)
                        for variant in cfg['variants']:
                            out=parent/variant;motifs=read_meme(out/'motifs.meme');frame=pd.read_csv(out/'instances.tsv',sep='\t')
                            support=frame.groupby('motif').size().to_dict()
                            ranking=rank_candidates(motifs,support,data['discovery'],data['discovery_control'],cfg)
                            save(out/'ranking.json',ranking);rankings[variant][name]=ranking;manifest['pwm_pools']+=1
                    for variant in cfg['variants']:
                        views=memberships(rankings[variant],cfg['budget_cap'])
                        save(folder/f'views__a{algorithm_seed}__{variant}.json',views);manifest['evaluation_sets']+=len(views)
                    save(SUMMARY/'manifest.json',manifest)
                manifest['conditions'].append(tag);save(SUMMARY/'manifest.json',manifest)
                print(f'completed {tag}: {manifest["cluster_fits"]}/135 fits, {manifest["evaluation_sets"]}/900 views ({time.time()-began:.1f}s)',flush=True)
                del features,pca,emb,oh,oh24,oh_pca,reduced;gc.collect()
        save(SUMMARY/'fit_lock.json',{'protocol_sha256':sha(PROTOCOL),'source_sha256':sources(),'before_test_scoring':True,
            'outputs':{p.relative_to(HERE).as_posix():sha(p) for p in ART.rglob('*') if p.is_file()}})
        manifest.update(status='FITS_COMPLETED',fit_elapsed_seconds=time.time()-manifest['started_unix']);save(SUMMARY/'manifest.json',manifest)
    except Exception as error:
        manifest.update(status='FAILED',error=f'{type(error).__name__}: {error}');save(SUMMARY/'manifest.json',manifest);raise


def assert_locked():
    lock=read(SUMMARY/'fit_lock.json');manifest=read(SUMMARY/'manifest.json')
    assert lock['protocol_sha256']==sha(PROTOCOL)==manifest['protocol_sha256']
    assert lock['source_sha256']==sources()==manifest['source_sha256']
    assert all(sha(HERE/name)==digest for name,digest in lock['outputs'].items())


def candidate_cache(motifs,data,cfg,prepared):
    """Test scans cached once per PWM, then reused by pre-locked budget memberships."""
    bg=background(data['discovery_control']);width=cfg['width'];mutation=data['scenario']['mutation']
    reference=np.full((width,4),mutation/3)
    for i,char in enumerate(cfg['synthetic_pattern']):reference[i,ALPHABET.index(char)]=1-mutation
    cache={};nulls=[]
    for name,pwm in motifs.items():
        match=similarity(pwm,reference,minimum_overlap=width)
        calibration,_,_=scan_windows(prepared['calibration_control'],pwm,bg)
        threshold=float(np.quantile(calibration,cfg['evaluation']['calibration_quantile'],method='higher'))
        record={'cosine':match[0],'orientation':match[2],'threshold':threshold}
        for split in ['development','test']:
            scores,starts,forward=scan_windows(prepared[split],pwm,bg)
            blank,_,_=scan_windows(prepared[split+'_control'],pwm,bg)
            record[split]={'present':scores>threshold,'start':starts,'forward':forward,'blank':blank>threshold}
            nulls.append({'split':split,'candidate':name,'threshold':threshold,'raw_blank_fpr':float((blank>threshold).mean())})
        cache[name]=record
    return cache,nulls


def evaluate_view(cache,names,data,cfg,split):
    family=np.zeros(len(data[split+'_control']),dtype=bool)
    for name in names:family|=cache[name][split]['blank']
    common={'split':split,'candidate_count':len(names),'family_blank_fpr':float(family.mean())}
    winner=max(names,key=lambda n:(cache[n]['cosine'],n)) if names else None
    if winner is None:
        return [{**common,'candidate':'','pwm_cosine':0.,'reference_eligible':False,'tolerance_bp':t,
                 'f1':0.,'precision':0.,'recall':0.,'tp':0,'prediction_count':0,'truth_n':len(data[split])}
                for t in cfg['evaluation']['position_tolerances']]
    record=cache[winner];scans=record[split];eligible=record['cosine']>=cfg['evaluation']['reference_cosine_gate']
    actual=np.array([r['start0'] for r in data[split+'_truth']]);actual_forward=np.array([r['strand']=='+' for r in data[split+'_truth']])
    canonical=scans['forward'] if record['orientation']=='+' else ~scans['forward']
    predictions=int(scans['present'].sum()+scans['blank'].sum()) if eligible else 0
    rows=[]
    for tolerance in cfg['evaluation']['position_tolerances']:
        tp=int((scans['present']&(canonical==actual_forward)&(np.abs(scans['start']-actual)<=tolerance)).sum()) if eligible else 0
        precision=tp/predictions if predictions else 0.;recall=tp/len(actual)
        rows.append({**common,'candidate':winner,'pwm_cosine':record['cosine'],'reference_eligible':eligible,
            'tolerance_bp':tolerance,'f1':2*precision*recall/(precision+recall) if precision+recall else 0.,
            'precision':precision,'recall':recall,'tp':tp,'prediction_count':predictions,'truth_n':len(actual)})
    return rows


def sign_flip(values):
    values=np.asarray(values,dtype=float);observed=abs(values.mean())
    means=np.array([(values*np.array(signs)).mean() for signs in itertools.product([-1,1],repeat=len(values))])
    return float((np.abs(means)>=observed-1e-12).mean())


def score():
    assert_locked();cfg=read(PROTOCOL);manifest=read(SUMMARY/'manifest.json');metrics=[];nulls=[]
    for tag in manifest['conditions']:
        folder=ART/tag;data=read(folder/'data.json')
        prepared={s:sequence_windows(data[s],cfg['width']) for s in ['calibration_control','development','development_control','test','test_control']}
        for algorithm_seed in cfg['algorithm_seeds']:
            for variant in cfg['variants']:
                caches={}
                for representation in cfg['representations']:
                    name=representation['id'];pool=folder/f'{name}__a{algorithm_seed}'/variant
                    caches[name],records=candidate_cache(read_meme(pool/'motifs.meme'),data,cfg,prepared)
                    nulls.extend({'condition':tag,'scenario':data['scenario']['id'],'data_seed':data['seed'],
                        'algorithm_seed':algorithm_seed,'representation':name,'variant':variant,**r} for r in records)
                for view in read(folder/f'views__a{algorithm_seed}__{variant}.json'):
                    common={'condition':tag,'scenario':data['scenario']['id'],'data_seed':data['seed'],
                        'algorithm_seed':algorithm_seed,'variant':variant,**{k:v for k,v in view.items() if k!='candidates'}}
                    for split in ['development','test']:
                        metrics.extend({**common,**r} for r in evaluate_view(caches[view['representation']],view['candidates'],data,cfg,split))
        print('scored '+tag,flush=True)
    frame=pd.DataFrame(metrics);frame.to_csv(SUMMARY/'metrics.tsv',sep='\t',index=False)
    pd.DataFrame(nulls).to_csv(SUMMARY/'candidate_null_fpr.tsv',sep='\t',index=False)
    primary=frame[(frame.split=='test')&(frame.tolerance_bp==2)]
    seeds=primary.groupby(['scenario','data_seed','variant','view','representation']).agg(
        f1=('f1','mean'),algorithm_min=('f1','min'),algorithm_max=('f1','max'),success_fraction=('reference_eligible','mean'),
        family_blank_fpr=('family_blank_fpr','mean'),budget=('budget','mean'),available=('available','mean'),
        zero_budget_count=('zero_budget','sum'),underfill_count=('underfilled','sum')).reset_index()
    seeds.to_csv(SUMMARY/'per_data_seed.tsv',sep='\t',index=False)
    seeds.groupby(['scenario','variant','view','representation']).agg(f1_mean=('f1','mean'),f1_data_min=('f1','min'),
        f1_data_max=('f1','max'),data_seeds=('data_seed','nunique'),success_fraction=('success_fraction','mean'),
        family_blank_fpr=('family_blank_fpr','mean'),budget=('budget','mean'),zero_budget_count=('zero_budget_count','sum'),
        underfill_count=('underfill_count','sum')).reset_index().to_csv(SUMMARY/'test_summary.tsv',sep='\t',index=False)
    pairs=[]
    for view in ['FULL','CAP5','MATCH_OH24','MATCH_OH48']:
        for comparator in ['OH24','OH48']:
            if view.startswith('MATCH_') and view!='MATCH_'+comparator:continue
            left=seeds[(seeds.view==view)&(seeds.representation=='EMB128')]
            right=seeds[(seeds.view==view)&(seeds.representation==comparator)]
            merged=left.merge(right,on=['scenario','data_seed','variant','view'],suffixes=('_emb','_baseline'),validate='one_to_one')
            for r in merged.to_dict('records'):
                pairs.append({'scenario':r['scenario'],'data_seed':r['data_seed'],'variant':r['variant'],'view':view,
                    'comparator':comparator,'f1_emb':r['f1_emb'],'f1_baseline':r['f1_baseline'],'delta':r['f1_emb']-r['f1_baseline'],
                    'family_fpr_emb':r['family_blank_fpr_emb'],'family_fpr_baseline':r['family_blank_fpr_baseline'],
                    'budget_emb':r['budget_emb'],'budget_baseline':r['budget_baseline'],
                    'zero_budget_count':r['zero_budget_count_emb']+r['zero_budget_count_baseline']})
    paired=pd.DataFrame(pairs);paired.to_csv(SUMMARY/'paired_deltas.tsv',sep='\t',index=False)
    paired.groupby(['scenario','variant','view','comparator']).agg(mean_delta=('delta','mean'),
        wins=('delta',lambda x:int((x>0).sum())),losses=('delta',lambda x:int((x<0).sum())),
        min_delta=('delta','min'),max_delta=('delta','max'),family_fpr_emb=('family_fpr_emb','mean'),
        family_fpr_baseline=('family_fpr_baseline','mean')).reset_index().to_csv(SUMMARY/'paired_summary.tsv',sep='\t',index=False)
    chosen=paired[(paired.scenario=='mut30_random')&(paired.variant=='BG_ONLY')&(paired.view=='MATCH_OH48')&(paired.comparator=='OH48')]
    capped=paired[(paired.scenario=='mut30_random')&(paired.variant=='BG_ONLY')&(paired.view=='CAP5')&(paired.comparator=='OH48')]
    assert len(chosen)==len(capped)==len(cfg['data_seeds'])
    gate={'matched_mean_positive':bool(chosen.delta.mean()>0),'cap5_mean_positive':bool(capped.delta.mean()>0),
          'at_least4_of5_wins':int((chosen.delta>0).sum())>=4,'no_zero_common_budget':int(chosen.zero_budget_count.sum())==0,
          'family_fpr_not_above_OH48':bool(chosen.family_fpr_emb.mean()<=chosen.family_fpr_baseline.mean())}
    decision={'endpoint':cfg['primary_endpoint'],'data_seeds':len(chosen),'mean_delta':float(chosen.delta.mean()),
        'wins':int((chosen.delta>0).sum()),'cap5_mean_delta':float(capped.delta.mean()),
        'paired_sign_flip_p':sign_flip(chosen.delta.to_numpy()),'gate_checks':gate,'advance_gate_passed':all(gate.values()),
        'interpretation':'gate only prioritizes further study, never functional/general superiority proof; n5 two-sided p minimum .0625'}
    save(SUMMARY/'primary_decision.json',decision)
    manifest.update(status='COMPLETED',elapsed_seconds=time.time()-manifest['started_unix']);save(SUMMARY/'manifest.json',manifest)


def verify():
    assert_locked();cfg=read(PROTOCOL);manifest=read(SUMMARY/'manifest.json')
    assert manifest['status']=='COMPLETED'
    assert manifest['cluster_fits']==cfg['expected_cluster_fits'] and manifest['pwm_pools']==cfg['expected_pwm_pools']
    assert manifest['evaluation_sets']==cfg['expected_evaluation_sets']
    protected=read(SUMMARY/'protected_before.json');assert all(sha(ROOT/n)==v for n,v in protected.items())
    instances=0;pwms=0;subset_pairs=0;views_checked=0
    for folder in ART.iterdir():
        data=read(folder/'data.json')
        for algorithm_seed in cfg['algorithm_seeds']:
            for representation in cfg['representations']:
                parent=folder/f'{representation["id"]}__a{algorithm_seed}';frames={}
                for variant in cfg['variants']:
                    out=parent/variant;frame=pd.read_csv(out/'instances.tsv',sep='\t');frames[variant]=frame
                    motifs=read_meme(out/'motifs.meme');assert set(frame.motif)==set(motifs)
                    assert not frame.duplicated(['motif','sequence_index']).any()
                    for r in frame.itertuples():
                        assert r.strand=='+' and r.shift==0 and r.start0==r.original_start0
                        assert 0<=r.start0<r.end0<=81 and r.end0-r.start0==12
                        assert data['discovery'][r.sequence_index][r.start0:r.end0]==r.fragment;instances+=1
                    for name,pwm in motifs.items():
                        fragments=frame[frame.motif==name].fragment.tolist();assert len(fragments)>=cfg['minimum_promoters']
                        _,rebuilt=make_pwm(fragments);assert np.allclose(rebuilt,pwm,rtol=0,atol=1e-8);pwms+=1
                    ranking=read(out/'ranking.json');assert {r['candidate'] for r in ranking}==set(motifs)
                    assert ranking==sorted(ranking,key=lambda r:(-r['score'],-r['information_bits'],-r['support'],r['candidate']))
                keys=lambda f:{(r.motif,r.sequence_index,r.start0,r.fragment) for r in f.itertuples()}
                assert keys(frames['BG_ONLY'])<=keys(frames['RAW']);subset_pairs+=1
            for variant in cfg['variants']:
                rankings={r['id']:read(folder/f'{r["id"]}__a{algorithm_seed}'/variant/'ranking.json') for r in cfg['representations']}
                expected=memberships(rankings,cfg['budget_cap']);actual=read(folder/f'views__a{algorithm_seed}__{variant}.json')
                assert actual==expected;views_checked+=len(actual)
                for comparator in ['OH24','OH48']:
                    pair=[v for v in actual if v['view']=='MATCH_'+comparator]
                    assert len(pair)==2 and pair[0]['budget']==pair[1]['budget']
    assert views_checked==cfg['expected_evaluation_sets']
    inventory={p.relative_to(HERE).as_posix():sha(p) for p in ART.rglob('*') if p.is_file()}
    inventory.update({p.relative_to(HERE).as_posix():sha(p) for p in SUMMARY.iterdir() if p.is_file() and p.name not in ['verification.json','artifact_hashes.json']})
    save(SUMMARY/'artifact_hashes.json',inventory)
    result={'status':'PASSED','hashes_checked':len(inventory),'protected_files_unchanged':len(protected),
        'cluster_fits':manifest['cluster_fits'],'pwm_pools':manifest['pwm_pools'],'evaluation_sets':views_checked,
        'instances_checked':instances,'pwms_rebuilt':pwms,'background_subset_pairs':subset_pairs,'real_holdout_used':False,
        'limits':['five synthetic data seeds; algorithm seeds not extra independent samples','artifact consistency, not independent model rerun']}
    save(SUMMARY/'verification.json',result);print(result,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['fit','score','verify'])
    {'fit':fit,'score':score,'verify':verify}[parser.parse_args().action]()
