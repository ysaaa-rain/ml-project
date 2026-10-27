"""Fixed, reference-free global window shifts and CNN direction ablation."""
import argparse
import copy
import hashlib
import importlib.metadata
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd
import torch
from sklearn.metrics import average_precision_score,roc_auc_score

from ..cnn_features_20261026.pipeline import (LocalCNN,background,candidate_cache,create_data,evaluate_view,
    make_pwm,maxima,oh_pool,rank_candidates,read,read_meme,reverse_complement,save,sequence_windows,
    setup,sha,tensor,train,write_meme)
from ..budget_validation_20261026.pipeline import sign_flip

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];PREVIOUS=HERE.parent/'cnn_features_20261026'
ART=HERE/'artifacts';SUMMARY=HERE/'summary';PROTOCOL=HERE/'protocol.json'
METHODS=['TIED_RAW','TIED_SHIFT','FWD_RAW','FWD_SHIFT','OH48_BG']
SHIFTS=[-2,-1,0,1,2]
COLUMNS=['motif','sequence_index','original_start0','start0','end0','strand','shift','fragment','activation']


class ForwardCNN(LocalCNN):
    def forward(self,x):
        plus,_=self.activations(x)
        return self.head(plus.amax(-1)).squeeze(-1)


def direction_maxima(model,sequences,direction):
    if direction=='TIED':return maxima(model,sequences)
    with torch.no_grad():plus,_=model.activations(tensor(sequences));values=plus.numpy()
    starts=values.argmax(-1);scores=np.take_along_axis(values,starts[...,None],axis=-1)[...,0]
    return scores,starts,np.ones(starts.shape,dtype=bool)


def shifted_fragment(sequence,start,forward,shift,width):
    new_start=int(start)+(shift if forward else -shift)
    if new_start<0 or new_start+width>len(sequence):return None
    fragment=sequence[new_start:new_start+width]
    return new_start,fragment if forward else reverse_complement(fragment)


def choose_shift(sequences,instances,bg,width):
    """Single Discovery-only shift per filter; identical eligible instances for every shift."""
    eligible=[r for r in instances if all(shifted_fragment(sequences[r['sequence_index']],r['start0'],r['forward'],s,width) is not None for s in SHIFTS)]
    if not eligible:return eligible,0,{str(s):0. for s in SHIFTS}
    values={}
    for shift in SHIFTS:
        fragments=[shifted_fragment(sequences[r['sequence_index']],r['start0'],r['forward'],shift,width)[1] for r in eligible]
        _,pwm=make_pwm(fragments);values[str(shift)]=float(np.sum(pwm*np.log2(pwm/bg)))
    best=min(SHIFTS,key=lambda s:(-values[str(s)],abs(s),s))
    return eligible,best,values


def save_pool(out,motifs,rows,positive,blank,cfg):
    write_meme(out/'motifs.meme',motifs,background(blank))
    frame=pd.DataFrame(rows,columns=COLUMNS);frame.to_csv(out/'instances.tsv',sep='\t',index=False)
    # Ranking is deliberately recomputed from the serialized PWM used by scoring.
    ranking=rank_candidates(read_meme(out/'motifs.meme'),frame.groupby('motif').size().to_dict(),positive,blank,cfg)
    save(out/'ranking.json',ranking)
    return ranking


def cnn_pools(model,positive,blank,cfg,direction,parent):
    score,start,forward=direction_maxima(model,positive,direction);null,_,_=direction_maxima(model,blank,direction)
    thresholds=np.quantile(null,.95,axis=0,method='higher');bg=background(blank)
    np.savez_compressed(parent/f'{direction}_maxima.npz',score=score,start=start,forward=forward,threshold=thresholds)
    motifs={'RAW':{},'SHIFT':{}};rows={'RAW':[],'SHIFT':[]};decisions=[]
    for j in range(score.shape[1]):
        selected=np.flatnonzero(score[:,j]>thresholds[j]);name=f'F{j:02d}'
        initial=[{'sequence_index':int(i),'start0':int(start[i,j]),'forward':bool(forward[i,j]),'activation':float(score[i,j])} for i in selected]
        eligible,best,values=choose_shift(positive,initial,bg,cfg['width'])
        retained=len(eligible)>=cfg['minimum_promoters']
        decisions.append({'filter':name,'selected':len(initial),'eligible':len(eligible),'excluded_edges':len(initial)-len(eligible),
            'retained':retained,'shift':best,'scores':values,'threshold':float(thresholds[j])})
        if not retained:continue
        for variant,shift in [('RAW',0),('SHIFT',best)]:
            for r in eligible:
                new_start,fragment=shifted_fragment(positive[r['sequence_index']],r['start0'],r['forward'],shift,cfg['width'])
                rows[variant].append({'motif':name,'sequence_index':r['sequence_index'],'original_start0':r['start0'],
                    'start0':new_start,'end0':new_start+cfg['width'],'strand':'+' if r['forward'] else '-',
                    'shift':shift,'fragment':fragment,'activation':r['activation']})
            _,motifs[variant][name]=make_pwm([r['fragment'] for r in rows[variant] if r['motif']==name])
    save(parent/f'{direction}_decisions.json',decisions)
    output={}
    for variant in ['RAW','SHIFT']:
        method=f'{direction}_{variant}';out=parent/method;out.mkdir()
        output[method]=save_pool(out,motifs[variant],rows[variant],positive,blank,cfg)
    return output


def views(rankings,cap=5):
    result=[]
    def add(method,view,names):
        result.append({'representation':method,'view':view,'candidates':names,'budget':len(names),
            'available':len(rankings[method]),'zero_budget':not names,'underfilled':len(names)<cap})
    for method in METHODS:
        names=[r['candidate'] for r in rankings[method]];add(method,'FULL',names);add(method,'CAP5',names[:cap])
    for view,pair in [('MATCH_SHIFT',['TIED_RAW','TIED_SHIFT']),('MATCH_OH',['TIED_SHIFT','OH48_BG'])]:
        m=min(cap,*(len(rankings[method]) for method in pair))
        for method in pair:add(method,view,[r['candidate'] for r in rankings[method][:m]])
    return result


def sources():
    names=list(read(PREVIOUS/'summary/manifest.json')['source_sha256'])
    paths=[ROOT/name for name in names]+[HERE/'pipeline.py',HERE/'__init__.py',PROTOCOL]
    return {p.relative_to(ROOT).as_posix():sha(p) for p in paths}


def protected_files():
    protected=read(PREVIOUS/'summary/protected_before.json')
    for name,digest in read(PREVIOUS/'summary/artifact_hashes.json').items():
        protected[(PREVIOUS/name).relative_to(ROOT).as_posix()]=digest
    for folder in [PREVIOUS,PREVIOUS/'summary']:
        for p in folder.iterdir():
            if p.is_file() and p.suffix in ['.py','.json','.md','.png','.tsv']:
                protected[p.relative_to(ROOT).as_posix()]=sha(p)
    assert all(sha(ROOT/name)==digest for name,digest in protected.items())
    return protected


def fit():
    cfg=read(PROTOCOL);setup(cfg);ART.mkdir(exist_ok=False);SUMMARY.mkdir(exist_ok=True)
    manifest={'status':'RUNNING','started_unix':time.time(),'source_sha256':sources(),
        'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'python':platform.python_version(),'dependencies':{k:importlib.metadata.version(k) for k in ['torch','numpy','pandas','scikit-learn']},
        'cnn_trains':0,'cluster_fits':0,'pwm_pools':0,'evaluation_sets':0,'real_holdout_used':False}
    save(SUMMARY/'manifest.json',manifest);save(SUMMARY/'protected_before.json',protected_files());curves=[]
    try:
        for scenario in cfg['scenarios']:
            for seed in cfg['data_seeds']:
                tag=f'{scenario["id"]}__s{seed}';folder=ART/tag;folder.mkdir();began=time.time()
                data=create_data(cfg,scenario,seed);save(folder/'data.json',data)
                for algorithm_seed in cfg['algorithm_seeds']:
                    parent=folder/f'a{algorithm_seed}';parent.mkdir();rankings={}
                    torch.manual_seed(algorithm_seed);initial=LocalCNN(cfg['cnn']['filters'],cfg['width'])
                    info={'initial_state_sha256':hashlib.sha256(b''.join(p.detach().numpy().tobytes() for p in initial.parameters())).hexdigest(),
                        'parameters':sum(p.numel() for p in initial.parameters())}
                    for direction in ['TIED','FWD']:
                        model=copy.deepcopy(initial) if direction=='TIED' else ForwardCNN(cfg['cnn']['filters'],cfg['width'])
                        if direction=='FWD':model.load_state_dict(initial.state_dict())
                        curve,label_hash=train(model,data['discovery']+data['discovery_control'],
                            [1]*len(data['discovery'])+[0]*len(data['discovery_control']),cfg['cnn'],algorithm_seed)
                        curves.extend({'condition':tag,'data_seed':seed,'algorithm_seed':algorithm_seed,'direction':direction,**r} for r in curve)
                        info[direction+'_labels_sha256']=label_hash;torch.save(model.state_dict(),parent/f'{direction}_model.pt')
                        rankings.update(cnn_pools(model,data['discovery'],data['discovery_control'],cfg,direction,parent));manifest['cnn_trains']+=1
                        manifest['pwm_pools']+=2
                    save(parent/'model_info.json',info)
                    out=parent/'OH48_BG';out.mkdir();oh_pool(data['discovery'],data['discovery_control'],cfg,algorithm_seed,out)
                    frame=pd.read_csv(out/'instances.tsv',sep='\t');rankings['OH48_BG']=rank_candidates(read_meme(out/'motifs.meme'),
                        frame.groupby('motif').size().to_dict(),data['discovery'],data['discovery_control'],cfg)
                    save(out/'ranking.json',rankings['OH48_BG']);manifest['cluster_fits']+=1;manifest['pwm_pools']+=1
                    selected=views(rankings,cfg['budget_cap']);save(parent/'views.json',selected);manifest['evaluation_sets']+=len(selected)
                    save(SUMMARY/'manifest.json',manifest)
                print(f'completed {tag}: {manifest["cnn_trains"]}/60 CNN, {manifest["pwm_pools"]}/150 pools ({time.time()-began:.1f}s)',flush=True)
        pd.DataFrame(curves).to_csv(SUMMARY/'training_curves.tsv',sep='\t',index=False)
        save(SUMMARY/'fit_lock.json',{'source_sha256':sources(),'before_development_test_scoring':True,
            'outputs':{p.relative_to(HERE).as_posix():sha(p) for p in ART.rglob('*') if p.is_file()}})
        manifest.update(status='FITS_COMPLETED',fit_seconds=time.time()-manifest['started_unix']);save(SUMMARY/'manifest.json',manifest)
    except Exception as e:
        manifest.update(status='FAILED',error=f'{type(e).__name__}: {e}');save(SUMMARY/'manifest.json',manifest);raise


def assert_locked():
    lock=read(SUMMARY/'fit_lock.json');manifest=read(SUMMARY/'manifest.json')
    assert lock['source_sha256']==sources()==manifest['source_sha256']
    assert all(sha(HERE/name)==digest for name,digest in lock['outputs'].items())


def score():
    cfg=read(PROTOCOL);setup(cfg);assert_locked();manifest=read(SUMMARY/'manifest.json');metrics=[];classifications=[]
    for folder in sorted(ART.iterdir()):
        data=read(folder/'data.json');prepared={k:sequence_windows(data[k],cfg['width']) for k in
            ['discovery_control','calibration_control','development','development_control','test','test_control']}
        for parent in sorted(folder.glob('a*')):
            algorithm_seed=int(parent.name[1:]);caches={}
            for method in METHODS:
                caches[method],_=candidate_cache(read_meme(parent/method/'motifs.meme'),data,cfg,prepared)
            for direction in ['TIED','FWD']:
                model=(LocalCNN if direction=='TIED' else ForwardCNN)(cfg['cnn']['filters'],cfg['width'])
                model.load_state_dict(torch.load(parent/f'{direction}_model.pt',weights_only=True));model.eval()
                for split in ['development','test']:
                    y=np.array([1]*len(data[split])+[0]*len(data[split+'_control']))
                    with torch.no_grad():pred=model(tensor(data[split]+data[split+'_control'])).numpy()
                    classifications.append({'scenario':data['scenario']['id'],'data_seed':data['seed'],'algorithm_seed':algorithm_seed,
                        'direction':direction,'split':split,'auroc':roc_auc_score(y,pred),'average_precision':average_precision_score(y,pred)})
            for view in read(parent/'views.json'):
                common={'scenario':data['scenario']['id'],'data_seed':data['seed'],'algorithm_seed':algorithm_seed,
                        **{k:v for k,v in view.items() if k!='candidates'}}
                for split in ['development','test']:
                    metrics.extend({**common,**r} for r in evaluate_view(caches[view['representation']],view['candidates'],data,cfg,split))
        print('scored '+folder.name,flush=True)
    frame=pd.DataFrame(metrics);frame.to_csv(SUMMARY/'metrics.tsv',sep='\t',index=False)
    pd.DataFrame(classifications).to_csv(SUMMARY/'classification.tsv',sep='\t',index=False)
    test=frame[(frame.split=='test')&(frame.tolerance_bp==2)]
    seeds=test.groupby(['scenario','data_seed','view','representation']).agg(f1=('f1','mean'),
        family_blank_fpr=('family_blank_fpr','mean'),budget=('budget','mean'),success_fraction=('reference_eligible','mean'),
        zero_budget_count=('zero_budget','sum'),underfill_count=('underfilled','sum')).reset_index()
    seeds.to_csv(SUMMARY/'per_data_seed.tsv',sep='\t',index=False)
    summary=seeds.groupby(['scenario','view','representation']).agg(f1_mean=('f1','mean'),f1_min=('f1','min'),f1_max=('f1','max'),
        family_blank_fpr=('family_blank_fpr','mean'),budget=('budget','mean'),success_fraction=('success_fraction','mean'),
        zero_budget_count=('zero_budget_count','sum'),underfill_count=('underfill_count','sum')).reset_index()
    summary.to_csv(SUMMARY/'test_summary.tsv',sep='\t',index=False)
    pairs=[]
    for view,comparator in [('MATCH_SHIFT','TIED_RAW'),('MATCH_OH','OH48_BG'),('CAP5','TIED_RAW'),('FULL','TIED_RAW')]:
        left=seeds[(seeds.view==view)&(seeds.representation=='TIED_SHIFT')]
        right=seeds[(seeds.view==view)&(seeds.representation==comparator)]
        merged=left.merge(right,on=['scenario','data_seed','view'],suffixes=('_shift','_baseline'),validate='one_to_one')
        for r in merged.to_dict('records'):
            pairs.append({'scenario':r['scenario'],'data_seed':r['data_seed'],'view':view,'comparator':comparator,
                'f1_shift':r['f1_shift'],'f1_baseline':r['f1_baseline'],'delta':r['f1_shift']-r['f1_baseline'],
                'family_fpr_shift':r['family_blank_fpr_shift'],'family_fpr_baseline':r['family_blank_fpr_baseline'],
                'budget_shift':r['budget_shift'],'budget_baseline':r['budget_baseline'],
                'zero_budget_count':r['zero_budget_count_shift']+r['zero_budget_count_baseline']})
    paired=pd.DataFrame(pairs);paired.to_csv(SUMMARY/'paired_deltas.tsv',sep='\t',index=False)
    chosen=paired[(paired.scenario=='mut30_random')&(paired.view=='MATCH_SHIFT')]
    cap=paired[(paired.scenario=='mut30_random')&(paired.view=='CAP5')]
    strong=paired[(paired.scenario=='strong_random')&(paired.view=='MATCH_SHIFT')]
    assert len(chosen)==len(cap)==len(strong)==len(cfg['data_seeds'])
    gate={'matched_mean_positive':bool(chosen.delta.mean()>0),'at_least4_of5_wins':int((chosen.delta>0).sum())>=4,
        'cap5_mean_positive':bool(cap.delta.mean()>0),'no_zero_budget':int(chosen.zero_budget_count.sum())==0,
        'family_fpr_not_above_raw':bool(chosen.family_fpr_shift.mean()<=chosen.family_fpr_baseline.mean()),
        'strong_mean_delta_at_least_minus_0_05':bool(strong.delta.mean()>=-.05)}
    save(SUMMARY/'decision.json',{'mean_delta':float(chosen.delta.mean()),'wins':int((chosen.delta>0).sum()),
        'losses':int((chosen.delta<0).sum()),'cap5_mean_delta':float(cap.delta.mean()),
        'sign_flip_p':sign_flip(chosen.delta.to_numpy()),'gate_checks':gate,'advance_gate_passed':all(gate.values()),
        'limits':'independent synthetic screening, five data seeds; no claim of significance, biological function or general model superiority'})
    manifest.update(status='COMPLETED',elapsed_seconds=time.time()-manifest['started_unix']);save(SUMMARY/'manifest.json',manifest)


def verify():
    cfg=read(PROTOCOL);setup(cfg);assert_locked();manifest=read(SUMMARY/'manifest.json')
    assert manifest['status']=='COMPLETED'
    for field in ['cnn_trains','cluster_fits','pwm_pools','evaluation_sets']:assert manifest[field]==cfg['expected_'+field]
    protected=read(SUMMARY/'protected_before.json');assert all(sha(ROOT/n)==d for n,d in protected.items())
    instances=pwms=rank_count=views_count=maps=truths=decisions_checked=0;edge_selected=edge_excluded=0
    for folder in sorted(ART.iterdir()):
        data=read(folder/'data.json');seen=set()
        for split in ['discovery','discovery_control','development','development_control','test','test_control','calibration_control']:
            canonical={min(s,reverse_complement(s)) for s in data[split]}
            assert len(canonical)==len(data[split]) and not canonical&seen;seen|=canonical
        for split in ['discovery','development','test']:
            for r in data[split+'_truth']:
                assert data[split][r['sequence_index']][r['start0']:r['end0']]==r['fragment'];truths+=1
        for parent in sorted(folder.glob('a*')):
            rankings={};bg=background(data['discovery_control'])
            for direction in ['TIED','FWD']:
                model=(LocalCNN if direction=='TIED' else ForwardCNN)(cfg['cnn']['filters'],cfg['width'])
                model.load_state_dict(torch.load(parent/f'{direction}_model.pt',weights_only=True));model.eval()
                score,start,forward=direction_maxima(model,data['discovery'],direction)
                null,_,_=direction_maxima(model,data['discovery_control'],direction);threshold=np.quantile(null,.95,axis=0,method='higher')
                stored=np.load(parent/f'{direction}_maxima.npz')
                for key,value in [('score',score),('start',start),('forward',forward),('threshold',threshold)]:assert np.array_equal(stored[key],value)
                maps+=1;recorded=read(parent/f'{direction}_decisions.json')
                for j,decision in enumerate(recorded):
                    selected=np.flatnonzero(score[:,j]>threshold[j])
                    initial=[{'sequence_index':int(i),'start0':int(start[i,j]),'forward':bool(forward[i,j]),'activation':float(score[i,j])} for i in selected]
                    eligible,best,values=choose_shift(data['discovery'],initial,bg,cfg['width'])
                    assert decision['selected']==len(initial) and decision['eligible']==len(eligible)
                    assert decision['shift']==best and decision['scores']==values and decision['retained']==(len(eligible)>=cfg['minimum_promoters'])
                    edge_selected+=len(initial);edge_excluded+=len(initial)-len(eligible);decisions_checked+=1
                    name=f'F{j:02d}'
                    for variant,shift in [('RAW',0),('SHIFT',best)]:
                        frame=pd.read_csv(parent/f'{direction}_{variant}'/'instances.tsv',sep='\t');rows=frame[frame.motif==name]
                        if len(eligible)<cfg['minimum_promoters']:assert len(rows)==0;continue
                        assert set(rows.sequence_index)=={r['sequence_index'] for r in eligible}
                        lookup={r['sequence_index']:r for r in eligible}
                        for r in rows.itertuples():
                            old=lookup[r.sequence_index];new_start,fragment=shifted_fragment(data['discovery'][r.sequence_index],old['start0'],old['forward'],shift,cfg['width'])
                            assert r.original_start0==old['start0'] and r.start0==new_start and r.fragment==fragment and r.shift==shift
            for method in METHODS:
                out=parent/method;motifs=read_meme(out/'motifs.meme');frame=pd.read_csv(out/'instances.tsv',sep='\t')
                assert set(frame.motif)==set(motifs) and not frame.duplicated(['motif','sequence_index']).any()
                for r in frame.itertuples():
                    dna=data['discovery'][r.sequence_index][r.start0:r.end0]
                    assert r.end0-r.start0==cfg['width'] and 0<=r.start0<=81-cfg['width']
                    assert r.fragment==(dna if r.strand=='+' else reverse_complement(dna));instances+=1
                for name,pwm in motifs.items():
                    rows=frame[frame.motif==name];assert len(rows)>=cfg['minimum_promoters']
                    _,rebuilt=make_pwm(rows.fragment.tolist());assert np.allclose(rebuilt,pwm,rtol=0,atol=1e-8);pwms+=1
                rankings[method]=rank_candidates(motifs,frame.groupby('motif').size().to_dict(),data['discovery'],data['discovery_control'],cfg)
                assert rankings[method]==read(out/'ranking.json');rank_count+=1
            actual=read(parent/'views.json');assert actual==views(rankings,cfg['budget_cap']);views_count+=len(actual)
    assert views_count==cfg['expected_evaluation_sets']
    inventory={p.relative_to(HERE).as_posix():sha(p) for p in ART.rglob('*') if p.is_file()}
    inventory.update({p.relative_to(HERE).as_posix():sha(p) for p in SUMMARY.iterdir() if p.is_file() and p.name not in ['artifact_hashes.json','verification.json']})
    save(SUMMARY/'artifact_hashes.json',inventory)
    result={'status':'PASSED','artifact_hashes':len(inventory),'historical_files_unchanged':len(protected),
        'instances_checked':instances,'pwms_rebuilt':pwms,'rankings_recomputed':rank_count,'cnn_activation_maps_recomputed':maps,
        'filter_shift_decisions_recomputed':decisions_checked,'selected_instances_before_edge':edge_selected,'edge_excluded':edge_excluded,
        'truth_sites_checked':truths,'evaluation_sets':views_count,'exact_RC_duplicates_per_condition':0,'real_holdout_used':False,
        'limits':'saved-model and artifact consistency; no independent retraining or near-homology audit'}
    save(SUMMARY/'verification.json',result);print(result,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['fit','score','verify'])
    {'fit':fit,'score':score,'verify':verify}[parser.parse_args().action]()
