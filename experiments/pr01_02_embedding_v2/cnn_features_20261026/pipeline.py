"""Sequence-supervised local CNN features; no planted-site access during fitting."""
import argparse
import copy
import hashlib
import importlib.metadata
import platform
import subprocess
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import torch
from sklearn.metrics import accuracy_score,average_precision_score,roc_auc_score

from ..budget_validation_20261026.pipeline import (ALPHABET,background,candidate_cache,cluster,create_data,
    evaluate_view,make_pwm,onehot,rank_candidates,raw_instances,read,read_meme,reject_background,
    save,select_instances,sequence_windows,sha,write_meme)
from ..bottleneck_20261025.pilot import reverse_complement

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
PREVIOUS=HERE.parent/'budget_validation_20261026'
ART=HERE/'artifacts';SUMMARY=HERE/'summary';PROTOCOL=HERE/'protocol.json'
METHODS=['CNN_LABEL','CNN_SHUFFLED','CNN_RANDOM','OH48_BG']
INSTANCE_COLUMNS=['motif','sequence_index','start0','end0','strand','fragment','selection_score']


class LocalCNN(torch.nn.Module):
    def __init__(self,filters,width):
        super().__init__()
        self.conv=torch.nn.Conv1d(4,filters,width)
        self.head=torch.nn.Linear(filters,1)

    def activations(self,x):
        return torch.relu(self.conv(x)),torch.relu(self.conv(x.flip((1,2))))

    def forward(self,x):
        plus,minus=self.activations(x)
        pooled=torch.maximum(plus.amax(-1),minus.amax(-1))
        return self.head(pooled).squeeze(-1)


def tensor(sequences):
    codes=np.array([[ALPHABET.index(c) for c in s] for s in sequences])
    return torch.from_numpy(np.eye(4,dtype=np.float32)[codes].transpose(0,2,1).copy())


def maxima(model,sequences):
    """Receptive-field coordinates in parent input; minus fragments are RC-normalized."""
    with torch.no_grad():
        plus,minus=model.activations(tensor(sequences))
        values=torch.cat((plus,minus),-1).numpy()
    nw=plus.shape[-1];index=values.argmax(-1)
    score=np.take_along_axis(values,index[...,None],axis=-1)[...,0]
    forward=index<nw
    start=np.where(forward,index,nw-1-(index-nw))
    return score,start,forward


def train(model,sequences,labels,cfg,seed,shuffled=False):
    """Only Discovery sequence labels. No truth, reference, Development or Test argument."""
    rng=np.random.default_rng(seed);x=tensor(sequences);y=np.array(labels,dtype=np.float32)
    if shuffled:y=y[rng.permutation(len(y))]
    label_hash=hashlib.sha256(y.tobytes()).hexdigest();y=torch.from_numpy(y)
    optimizer=torch.optim.Adam(model.parameters(),lr=cfg['learning_rate'],weight_decay=cfg['weight_decay'])
    rows=[];model.train()
    for epoch in range(cfg['epochs']):
        order=rng.permutation(len(y));loss_sum=0.
        for left in range(0,len(y),cfg['batch_size']):
            ids=order[left:left+cfg['batch_size']]
            optimizer.zero_grad(set_to_none=True)
            loss=torch.nn.functional.binary_cross_entropy_with_logits(model(x[ids]),y[ids])
            loss.backward();optimizer.step();loss_sum+=loss.item()*len(ids)
        if epoch==0 or (epoch+1)%10==0:rows.append({'epoch':epoch+1,'loss':loss_sum/len(y)})
    model.eval()
    return rows,label_hash


def cnn_pool(model,sequences,blank,cfg,out):
    score,start,forward=maxima(model,sequences);null,_,_=maxima(model,blank)
    thresholds=np.quantile(null,.95,axis=0,method='higher');rows=[];motifs={}
    for j in range(score.shape[1]):
        selected=np.flatnonzero(score[:,j]>thresholds[j])
        if len(selected)<cfg['minimum_promoters']:continue
        name=f'F{j:02d}'
        for i in selected:
            begin=int(start[i,j]);fragment=sequences[i][begin:begin+cfg['width']]
            rows.append({'motif':name,'sequence_index':int(i),'start0':begin,'end0':begin+cfg['width'],
                'strand':'+' if forward[i,j] else '-',
                'fragment':fragment if forward[i,j] else reverse_complement(fragment),
                'selection_score':float(score[i,j])})
        _,motifs[name]=make_pwm([r['fragment'] for r in rows if r['motif']==name])
    np.savez_compressed(out/'activation_maxima.npz',score=score,start=start,forward=forward,threshold=thresholds)
    save_pool(out,motifs,rows,sequences,blank,cfg)


def oh_pool(sequences,blank,cfg,seed,out):
    features=onehot(sequences,cfg['width']).reshape(-1,48)
    km,labels,distance=cluster(features,seed,cfg);joblib.dump(km,out/'kmeans.joblib')
    np.savez_compressed(out/'assignments.npz',labels=labels,distance=distance)
    selections,_=select_instances(labels,distance,82-cfg['width'],1.0)
    bg=background(blank);blank_windows=sequence_windows(blank,cfg['width']);motifs={};rows=[]
    for label,selected in selections.items():
        raw=raw_instances(sequences,selected,distance,cfg['width'])
        if len(raw)<cfg['minimum_promoters']:continue
        _,pwm=make_pwm([r['fragment'] for r in raw]);kept,_=reject_background(raw,pwm,blank_windows,cfg,bg)
        if len(kept)<cfg['minimum_promoters']:continue
        name=f'C{label:02d}';_,motifs[name]=make_pwm([r['fragment'] for r in kept])
        rows.extend({'motif':name,'sequence_index':r['sequence_index'],'start0':r['start0'],'end0':r['end0'],
                     'strand':'+','fragment':r['fragment'],'selection_score':r['distance']} for r in kept)
    save_pool(out,motifs,rows,sequences,blank,cfg)


def save_pool(out,motifs,rows,sequences,blank,cfg):
    write_meme(out/'motifs.meme',motifs,background(blank))
    frame=pd.DataFrame(rows,columns=INSTANCE_COLUMNS);frame.to_csv(out/'instances.tsv',sep='\t',index=False)
    support=frame.groupby('motif').size().to_dict()
    save(out/'ranking.json',rank_candidates(motifs,support,sequences,blank,cfg))


def views(rankings,cap=5):
    result=[]
    def add(method,view,names):
        result.append({'representation':method,'view':view,'candidates':names,'budget':len(names),
            'available':len(rankings[method]),'underfilled':len(names)<cap,'zero_budget':not names})
    for method,ranking in rankings.items():
        names=[r['candidate'] for r in ranking];add(method,'FULL',names);add(method,'CAP5',names[:cap])
    m=min(cap,len(rankings['CNN_LABEL']),len(rankings['OH48_BG']))
    for method in ['CNN_LABEL','OH48_BG']:add(method,'MATCH_OH48',[r['candidate'] for r in rankings[method][:m]])
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


def setup(cfg):
    torch.set_num_threads(cfg['threads']);torch.use_deterministic_algorithms(True)


def fit():
    cfg=read(PROTOCOL);setup(cfg);ART.mkdir(exist_ok=False);SUMMARY.mkdir(exist_ok=True)
    manifest={'status':'RUNNING','started_unix':time.time(),'source_sha256':sources(),
        'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'python':platform.python_version(),'dependencies':{k:importlib.metadata.version(k) for k in ['torch','numpy','pandas','scikit-learn']},
        'cnn_trains':0,'random_models':0,'cluster_fits':0,'pwm_pools':0,'evaluation_sets':0,'real_holdout_used':False}
    save(SUMMARY/'manifest.json',manifest);save(SUMMARY/'protected_before.json',protected_files());training=[]
    try:
        for scenario in cfg['scenarios']:
            for seed in cfg['data_seeds']:
                tag=f'{scenario["id"]}__s{seed}';folder=ART/tag;folder.mkdir();began=time.time()
                data=create_data(cfg,scenario,seed);save(folder/'data.json',data)
                for algorithm_seed in cfg['algorithm_seeds']:
                    parent=folder/f'a{algorithm_seed}';parent.mkdir();rankings={}
                    torch.manual_seed(algorithm_seed);initial=LocalCNN(cfg['cnn']['filters'],cfg['width'])
                    save(parent/'model_info.json',{'parameters':sum(p.numel() for p in initial.parameters()),
                        'initial_state_sha256':hashlib.sha256(b''.join(p.detach().numpy().tobytes() for p in initial.parameters())).hexdigest(),
                        'supervision':'sequence-level only; no coordinate labels; no external pretrained weights'})
                    for method in METHODS:
                        out=parent/method;out.mkdir()
                        if method.startswith('CNN'):
                            model=copy.deepcopy(initial)
                            if method!='CNN_RANDOM':
                                curve,label_hash=train(model,data['discovery']+data['discovery_control'],
                                    [1]*len(data['discovery'])+[0]*len(data['discovery_control']),cfg['cnn'],algorithm_seed,method=='CNN_SHUFFLED')
                                training.extend({'condition':tag,'method':method,'algorithm_seed':algorithm_seed,**r} for r in curve)
                                save(out/'training.json',{'labels_sha256':label_hash,'shuffled':method=='CNN_SHUFFLED',
                                    'epochs':cfg['cnn']['epochs'],'checkpoint':'fixed_final_epoch'})
                                manifest['cnn_trains']+=1
                            else:manifest['random_models']+=1
                            torch.save(model.state_dict(),out/'model.pt');cnn_pool(model,data['discovery'],data['discovery_control'],cfg,out)
                        else:
                            oh_pool(data['discovery'],data['discovery_control'],cfg,algorithm_seed,out);manifest['cluster_fits']+=1
                        rankings[method]=read(out/'ranking.json');manifest['pwm_pools']+=1
                    selected=views(rankings,cfg['budget_cap']);save(parent/'views.json',selected);manifest['evaluation_sets']+=len(selected)
                    save(SUMMARY/'manifest.json',manifest)
                print(f'completed {tag}: {manifest["pwm_pools"]}/72 pools ({time.time()-began:.1f}s)',flush=True)
        pd.DataFrame(training).to_csv(SUMMARY/'training_curves.tsv',sep='\t',index=False)
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
                out=parent/method;motifs=read_meme(out/'motifs.meme');caches[method],_=candidate_cache(motifs,data,cfg,prepared)
                if method.startswith('CNN'):
                    model=LocalCNN(cfg['cnn']['filters'],cfg['width']);model.load_state_dict(torch.load(out/'model.pt',weights_only=True));model.eval()
                    for split in ['development','test']:
                        y=np.array([1]*len(data[split])+[0]*len(data[split+'_control']))
                        with torch.no_grad():pred=model(tensor(data[split]+data[split+'_control'])).numpy()
                        classifications.append({'scenario':data['scenario']['id'],'data_seed':data['seed'],'algorithm_seed':algorithm_seed,
                            'method':method,'split':split,'auroc':roc_auc_score(y,pred),'average_precision':average_precision_score(y,pred),
                            'accuracy_logit0':accuracy_score(y,pred>0)})
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
    primary=seeds[(seeds.scenario=='mut30_random')&(seeds.view=='MATCH_OH48')]
    cnn=primary[primary.representation=='CNN_LABEL'];oh=primary[primary.representation=='OH48_BG']
    paired=cnn.merge(oh,on=['scenario','data_seed','view'],suffixes=('_cnn','_oh'),validate='one_to_one')
    paired['delta']=paired.f1_cnn-paired.f1_oh;paired.to_csv(SUMMARY/'primary_pairs.tsv',sep='\t',index=False)
    cap=summary[(summary.scenario=='mut30_random')&(summary.view=='CAP5')].set_index('representation')
    checks={'matched_mean_positive':bool(paired.delta.mean()>0),'all3_data_seeds_positive':bool((paired.delta>0).all()),
        'cap5_exceeds_controls_and_oh':bool(all(cap.loc['CNN_LABEL','f1_mean']>cap.loc[m,'f1_mean'] for m in METHODS[1:])),
        'no_zero_common_budget':bool((paired.zero_budget_count_cnn==0).all()),
        'family_fpr_not_above_oh':bool(paired.family_blank_fpr_cnn.mean()<=paired.family_blank_fpr_oh.mean())}
    save(SUMMARY/'pilot_decision.json',{'mean_delta':float(paired.delta.mean()),'wins':int((paired.delta>0).sum()),
        'data_seeds':len(paired),'checks':checks,'advance_gate_passed':all(checks.values()),
        'limits':'screening only; supervised vs unsupervised pipelines; not significance, standalone embedding advantage or functional motif proof'})
    manifest.update(status='COMPLETED',elapsed_seconds=time.time()-manifest['started_unix']);save(SUMMARY/'manifest.json',manifest)


def verify():
    cfg=read(PROTOCOL);setup(cfg);assert_locked();manifest=read(SUMMARY/'manifest.json')
    assert manifest['status']=='COMPLETED'
    for field in ['cnn_trains','random_models','cluster_fits','pwm_pools','evaluation_sets']:
        assert manifest[field]==cfg['expected_'+field]
    protected=read(SUMMARY/'protected_before.json');assert all(sha(ROOT/n)==d for n,d in protected.items())
    instances=0;pwms=0;rankings_checked=0;views_checked=0;truths=0;cnn_maps=0
    for folder in sorted(ART.iterdir()):
        data=read(folder/'data.json');seen=set()
        for split in ['discovery','discovery_control','development','development_control','test','test_control','calibration_control']:
            canonical={min(s,reverse_complement(s)) for s in data[split]}
            assert len(canonical)==len(data[split]) and not canonical&seen;seen|=canonical
        for split in ['discovery','development','test']:
            for r in data[split+'_truth']:
                assert data[split][r['sequence_index']][r['start0']:r['end0']]==r['fragment'];truths+=1
        for parent in sorted(folder.glob('a*')):
            ranking={}
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
                ranking[method]=rank_candidates(motifs,frame.groupby('motif').size().to_dict(),data['discovery'],data['discovery_control'],cfg)
                assert ranking[method]==read(out/'ranking.json');rankings_checked+=1
                if method.startswith('CNN'):
                    model=LocalCNN(cfg['cnn']['filters'],cfg['width']);model.load_state_dict(torch.load(out/'model.pt',weights_only=True));model.eval()
                    score,start,forward=maxima(model,data['discovery']);null,_,_=maxima(model,data['discovery_control'])
                    saved=np.load(out/'activation_maxima.npz');threshold=np.quantile(null,.95,axis=0,method='higher')
                    assert np.array_equal(score,saved['score']) and np.array_equal(start,saved['start']) and np.array_equal(forward,saved['forward'])
                    assert np.array_equal(threshold,saved['threshold'])
                    expected={f'F{j:02d}' for j in range(score.shape[1]) if (score[:,j]>threshold[j]).sum()>=cfg['minimum_promoters']}
                    assert expected==set(motifs)
                    for name in expected:
                        j=int(name[1:]);selected=frame[frame.motif==name]
                        assert set(selected.sequence_index)==set(np.flatnonzero(score[:,j]>threshold[j]))
                        for r in selected.itertuples():
                            assert r.start0==start[r.sequence_index,j] and (r.strand=='+')==forward[r.sequence_index,j]
                    cnn_maps+=1
                else:assert joblib.load(out/'kmeans.joblib').n_features_in_==48
            actual=read(parent/'views.json');assert actual==views(ranking,cfg['budget_cap']);views_checked+=len(actual)
    assert views_checked==cfg['expected_evaluation_sets']
    inventory={p.relative_to(HERE).as_posix():sha(p) for p in ART.rglob('*') if p.is_file()}
    inventory.update({p.relative_to(HERE).as_posix():sha(p) for p in SUMMARY.iterdir() if p.is_file() and p.name not in ['artifact_hashes.json','verification.json']})
    save(SUMMARY/'artifact_hashes.json',inventory)
    result={'status':'PASSED','artifact_hashes':len(inventory),'historical_files_unchanged':len(protected),
        'instances_rebuilt':instances,'pwms_rebuilt':pwms,'rankings_recomputed':rankings_checked,
        'cnn_activation_maps_recomputed':cnn_maps,'truth_sites_checked':truths,'evaluation_sets':views_checked,
        'exact_RC_duplicates_per_condition':0,'real_holdout_used':False,
        'limits':'saved-model activation recomputation and artifact consistency, not independent retraining or near-homology audit'}
    save(SUMMARY/'verification.json',result);print(result,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('action',choices=['fit','score','verify'])
    {'fit':fit,'score':score,'verify':verify}[parser.parse_args().action]()
