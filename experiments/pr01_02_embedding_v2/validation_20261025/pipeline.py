"""New-seed validation; fixed nested PCA and paired instance selection."""
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
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from threadpoolctl import threadpool_limits

from ..bottleneck_20261025.pilot import (ALPHABET, background, crop_spec, download,
    encode, make_pwm, read, read_meme, reverse_complement, save, scan, sha, similarity, write_meme)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
OLD = HERE.parent / 'bottleneck_20261025'
PROTOCOL = HERE / 'protocol.json'
ART = HERE / 'artifacts'
SUMMARY = HERE / 'summary'


def sources():
    paths = [PROTOCOL, HERE / 'pipeline.py', HERE / '__init__.py', OLD / 'pilot.py',
             HERE.parent / 'encoder.py', HERE.parent / 'discovery.py']
    return {p.relative_to(ROOT).as_posix(): sha(p) for p in paths}


def load_encoder(cfg):
    """Read historical cache or download to this version; never write old evidence."""
    import torch
    torch.set_num_threads(cfg['threads'])
    info = read(OLD / 'model_info.json')['models']['mini']
    folder = OLD / 'model_cache/mini'
    if not all((folder / name).exists() for name in info['files']):
        folder = HERE / 'model_cache'
        for name, digest in info['files'].items():
            if name in ['models.py', 'tokenizer.py']:
                url = f'https://huggingface.co/neuralbioinfo/nbrg-transformers/resolve/{info["code_revision"]}/{name}'
            else:
                url = f'https://huggingface.co/{info["repository"]}/resolve/{info["revision"]}/{name}'
            download(url, folder / name, digest)
    for name, digest in info['files'].items():
        if sha(folder / name) != digest:
            raise ValueError('Encoder hash differs from frozen source')
    modules = {}
    for name in ['models', 'tokenizer']:
        spec = importlib.util.spec_from_file_location('emb_validation_' + name, folder / (name + '.py'))
        module = importlib.util.module_from_spec(spec); sys.modules[spec.name] = module
        spec.loader.exec_module(module); modules[name] = module
    tokenizer = modules['tokenizer'].LCATokenizer(vocab_file=str(folder / 'vocab.txt'), kmer=6, shift=1)
    config = modules['models'].ProkBertConfig.from_pretrained(str(folder), local_files_only=True)
    model, loading = modules['models'].ProkBertModel.from_pretrained(str(folder), config=config,
        local_files_only=True, output_loading_info=True)
    if any(loading.get(key) for key in ['missing_keys', 'unexpected_keys', 'mismatched_keys', 'error_msgs']):
        raise ValueError(f'Incomplete encoder: {loading}')
    model.eval().requires_grad_(False)
    save(SUMMARY / 'encoder.json', {**info, 'loading': loading,
        'parameters': sum(p.numel() for p in model.parameters()), 'device': 'cpu'})
    return model, tokenizer


def generate(n, scenario, stream, planted, cfg):
    rng = np.random.default_rng(stream); sequences, truth = [], []
    fraction = scenario['gc']; prob = [(1-fraction)/2, fraction/2, fraction/2, (1-fraction)/2]
    pattern = cfg['synthetic_pattern']; length = cfg['sequence_length']
    for index in range(n):
        bases = rng.choice(list(ALPHABET), length, p=prob)
        if scenario.get('persistence', 0):
            for position in range(1, length):
                if rng.random() < scenario['persistence']:
                    bases[position] = bases[position-1]
        if planted:
            start = int(rng.integers(0, length-len(pattern)+1)); fragment = list(pattern)
            for i, char in enumerate(fragment):
                if rng.random() < scenario['mutation']:
                    fragment[i] = rng.choice([c for c in ALPHABET if c != char])
            strand = '-' if rng.random() < scenario.get('reverse_probability', 0) else '+'
            fragment = ''.join(fragment)
            if strand == '-':
                fragment = reverse_complement(fragment)
            bases[start:start+len(pattern)] = list(fragment)
            truth.append({'sequence_index': index, 'start0': start, 'end0': start+len(pattern),
                          'strand': strand, 'fragment': fragment})
        sequences.append(''.join(bases))
    return sequences, truth


def create_data(cfg, scenario, seed):
    streams = np.random.SeedSequence(seed).spawn(7)
    data = {'scenario': scenario, 'seed': seed}
    for index, split in enumerate(['discovery', 'development', 'test']):
        data[split], data[split+'_truth'] = generate(cfg[split+'_n'], scenario, streams[index*2], True, cfg)
        data[split+'_control'], _ = generate(cfg[split+'_n'], scenario, streams[index*2+1], False, cfg)
    data['calibration_control'], _ = generate(cfg['calibration_n'], scenario, streams[6], False, cfg)
    seen = set()
    for split in ['discovery', 'discovery_control', 'development', 'development_control', 'test', 'test_control', 'calibration_control']:
        canonical = {min(s, reverse_complement(s)) for s in data[split]}
        if len(canonical) != len(data[split]) or seen & canonical:
            raise ValueError('Exact/RC duplicate within or across synthetic partitions')
        seen.update(canonical)
    return data


def local_order(sequences, model, tokenizer, cfg, target):
    """Encode each unique crop once; shared edge crops preserve exact core mapping."""
    width, context = cfg['width'], cfg['local_context']; nw = 82-width; nt = width-5
    result = np.lib.format.open_memmap(target, mode='w+', dtype='float32',
        shape=(len(sequences), nw, nt*model.config.hidden_size))
    mapping = {}
    for start in range(nw):
        left, offset = crop_spec(start, width, context)
        mapping.setdefault(left, []).append((start, offset))
    batch, slots = [], []; began = time.time()

    def flush():
        if not batch:
            return
        hidden = encode(batch, model, tokenizer, 6, cfg['batch_size'], cfg['layer'])
        for block, (row, left) in zip(hidden, slots):
            for start, offset in mapping[left]:
                result[row, start] = block[offset:offset+nt].reshape(-1)
        batch.clear(); slots.clear()

    for row, sequence in enumerate(sequences):
        for left in mapping:
            batch.append(sequence[left:left+context]); slots.append((row, left))
            if len(batch) == cfg['batch_size']:
                flush()
        if (row+1) % 64 == 0:
            print(f'local encoding {target.parent.name}/{target.name}: {row+1}/{len(sequences)} ({time.time()-began:.1f}s)', flush=True)
    flush(); result.flush()
    return result


def onehot(sequences, width):
    encoded = np.array([[ALPHABET.index(c) for c in s] for s in sequences])
    window = np.lib.stride_tricks.sliding_window_view(encoded, width, axis=1)
    return np.eye(4, dtype=np.float32)[window].reshape(len(sequences), 82-width, -1)


def nested_pca(features, components, cfg, path):
    with threadpool_limits(limits=cfg['threads']):
        pca = PCA(n_components=components, svd_solver='randomized', random_state=cfg['pca_seed'])
        reduced = pca.fit_transform(features.reshape(-1, features.shape[-1]))
    joblib.dump(pca, path)
    return pca, reduced


def cluster(reduced, algorithm_seed, cfg):
    with threadpool_limits(limits=cfg['threads']):
        km = MiniBatchKMeans(n_clusters=cfg['clusters'], n_init=cfg['n_init'], batch_size=2048,
                            random_state=algorithm_seed)
        labels = km.fit_predict(reduced)
        distance = np.linalg.norm(reduced-km.cluster_centers_[labels], axis=1)
    return km, labels, distance


def select_instances(labels, distances, nw, quantile):
    selections, statistics = {}, []
    for label in sorted(set(labels.tolist())):
        indices = np.flatnonzero(labels == label)
        threshold = float(np.quantile(distances[indices], quantile, method='linear'))
        eligible = indices[distances[indices] <= threshold]
        selected = {}
        for index in eligible[np.argsort(distances[eligible], kind='stable')]:
            selected.setdefault(int(index//nw), int(index%nw))
        selections[int(label)] = selected
        statistics.append({'cluster': int(label), 'distance_threshold': threshold,
            'windows_before': len(indices), 'windows_after': len(eligible), 'promoters_after': len(selected)})
    return selections, statistics


def build_pwm(sequences, labels, distances, quantile, cfg, destination):
    """No motif truth/reference input; paired filters reuse exact cluster membership."""
    width = cfg['width']; selections, statistics = select_instances(labels, distances, 82-width, quantile)
    motifs, instances = {}, []
    for label, selected in selections.items():
        if len(selected) < cfg['minimum_promoters']:
            continue
        name = f'C{label:02d}'
        _, motifs[name] = make_pwm([sequences[row][start:start+width] for row, start in selected.items()])
        instances.extend({'motif': name, 'sequence_index': row, 'start0': start, 'end0': start+width,
                          'strand': '+', 'fragment': sequences[row][start:start+width]}
                         for row, start in selected.items())
    pd.DataFrame(instances, columns=['motif','sequence_index','start0','end0','strand','fragment']).to_csv(destination/'instances.tsv', sep='\t', index=False)
    pd.DataFrame(statistics).to_csv(destination/'selection.tsv', sep='\t', index=False)
    return motifs


def protected_files():
    protected = read(OLD/'summary/protected_before.json')
    for name, digest in read(OLD/'summary/artifact_hashes.json').items():
        protected[(OLD/ name).relative_to(ROOT).as_posix()] = digest
    for path in OLD.glob('*'):
        if path.is_file() and path.suffix in ['.py','.json','.md']:
            protected[path.relative_to(ROOT).as_posix()] = sha(path)
    for name, digest in protected.items():
        if sha(ROOT/name) != digest:
            raise ValueError(f'Historical evidence already differs: {name}')
    return protected


def all_conditions(cfg):
    return [(s, seed, 'main') for s in cfg['scenarios'] for seed in cfg['data_seeds']] + [
        (s, seed, 'stress') for s in cfg['stress_scenarios'] for seed in cfg['stress_seeds']]


def fit_all():
    cfg = read(PROTOCOL); ART.mkdir(exist_ok=False); SUMMARY.mkdir(exist_ok=True)
    manifest = {'status':'RUNNING', 'started_unix':time.time(), 'protocol_sha256':sha(PROTOCOL),
        'source_sha256':sources(), 'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'python':platform.python_version(), 'dependencies':{k:importlib.metadata.version(k) for k in ['torch','transformers','numpy','pandas','scipy','scikit-learn']},
        'real_holdout_used':False, 'cluster_fits':0, 'runs':[]}
    save(SUMMARY/'manifest.json',manifest); save(SUMMARY/'protected_before.json',protected_files())
    try:
        model, tokenizer = load_encoder(cfg)
        for scenario, seed, tier in all_conditions(cfg):
            tag = f'{scenario["id"]}__s{seed}'; folder=ART/tag; folder.mkdir()
            began=time.time(); data=create_data(cfg,scenario,seed); data['tier']=tier; save(folder/'data.json',data)
            local=local_order(data['discovery'],model,tokenizer,cfg,folder/'local_discovery.npy')
            local_pca, local_z=nested_pca(local,128,cfg,folder/'local_pca.joblib')
            oh=onehot(data['discovery'],cfg['width']); oh_pca, oh_z=nested_pca(oh,24,cfg,folder/'onehot_pca.joblib')
            methods=cfg['methods'] if tier=='main' else [m for m in cfg['methods'] if m['id'] in cfg['stress_methods']]
            keys=sorted({(m['representation'],m['pca']) for m in methods})
            bg=background(data['discovery_control'])
            for representation, dimension in keys:
                reduced=oh_z if representation=='onehot' else local_z[:,:dimension]
                for algorithm_seed in cfg['algorithm_seeds']:
                    km, labels, distance=cluster(reduced,algorithm_seed,cfg); manifest['cluster_fits']+=1
                    assignment=f'{representation}_P{dimension}__a{algorithm_seed}'
                    joblib.dump(km,folder/(assignment+'.joblib'))
                    np.savez_compressed(folder/(assignment+'.npz'),labels=labels,distance=distance)
                    for method in [m for m in methods if m['representation']==representation and m['pca']==dimension]:
                        out=folder/f'{method["id"]}__a{algorithm_seed}'; out.mkdir()
                        motifs=build_pwm(data['discovery'],labels,distance,method['distance_quantile'],cfg,out)
                        write_meme(out/'motifs.meme',motifs,bg)
                        item={'condition':tag,'scenario':scenario['id'],'tier':tier,'data_seed':seed,
                              'algorithm_seed':algorithm_seed,'method':method['id'],'candidate_count':len(motifs),
                              'assignment':assignment,'pca_dimensions':dimension,'status':'COMPLETED'}
                        manifest['runs'].append(item); save(out/'run.json',item)
                    save(SUMMARY/'manifest.json',manifest)
            print(f'condition completed {tag}: {len(manifest["runs"])}/{cfg["expected_discovery_outputs"]} outputs ({time.time()-began:.1f}s)',flush=True)
            del local, local_z, local_pca, oh, oh_z, oh_pca, reduced; gc.collect()
        lock={'protocol_sha256':sha(PROTOCOL),'source_sha256':sources(),'locked_before_test':True,
              'outputs':{p.relative_to(HERE).as_posix():sha(p) for p in ART.rglob('*')
                         if p.is_file() and p.name not in ['local_discovery.npy']}}
        save(SUMMARY/'fit_lock.json',lock)
        manifest.update(status='FITS_COMPLETED',fit_elapsed_seconds=time.time()-manifest['started_unix'])
        save(SUMMARY/'manifest.json',manifest)
    except Exception as error:
        manifest.update(status='FAILED',error=f'{type(error).__name__}: {error}')
        save(SUMMARY/'manifest.json',manifest); raise


def assert_locked():
    lock=read(SUMMARY/'fit_lock.json')
    if lock['protocol_sha256']!=sha(PROTOCOL) or lock['source_sha256']!=sources():
        raise ValueError('Frozen protocol/source changed; use a new version')
    for name,digest in lock['outputs'].items():
        if sha(HERE/name)!=digest:
            raise ValueError(f'Frozen fitted output changed: {name}')


def evaluate(motifs,data,cfg,split):
    width=cfg['width']; pattern=cfg['synthetic_pattern']; mutation=data['scenario']['mutation']
    reference=np.full((width,4),mutation/3)
    for i,char in enumerate(pattern): reference[i,ALPHABET.index(char)]=1-mutation
    bg=background(data['discovery_control']); rows=[]; nulls=[]
    matches=[(similarity(pwm,reference,minimum_overlap=width),name) for name,pwm in motifs.items()]
    match,winner=max(matches,key=lambda x:(x[0][0],x[1])) if matches else ((0.,0,'+'),None)
    for name,pwm in motifs.items():
        calibration,_,_=scan(data['calibration_control'],pwm,bg)
        threshold=float(np.quantile(calibration,cfg['evaluation']['calibration_quantile'],method='higher'))
        scores,starts,forward=scan(data[split],pwm,bg); blank,_,_=scan(data[split+'_control'],pwm,bg)
        present=scores>threshold; false=blank>threshold
        nulls.append({'split':split,'candidate':name,'raw_blank_fpr':float(false.mean()),
                      'raw_blank_hits':int(false.sum()),'blank_n':len(blank),'threshold':threshold})
        if name!=winner: continue
        eligible=match[0]>=cfg['evaluation']['reference_cosine_gate']
        actual=np.array([r['start0'] for r in data[split+'_truth']])
        actual_forward=np.array([r['strand']=='+' for r in data[split+'_truth']])
        canonical_forward=forward if match[2]=='+' else ~forward
        for tolerance in cfg['evaluation']['position_tolerances']:
            tp=int((present & (canonical_forward==actual_forward) & (np.abs(starts-actual)<=tolerance)).sum()) if eligible else 0
            predictions=int(present.sum()+false.sum()) if eligible else 0
            precision=tp/predictions if predictions else 0.; recall=tp/len(actual)
            rows.append({'split':split,'tolerance_bp':tolerance,'candidate':winner,'candidate_count':len(motifs),
                'pwm_cosine':match[0],'reference_eligible':eligible,'precision':precision,'recall':recall,
                'f1':2*precision*recall/(precision+recall) if precision+recall else 0.,'tp':tp,
                'prediction_count':predictions,'truth_n':len(actual),'raw_blank_fpr':float(false.mean())})
    if winner is None:
        rows=[{'split':split,'tolerance_bp':t,'candidate':'','candidate_count':0,'pwm_cosine':0.,
               'reference_eligible':False,'precision':0.,'recall':0.,'f1':0.,'tp':0,'prediction_count':0,
               'truth_n':len(data[split]),'raw_blank_fpr':0.} for t in cfg['evaluation']['position_tolerances']]
    return rows,nulls


def score():
    assert_locked(); cfg=read(PROTOCOL); metrics=[]; nulls=[]; instance_diagnostics=[]
    manifest=read(SUMMARY/'manifest.json'); match_rows=[]
    for item in manifest['runs']:
        folder=ART/item['condition']; out=folder/f'{item["method"]}__a{item["algorithm_seed"]}'
        data=read(folder/'data.json'); motifs=read_meme(out/'motifs.meme')
        for split in ['development','test']:
            rows,blank=evaluate(motifs,data,cfg,split)
            metrics.extend({**item,**r} for r in rows); nulls.extend({**item,**r} for r in blank)
        instances=pd.read_csv(out/'instances.tsv',sep='\t')
        truth={r['sequence_index']:r for r in data['discovery_truth']}
        for name,group in instances.groupby('motif'):
            exact=sum(row.start0==truth[row.sequence_index]['start0'] for row in group.itertuples())
            near=sum(abs(row.start0-truth[row.sequence_index]['start0'])<=2 for row in group.itertuples())
            instance_diagnostics.append({**item,'candidate':name,'instances':len(group),
                'exact_core_fraction':exact/len(group),'within2_core_fraction':near/len(group),
                'note':'post-fit diagnostic only; mixed-strand instances not aligned'})
    frame=pd.DataFrame(metrics); frame.to_csv(SUMMARY/'metrics.tsv',sep='\t',index=False)
    pd.DataFrame(nulls).to_csv(SUMMARY/'candidate_null_fpr.tsv',sep='\t',index=False)
    pd.DataFrame(instance_diagnostics).to_csv(SUMMARY/'instance_diagnostics.tsv',sep='\t',index=False)
    primary=frame[(frame.split=='test') & (frame.tolerance_bp==2)]
    per_seed=primary.groupby(['tier','scenario','data_seed','method']).agg(
        f1=('f1','mean'),f1_algorithm_min=('f1','min'),f1_algorithm_max=('f1','max'),
        success_fraction=('reference_eligible','mean'),pwm_cosine=('pwm_cosine','mean')).reset_index()
    per_seed.to_csv(SUMMARY/'per_data_seed.tsv',sep='\t',index=False)
    per_seed.groupby(['tier','scenario','method']).agg(f1_mean=('f1','mean'),
        f1_data_min=('f1','min'),f1_data_max=('f1','max'),data_seeds=('data_seed','nunique'),
        success_fraction=('success_fraction','mean')).reset_index().to_csv(SUMMARY/'test_summary.tsv',sep='\t',index=False)
    pairs=per_seed.pivot(index=['tier','scenario','data_seed'],columns='method',values='f1')
    pairs.subtract(pairs['OH_P24'],axis=0).reset_index().to_csv(SUMMARY/'paired_delta_vs_OH.tsv',sep='\t',index=False)
    pairs.subtract(pairs['LOCAL_P24'],axis=0).reset_index().to_csv(SUMMARY/'paired_delta_vs_LOCAL_P24.tsv',sep='\t',index=False)
    # Algorithm stability: same-data PWM family matching, not independent sample replication.
    for scenario,seed,tier in all_conditions(cfg):
        condition=f'{scenario["id"]}__s{seed}'; folder=ART/condition
        methods=[m for m in cfg['methods'] if tier=='main' or m['id'] in cfg['stress_methods']]
        for method in methods:
            for a,b in itertools.combinations(cfg['algorithm_seeds'],2):
                names_a=read_meme(folder/f'{method["id"]}__a{a}'/'motifs.meme')
                names_b=read_meme(folder/f'{method["id"]}__a{b}'/'motifs.meme')
                if not names_a or not names_b: continue
                an=list(names_a); bn=list(names_b)
                sims=np.array([[similarity(names_a[x],names_b[y],minimum_overlap=cfg['width'])[0] for y in bn] for x in an])
                ia,ib=linear_sum_assignment(-sims)
                af=pd.read_csv(folder/f'{method["id"]}__a{a}'/'instances.tsv',sep='\t')
                bf=pd.read_csv(folder/f'{method["id"]}__a{b}'/'instances.tsv',sep='\t')
                for i,j in zip(ia,ib):
                    aset={(r.sequence_index,r.start0) for r in af[af.motif==an[i]].itertuples()}
                    bset={(r.sequence_index,r.start0) for r in bf[bf.motif==bn[j]].itertuples()}
                    match_rows.append({'condition':condition,'tier':tier,'method':method['id'],'algorithm_a':a,
                        'algorithm_b':b,'candidate_a':an[i],'candidate_b':bn[j],'cosine':sims[i,j],
                        'instance_jaccard':len(aset & bset)/len(aset | bset) if aset|bset else 0.,
                        'position_rule':'exact start, physical coordinates; no claim of direction alignment'})
    pd.DataFrame(match_rows).to_csv(SUMMARY/'algorithm_stability.tsv',sep='\t',index=False)
    manifest.update(status='SCORED',score_elapsed_seconds=time.time()-manifest['started_unix'])
    save(SUMMARY/'manifest.json',manifest); print('Locked unsupervised Test evaluation completed',flush=True)


def probe_examples(features,truth,cfg):
    rows=[]; labels=[]; rng=np.random.default_rng(cfg['pca_seed'])
    for item in truth:
        row=item['sequence_index']; start=item['start0']; width=cfg['width']
        candidates=[i for i in range(82-width) if i+width<=start or i>=start+width]
        negative=rng.choice(candidates,cfg['probe']['negatives_per_sequence'],replace=False)
        rows.append(features[row,start]); labels.append(1)
        rows.extend(features[row,i] for i in negative); labels.extend([0]*len(negative))
    return np.stack(rows),np.array(labels)


def probes():
    assert_locked(); cfg=read(PROTOCOL); model,tokenizer=load_encoder(cfg); records=[]
    for scenario in cfg['scenarios']:
        if scenario['id'] not in cfg['probe']['scenarios']: continue
        for seed in cfg['data_seeds']:
            folder=ART/f'{scenario["id"]}__s{seed}'; data=read(folder/'data.json')
            train=np.load(folder/'local_discovery.npy',mmap_mode='r')
            target=folder/'local_probe_test.npy'
            if target.exists(): raise FileExistsError(target)
            test=local_order(data['test'],model,tokenizer,cfg,target)
            pca=joblib.load(folder/'local_pca.joblib')
            with threadpool_limits(limits=cfg['threads']):
                train_z=pca.transform(train.reshape(-1,train.shape[-1])).reshape(len(train),70,128)
                test_z=pca.transform(test.reshape(-1,test.shape[-1])).reshape(len(test),70,128)
            oh_train=onehot(data['discovery'],cfg['width']); oh_test=onehot(data['test'],cfg['width'])
            for representation in cfg['probe']['representations']:
                if representation=='OH_RAW': x,y=oh_train,oh_test
                elif representation=='LOCAL_RAW': x,y=train,test
                else:
                    dimensions=int(representation.removeprefix('LOCAL_P')); x,y=train_z[:,:,:dimensions],test_z[:,:,:dimensions]
                xx,yy=probe_examples(x,data['discovery_truth'],cfg)
                classifier=make_pipeline(StandardScaler(),LogisticRegression(C=1,class_weight='balanced',
                    max_iter=1000,random_state=cfg['pca_seed']))
                with threadpool_limits(limits=cfg['threads']):
                    classifier.fit(xx,yy)
                    scores=classifier.decision_function(y.reshape(-1,y.shape[-1])).reshape(len(y),70)
                predicted=scores.argmax(axis=1); actual=np.array([r['start0'] for r in data['test_truth']])
                record={'scenario':scenario['id'],'data_seed':seed,'representation':representation,
                    'top1_exact':float((predicted==actual).mean()),'top1_within2':float((np.abs(predicted-actual)<=2).mean()),
                    'test_sequences':len(y),'training_positive_windows':int(yy.sum()),'training_negative_windows':int((yy==0).sum()),
                    'iterations':int(classifier[-1].n_iter_[0]),'note':'supervised ranking diagnostic; no blank/detection FPR'}
                records.append(record); joblib.dump(classifier,folder/f'probe_{representation}.joblib')
                print(f'probe {seed} {representation}: top1±2={record["top1_within2"]:.3f}',flush=True)
            del train,test,train_z,test_z,oh_train,oh_test,x,y; gc.collect()
    pd.DataFrame(records).to_csv(SUMMARY/'supervised_probes.tsv',sep='\t',index=False)
    manifest=read(SUMMARY/'manifest.json'); manifest.update(status='COMPLETED',supervised_probes=len(records),
        elapsed_seconds=time.time()-manifest['started_unix']); save(SUMMARY/'manifest.json',manifest)


def verify():
    assert_locked(); cfg=read(PROTOCOL); manifest=read(SUMMARY/'manifest.json')
    assert len(manifest['runs'])==cfg['expected_discovery_outputs']
    assert manifest['cluster_fits']==cfg['expected_cluster_fits']
    assert manifest['status']=='COMPLETED' and manifest['supervised_probes']==cfg['expected_supervised_probes']
    protected=read(SUMMARY/'protected_before.json')
    assert all(sha(ROOT/name)==digest for name,digest in protected.items())
    instances=0;truths=0;cores=0
    for folder in ART.iterdir():
        data=read(folder/'data.json')
        for split in ['discovery','development','test']:
            for item in data[split+'_truth']:
                assert data[split][item['sequence_index']][item['start0']:item['end0']]==item['fragment']; truths+=1
        for sequence in data['discovery']:
            for start in range(70):
                left,offset=crop_spec(start,12,24)
                assert sequence[left:left+24][offset:offset+12]==sequence[start:start+12]; cores+=1
        for path in folder.glob('*/instances.tsv'):
            for item in pd.read_csv(path,sep='\t').itertuples():
                assert data['discovery'][item.sequence_index][item.start0:item.end0]==item.fragment; instances+=1
    inventory={p.relative_to(HERE).as_posix():sha(p) for p in ART.rglob('*') if p.is_file()}
    inventory.update({p.relative_to(HERE).as_posix():sha(p) for p in SUMMARY.iterdir()
                      if p.is_file() and p.name not in ['artifact_hashes.json','verification.json']})
    save(SUMMARY/'artifact_hashes.json',inventory)
    result={'status':'PASSED','discovery_outputs':len(manifest['runs']),'cluster_fits':manifest['cluster_fits'],
            'supervised_probes':manifest['supervised_probes'],'hashes_checked':len(inventory),
            'protected_files_unchanged':len(protected),'instances_checked':instances,'truth_sites_checked':truths,
            'crop_core_checks':cores,'real_holdout_used':False,
            'limits':['hash consistency, not independent historical reexecution','no exhaustive near-homology audit',
                      'main independent units are five data seeds; stress has only two','STREME not available in current runtime']}
    save(SUMMARY/'verification.json',result); print(result,flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action',choices=['fit','score','probes','verify'])
    action=parser.parse_args().action
    {'fit':fit_all,'score':score,'probes':probes,'verify':verify}[action]()
