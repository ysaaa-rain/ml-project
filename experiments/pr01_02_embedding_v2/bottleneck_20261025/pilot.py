"""Bounded, versioned pilot. Fit uses sequences only; references are evaluation-only."""
import argparse
import gc
import hashlib
import importlib.metadata
import importlib.util
import json
import platform
import shutil
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import requests
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits

from ..discovery import background, make_pwm, read_meme, similarity, write_meme
from ..encoder import token_to_base, windows

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
PROTOCOL = HERE / 'protocol.json'
ARTIFACTS = HERE / 'artifacts'
SUMMARY = HERE / 'summary'
ALPHABET = 'ACGT'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def save(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def read(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def download(url, path, expected=None):
    """Bounded downloads, verified ranges and an optional publisher LFS hash."""
    path = Path(path)
    if path.exists():
        if expected and sha(path) != expected:
            raise ValueError(f'Existing cache hash mismatch: {path.name}')
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    part = path.with_suffix(path.suffix + '.part')
    first = requests.get(url, headers={'Range': 'bytes=0-1048575'}, timeout=(20, 90))
    first.raise_for_status()
    if first.status_code == 206:
        size = int(first.headers['Content-Range'].split('/')[-1])
        end = min(1048575, size - 1)
        if first.headers['Content-Range'] != f'bytes 0-{end}/{size}' or len(first.content) != end + 1:
            raise ValueError('Invalid initial HTTP range')

        def chunk(start):
            stop = min(start + 1048575, size - 1)
            for attempt in range(3):
                try:
                    response = requests.get(url + f'?range_start={start}',
                        headers={'Range': f'bytes={start}-{stop}'}, timeout=(20, 90))
                    response.raise_for_status()
                    if (response.status_code != 206 or
                        response.headers.get('Content-Range') != f'bytes {start}-{stop}/{size}' or
                        len(response.content) != stop - start + 1):
                        raise ValueError('Invalid HTTP range')
                    return response.content
                except (requests.RequestException, ValueError):
                    if attempt == 2:
                        raise
        with part.open('wb') as stream, ThreadPoolExecutor(max_workers=4) as executor:
            stream.write(first.content)
            for content in executor.map(chunk, range(1048576, size, 1048576)):
                stream.write(content)
    else:
        part.write_bytes(first.content)
    if expected and sha(part) != expected:
        raise ValueError(f'Download hash mismatch: {path.name}')
    part.replace(path)


def prepare():
    cfg = read(PROTOCOL)
    info = {'protocol_sha256': sha(PROTOCOL), 'models': {}}
    for key in ['mini', 'mini_c']:
        spec = cfg[key]
        folder = HERE / 'model_cache' / key
        response = requests.get(f'https://huggingface.co/api/models/{spec["repository"]}/revision/{spec["revision"]}',
                                params={'blobs': 'true'}, timeout=30)
        response.raise_for_status()
        model_api = response.json()
        if model_api['sha'] != spec['revision']:
            raise ValueError('Model revision changed')
        weights = next(x for x in model_api['siblings'] if x['rfilename'] == 'model.safetensors')
        expected = weights['lfs']['sha256']
        if key == 'mini' and expected != spec['weights_sha256']:
            raise ValueError('Publisher weight SHA differs from historical frozen SHA')
        for name in ['config.json', 'vocab.txt', 'model.safetensors']:
            cached = HERE.parent / 'model_cache/pinned' / name
            target = folder / name
            if key == 'mini' and cached.exists() and not target.exists():
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(cached, target)
            download(f'https://huggingface.co/{spec["repository"]}/resolve/{spec["revision"]}/{name}',
                     target, expected if name == 'model.safetensors' else None)
        for name in ['models.py', 'tokenizer.py']:
            cached = HERE.parent / 'model_cache/pinned' / name
            target = folder / name
            url = f'https://huggingface.co/neuralbioinfo/nbrg-transformers/resolve/{spec["code_revision"]}/{name}'
            # Verify source against the fixed remote revision even if a historical copy exists.
            response = requests.get(url, timeout=30)
            response.raise_for_status()
            expected_code = hashlib.sha256(response.content).hexdigest()
            if not target.exists():
                target.write_bytes(response.content)
            if sha(target) != expected_code:
                raise ValueError('Pinned implementation mismatch')
        info['models'][key] = {**spec, 'weights_sha256': expected,
            'files': {p.name: sha(p) for p in sorted(folder.iterdir()) if p.is_file()},
            'config': read(folder / 'config.json')}
        print(f'prepared {key}: weights verified', flush=True)
    save(HERE / 'model_info.json', info)


def load_model(key, cfg):
    import torch
    torch.set_num_threads(cfg['threads'])
    folder = HERE / 'model_cache' / key
    info = read(HERE / 'model_info.json')['models'][key]
    for name, digest in info['files'].items():
        if sha(folder / name) != digest:
            raise ValueError('Model cache changed')
    modules = {}
    for name in ['models', 'tokenizer']:
        module_spec = importlib.util.spec_from_file_location(f'emb_pilot_{key}_{name}', folder / f'{name}.py')
        module = importlib.util.module_from_spec(module_spec)
        sys.modules[module_spec.name] = module
        module_spec.loader.exec_module(module)
        modules[name] = module
    tokenizer = modules['tokenizer'].LCATokenizer(vocab_file=str(folder / 'vocab.txt'),
                                                kmer=cfg[key]['kmer'], shift=1)
    config = modules['models'].ProkBertConfig.from_pretrained(str(folder), local_files_only=True)
    model, loading = modules['models'].ProkBertModel.from_pretrained(
        str(folder), config=config, local_files_only=True, output_loading_info=True)
    if loading.get('missing_keys') or loading.get('mismatched_keys') or loading.get('error_msgs'):
        raise ValueError(f'Incomplete pretrained encoder loading: {loading}')
    # Only the pretrained MLM head may be unused when loading the base encoder.
    if any(not key.startswith('cls.') for key in loading.get('unexpected_keys', [])):
        raise ValueError(f'Unexpected non-MLM weights: {loading}')
    model.eval().requires_grad_(False)
    save(SUMMARY / f'{key}_loading.json', {'loading': loading, 'parameters': sum(p.numel() for p in model.parameters()),
        'hidden_size': config.hidden_size, 'layers': config.num_hidden_layers,
        'kmer': cfg[key]['kmer'], 'device': 'cpu', 'implementation_revision': cfg[key]['code_revision']})
    return model, tokenizer


def encode(sequences, model, tokenizer, kmer, batch_size, layer):
    import torch
    blocks = []
    with torch.inference_mode():
        for offset in range(0, len(sequences), batch_size):
            batch = sequences[offset:offset + batch_size]
            inputs = tokenizer(batch, return_tensors='pt', padding=True)
            for ids, sequence in zip(inputs['input_ids'].tolist(), batch):
                expected = [tokenizer.convert_tokens_to_ids(sequence[i:i + kmer])
                            for i in range(len(sequence) - kmer + 1)]
                if ids[1:-1] != expected or tokenizer.unk_token_id in expected:
                    raise ValueError('Token/base mapping mismatch')
            hidden = model(**inputs, output_hidden_states=True).hidden_states[layer][:, 1:-1]
            blocks.append(hidden.numpy().copy())
    return np.concatenate(blocks)


def core_tokens(tokens, width, kmer):
    count = width - kmer + 1
    return np.stack([tokens[:, start:start + count]
                     for start in range(82 - width)], axis=1)


def crop_spec(start, width, context, length=81):
    left = max(0, min(start - (context - width) // 2, length - context))
    return left, start - left


def local_features(sequences, model, tokenizer, cfg):
    width, context = cfg['width'], cfg['local_context']
    nw, nt = 82 - width, width - 5
    ordered = np.empty((len(sequences), nw, nt * model.config.hidden_size), dtype=np.float32)
    mean = np.empty((len(sequences), nw, model.config.hidden_size), dtype=np.float32)
    crops, targets, offsets = [], [], []

    def flush():
        if not crops:
            return
        hidden = encode(crops, model, tokenizer, 6, cfg['batch_size'], cfg['layer'])
        for row, ((sequence_index, start), core_start) in enumerate(zip(targets, offsets)):
            core = hidden[row, core_start:core_start + nt]
            ordered[sequence_index, start] = core.reshape(-1)
            mean[sequence_index, start] = core.mean(axis=0)
        crops.clear(); targets.clear(); offsets.clear()

    began = time.time()
    for sequence_index, sequence in enumerate(sequences):
        for start in range(nw):
            left, offset = crop_spec(start, width, context)
            crops.append(sequence[left:left + context]); targets.append((sequence_index, start)); offsets.append(offset)
            if len(crops) == cfg['batch_size']:
                flush()
        if (sequence_index + 1) % 32 == 0:
            print(f'local crops {sequence_index + 1}/{len(sequences)} ({time.time() - began:.1f}s)', flush=True)
    flush()
    return mean, ordered


def generate(n, scenario, seed, planted, pattern, length=81):
    rng = np.random.default_rng(seed)
    gc_fraction = scenario['gc']
    probability = [(1 - gc_fraction) / 2, gc_fraction / 2, gc_fraction / 2, (1 - gc_fraction) / 2]
    sequences, truth = [], []
    for index in range(n):
        sequence = rng.choice(list(ALPHABET), length, p=probability)
        if planted:
            start = int(rng.integers(0, length - len(pattern) + 1))
            fragment = list(pattern)
            for position, char in enumerate(fragment):
                if rng.random() < scenario['mutation']:
                    fragment[position] = rng.choice([c for c in ALPHABET if c != char])
            sequence[start:start + len(pattern)] = fragment
            truth.append({'sequence_index': index, 'start0': start, 'end0': start + len(pattern),
                          'strand': '+', 'fragment': ''.join(fragment)})
        sequences.append(''.join(sequence))
    return sequences, truth


def create_data(cfg, scenario, seed):
    data = {'scenario': scenario, 'seed': seed}
    # Split streams do not depend on model output, preceding split size, or test results.
    streams = np.random.SeedSequence(seed).spawn(7)
    for index, split in enumerate(['discovery', 'development', 'test']):
        n = cfg[f'{split}_n']
        data[split], data[split + '_truth'] = generate(n, scenario, streams[2 * index], True, cfg['synthetic_pattern'])
        data[split + '_control'], _ = generate(n, scenario, streams[2 * index + 1], False, cfg['synthetic_pattern'])
    data['calibration_control'], _ = generate(cfg['calibration_n'], scenario, streams[6], False, cfg['synthetic_pattern'])
    seen = set()
    for split in ['discovery', 'discovery_control', 'development', 'development_control', 'test', 'test_control', 'calibration_control']:
        canonical = {min(s, reverse_complement(s)) for s in data[split]}
        if len(canonical) != len(data[split]) or seen & canonical:
            raise ValueError('Exact/RC duplicate within or across synthetic partitions')
        seen.update(canonical)
    return data


def reverse_complement(sequence):
    return sequence.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


def fit(features, sequences, cfg, destination):
    """No truth/reference argument; closest window per cluster per sequence."""
    width = cfg['width']; nw = 82 - width
    flat = features.reshape(-1, features.shape[-1])
    with threadpool_limits(limits=cfg['threads']):
        pca = PCA(n_components=min(cfg['pca_components'], flat.shape[-1]),
                  svd_solver='randomized', random_state=cfg['algorithm_seed'])
        reduced = pca.fit_transform(flat)
        km = MiniBatchKMeans(n_clusters=cfg['clusters'], n_init=cfg['n_init'],
            batch_size=2048, random_state=cfg['algorithm_seed'])
        labels = km.fit_predict(reduced)
        distance = np.linalg.norm(reduced - km.cluster_centers_[labels], axis=1)
    joblib.dump({'pca': pca, 'kmeans': km}, destination / 'fitted.joblib')
    motifs, rows = {}, []
    for cluster in range(cfg['clusters']):
        members = np.flatnonzero(labels == cluster)
        selected = {}
        for index in members[np.argsort(distance[members], kind='stable')]:
            selected.setdefault(int(index // nw), int(index % nw))
        if len(selected) < cfg['minimum_promoters']:
            continue
        name = f'C{cluster:02d}'
        fragments = [sequences[i][start:start + width] for i, start in selected.items()]
        _, motifs[name] = make_pwm(fragments)
        rows.extend({'motif': name, 'sequence_index': i, 'start0': start, 'end0': start + width,
                     'strand': '+', 'fragment': sequences[i][start:start + width]}
                    for i, start in selected.items())
    pd.DataFrame(rows).to_csv(destination / 'instances.tsv', sep='\t', index=False)
    return motifs, float(pca.explained_variance_ratio_.sum())


def scan(sequences, pwm, bg):
    encoded = np.array([[ALPHABET.index(c) for c in s] for s in sequences], dtype=np.int8)
    window = np.lib.stride_tricks.sliding_window_view(encoded, len(pwm), axis=1)
    forward = np.log2(pwm / bg)[np.arange(len(pwm)), window].sum(axis=-1)
    reverse = np.log2(pwm[::-1, ::-1] / bg)[np.arange(len(pwm)), window].sum(axis=-1)
    scores = np.concatenate([forward, reverse], axis=1)
    indices = scores.argmax(axis=1)
    return scores[np.arange(len(sequences)), indices], indices % forward.shape[1], indices < forward.shape[1]


def evaluate(motifs, data, cfg, split):
    """Reference-assisted localization plus unconditional per-candidate null FPR."""
    pattern = cfg['synthetic_pattern']; scenario = data['scenario']; width = cfg['width']
    reference = np.full((len(pattern), 4), scenario['mutation'] / 3)
    for i, char in enumerate(pattern):
        reference[i, ALPHABET.index(char)] = 1 - scenario['mutation']
    bg = background(data['discovery_control'])
    matched = [(similarity(pwm, reference, minimum_overlap=width), name) for name, pwm in motifs.items()]
    match, winner = max(matched, key=lambda item: (item[0][0], item[1])) if matched else ((0., 0, '+'), None)
    rows, candidate_rows = [], []
    for name, pwm in motifs.items():
        calibration, _, _ = scan(data['calibration_control'], pwm, bg)
        threshold = float(np.quantile(calibration, cfg['evaluation']['calibration_quantile'], method='higher'))
        scores, starts, forward = scan(data[split], pwm, bg)
        null_scores, _, _ = scan(data[split + '_control'], pwm, bg)
        raw_positive = scores > threshold; raw_null = null_scores > threshold
        candidate_rows.append({'candidate': name, 'split': split, 'threshold': threshold,
            'raw_positive_hits': int(raw_positive.sum()), 'raw_blank_hits': int(raw_null.sum()),
            'blank_n': len(null_scores), 'raw_blank_fpr': float(raw_null.mean())})
        if name != winner:
            continue
        eligible = match[0] >= cfg['evaluation']['reference_cosine_gate']
        actual = np.array([row['start0'] for row in data[split + '_truth']])
        # Full-width matching requires offset=0; direction tracks PWM/reference alignment.
        direction_correct = forward == (match[2] == '+')
        for tolerance in cfg['evaluation']['position_tolerances']:
            tp = int((raw_positive & direction_correct & (np.abs(starts - actual) <= tolerance)).sum()) if eligible else 0
            n_predictions = int(raw_positive.sum() + raw_null.sum()) if eligible else 0
            precision = tp / n_predictions if n_predictions else 0.
            recall = tp / len(actual)
            rows.append({'split': split, 'tolerance_bp': tolerance, 'candidate': winner,
                'candidate_count': len(motifs), 'pwm_cosine': match[0], 'reference_eligible': eligible,
                'precision': precision, 'recall': recall,
                'f1': 2 * precision * recall / (precision + recall) if precision + recall else 0.,
                'true_positives': tp, 'prediction_count': n_predictions, 'truth_n': len(actual),
                'raw_positive_hits': int(raw_positive.sum()), 'raw_blank_hits': int(raw_null.sum()),
                'raw_blank_fpr': float(raw_null.mean())})
    if winner is None:
        rows = [{'split': split, 'tolerance_bp': t, 'candidate': '', 'candidate_count': 0,
                 'pwm_cosine': 0., 'reference_eligible': False, 'precision': 0., 'recall': 0., 'f1': 0.,
                 'true_positives': 0, 'prediction_count': 0, 'truth_n': len(data[split]),
                 'raw_positive_hits': 0, 'raw_blank_hits': 0, 'raw_blank_fpr': 0.}
                for t in cfg['evaluation']['position_tolerances']]
    return rows, candidate_rows


def source_hashes():
    files = [PROTOCOL, HERE / 'pilot.py', HERE / '__init__.py', HERE.parent / 'encoder.py', HERE.parent / 'discovery.py']
    return {str(p.relative_to(ROOT)).replace('\\', '/'): sha(p) for p in files}


def run():
    cfg = read(PROTOCOL)
    ARTIFACTS.mkdir(exist_ok=False)
    SUMMARY.mkdir(exist_ok=True)
    manifest = {'status': 'RUNNING', 'protocol_sha256': sha(PROTOCOL), 'source_sha256': source_hashes(),
        'git_head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], cwd=ROOT, text=True).strip(),
        'started_unix': time.time(), 'real_holdout_used': False,
        'python': platform.python_version(), 'platform': platform.platform(),
        'dependencies': {k: importlib.metadata.version(k) for k in ['torch', 'transformers', 'numpy', 'pandas', 'scipy', 'scikit-learn']},
        'model_info_sha256': sha(HERE / 'model_info.json'), 'runs': []}
    save(SUMMARY / 'run_manifest.json', manifest)
    pd.DataFrame([{'window_start0': start, 'crop_start0': crop_spec(start, cfg['width'], cfg['local_context'])[0],
        'core_start_in_crop0': crop_spec(start, cfg['width'], cfg['local_context'])[1],
        'core_width': cfg['width'], 'crop_length': cfg['local_context']}
        for start in range(82 - cfg['width'])]).to_csv(SUMMARY / 'crop_map.tsv', sep='\t', index=False)
    frozen = ROOT / 'data/processed/pr01_02_data_v2'
    protected = [frozen / 'manifest.json', HERE.parent / 'REPORT.md', HERE.parent / 'next_round/REPORT.md',
                 HERE.parent / 'next_round/verification.json']
    protected += list((HERE.parent / 'runs').glob('*/manifest.json'))
    protected += [frozen / name for name in read(frozen / 'manifest.json')['output_sha256']]
    save(SUMMARY / 'protected_before.json', {str(p.relative_to(ROOT)).replace('\\', '/'): sha(p) for p in protected})
    models = {}
    try:
        for key in ['mini', 'mini_c']:
            models[key] = load_model(key, cfg)
        for scenario in cfg['scenarios']:
            for seed in cfg['data_seeds']:
                tag = f'{scenario["id"]}__s{seed}'
                condition = ARTIFACTS / tag; condition.mkdir()
                data = create_data(cfg, scenario, seed); save(condition / 'data.json', data)
                sequences = data['discovery']; width = cfg['width']; nw = 82 - width
                tokens = encode(sequences, *models['mini'], 6, cfg['batch_size'], cfg['layer'])
                full_core = core_tokens(tokens, width, 6)
                full_mean = full_core.mean(axis=2)
                full_order = full_core.reshape(len(sequences), nw, -1)
                a0 = windows(token_to_base(tokens), width)
                print(f'{tag}: full encoding complete', flush=True)
                local_mean, local_order = local_features(sequences, *models['mini'], cfg)
                minic = encode(sequences, *models['mini_c'], 1, cfg['batch_size'], cfg['layer'])
                minic_mean = core_tokens(minic, width, 1).mean(axis=2)
                onehot = np.array([np.eye(4, dtype=np.float32)[[ALPHABET.index(c) for c in s[i:i + width]]].ravel()
                    for s in sequences for i in range(nw)]).reshape(len(sequences), nw, -1)
                features = dict(zip(cfg['methods'], [onehot, a0, full_mean, full_order, local_mean, local_order, minic_mean]))
                bg = background(data['discovery_control'])
                for method in cfg['methods']:
                    began = time.time(); destination = condition / method; destination.mkdir()
                    motifs, variance = fit(features[method], sequences, cfg, destination)
                    write_meme(destination / 'motifs.meme', motifs, bg)
                    item = {'condition': tag, 'scenario': scenario['id'], 'data_seed': seed, 'method': method,
                        'algorithm_seed': cfg['algorithm_seed'], 'candidate_count': len(motifs),
                        'pca_explained_variance': variance, 'input_dimensions': features[method].shape[-1],
                        'fit_seconds': time.time() - began, 'status': 'COMPLETED'}
                    save(destination / 'fit.json', item); manifest['runs'].append(item)
                    save(SUMMARY / 'run_manifest.json', manifest)
                    print(f'{tag} {method}: {len(motifs)} PWMs ({item["fit_seconds"]:.1f}s)', flush=True)
                del features, onehot, tokens, full_core, full_mean, full_order, a0, local_mean, local_order, minic, minic_mean
                gc.collect()
        save(SUMMARY / 'fit_lock.json', {'protocol_sha256': sha(PROTOCOL), 'source_sha256': source_hashes(),
            'frozen_before_test_evaluation': True,
            'pwm_sha256': {str(p.relative_to(HERE)).replace('\\', '/'): sha(p) for p in ARTIFACTS.glob('*/*/motifs.meme')}})
        manifest.update(status='FITS_COMPLETED', elapsed_seconds=time.time() - manifest['started_unix'])
        save(SUMMARY / 'run_manifest.json', manifest)
        score_all()
    except Exception as error:
        manifest.update(status='FAILED', error=f'{type(error).__name__}: {error}', elapsed_seconds=time.time() - manifest['started_unix'])
        save(SUMMARY / 'run_manifest.json', manifest)
        raise


def score_all():
    cfg = read(PROTOCOL); lock = read(SUMMARY / 'fit_lock.json')
    if lock['protocol_sha256'] != sha(PROTOCOL) or lock['source_sha256'] != source_hashes():
        raise ValueError('Protocol/code changed after fitting; create a new version')
    for name, digest in lock['pwm_sha256'].items():
        if sha(HERE / name) != digest:
            raise ValueError('PWM changed after fit lock')
    metrics, nulls = [], []
    for condition in sorted(ARTIFACTS.iterdir()):
        data = read(condition / 'data.json')
        for method in cfg['methods']:
            destination = condition / method
            motifs = read_meme(destination / 'motifs.meme')
            for split in ['development', 'test']:
                rows, candidate_rows = evaluate(motifs, data, cfg, split)
                common = {'condition': condition.name, 'scenario': data['scenario']['id'], 'data_seed': data['seed'], 'method': method}
                metrics.extend({**common, **row} for row in rows)
                nulls.extend({**common, **row} for row in candidate_rows)
    frame = pd.DataFrame(metrics); frame.to_csv(SUMMARY / 'metrics.tsv', sep='\t', index=False)
    pd.DataFrame(nulls).to_csv(SUMMARY / 'candidate_null_fpr.tsv', sep='\t', index=False)
    primary = frame[(frame.split == 'test') & (frame.tolerance_bp == cfg['evaluation']['primary_tolerance'])]
    primary.groupby(['scenario', 'method']).agg(f1_mean=('f1', 'mean'), f1_min=('f1', 'min'),
        f1_max=('f1', 'max'), reference_passes=('reference_eligible', 'sum'), data_seeds=('data_seed', 'nunique'),
        blank_fpr_mean=('raw_blank_fpr', 'mean')).reset_index().to_csv(SUMMARY / 'test_summary.tsv', sep='\t', index=False)
    paired = primary.pivot(index=['scenario', 'data_seed'], columns='method', values='f1')
    for reference in ['OH', 'MINI_FULL_A0']:
        delta = paired.subtract(paired[reference], axis=0)
        delta.reset_index().to_csv(SUMMARY / f'paired_delta_vs_{reference}.tsv', sep='\t', index=False)
    manifest = read(SUMMARY / 'run_manifest.json')
    manifest.update(status='COMPLETED', evaluation='post-fit, protocol-locked, reference-assisted synthetic localization',
                    real_holdout_used=False, elapsed_seconds=time.time() - manifest['started_unix'])
    save(SUMMARY / 'run_manifest.json', manifest)
    inventory = {str(p.relative_to(HERE)).replace('\\', '/'): sha(p)
                 for p in ARTIFACTS.rglob('*') if p.is_file()}
    inventory.update({str(p.relative_to(HERE)).replace('\\', '/'): sha(p)
                      for p in SUMMARY.iterdir() if p.is_file() and p.name not in ['artifact_hashes.json', 'verification.json']})
    save(SUMMARY / 'artifact_hashes.json', inventory)
    print('All fits and locked evaluation completed', flush=True)


def verify():
    cfg = read(PROTOCOL); manifest = read(SUMMARY / 'run_manifest.json')
    checks = read(SUMMARY / 'artifact_hashes.json')
    mismatches = [name for name, digest in checks.items() if not (HERE / name).is_file() or sha(HERE / name) != digest]
    instance_count = 0; truth_count = 0; crop_count = 0
    for condition in sorted(ARTIFACTS.iterdir()):
        data = read(condition / 'data.json')
        for split in ['discovery', 'development', 'test']:
            for row in data[split + '_truth']:
                assert data[split][row['sequence_index']][row['start0']:row['end0']] == row['fragment']
                truth_count += 1
        for start in range(82 - cfg['width']):
            left, offset = crop_spec(start, cfg['width'], cfg['local_context'])
            for sequence in data['discovery']:
                assert sequence[left:left + cfg['local_context']][offset:offset + cfg['width']] == sequence[start:start + cfg['width']]
                crop_count += 1
        for method in cfg['methods']:
            rows = pd.read_csv(condition / method / 'instances.tsv', sep='\t')
            for row in rows.itertuples():
                assert data['discovery'][row.sequence_index][row.start0:row.end0] == row.fragment
                instance_count += 1
    expected = len(cfg['scenarios']) * len(cfg['data_seeds']) * len(cfg['methods'])
    assert manifest['status'] == 'COMPLETED' and len(manifest['runs']) == expected
    assert all(row['status'] == 'COMPLETED' for row in manifest['runs'])
    if mismatches:
        raise ValueError(f'Hash mismatches: {mismatches}')
    protected = read(SUMMARY / 'protected_before.json')
    assert all(sha(ROOT / name) == digest for name, digest in protected.items())
    result = {'status': 'PASSED', 'hashes_checked': len(checks), 'mismatches': mismatches,
        'fits': expected, 'instances_checked': instance_count, 'truth_sites_checked': truth_count,
        'crop_core_checks': crop_count, 'real_holdout_used': False, 'protected_files_unchanged': len(protected),
        'limits': ['artifact consistency is not an independent rerun', 'two data seeds; one algorithm seed',
                   'no exhaustive synthetic homology audit; exact/RC duplicates checked during generation']}
    save(SUMMARY / 'verification.json', result); print(json.dumps(result), flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['prepare', 'run', 'score', 'verify'])
    action = parser.parse_args().action
    {'prepare': prepare, 'run': run, 'score': score_all, 'verify': verify}[action]()
