"""Audit, complete and display existing motif methods on one frozen course dataset.

This module never fits BERT/CNN/One-hot, changes a split, or scores Holdout.
Native MEME output stays local; only matrices and aggregate display tables publish.
"""
import argparse
from collections import Counter
import hashlib
import json
import math
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

from ..pr01_02_embedding_v2.audit import fasta
from ..pr01_02_embedding_v2.discovery import similarity, write_meme
from ..pr01_02_embedding_v2.cnn_real_pilot_20261027.pipeline import (
    background, groups, load_split, rank, rates, read_meme, scan_windows,
    sequence_windows, validate_lock,
)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
BASE = ROOT / 'data/processed/pr01_02_data_v2'
DATA = BASE / 'tjupan_bacillus_subtilis/main'
BRANCH = ROOT / 'experiments/pr01_02_embedding_v2'
STREME = ROOT / 'results/motif/tju_streme_v2_20261008'
OUT = HERE / 'summary'
ART = HERE / 'artifacts'
PROTOCOL = HERE / 'protocol.json'
DISTRO = 'mlproject-meme-20261028'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def save(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False) + '\n', encoding='utf-8')


def sha(path):
    with path.open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def relative(path):
    return path.relative_to(ROOT).as_posix()


def tsv(name, rows):
    pd.DataFrame(rows).to_csv(OUT / name, sep='\t', index=False, lineterminator='\n')


def checked_pwm(path):
    motifs = read_meme(path)
    assert motifs, f'No PWM: {path}'
    for name, pwm in motifs.items():
        assert pwm.ndim == 2 and pwm.shape[1] == 4, (path, name)
        assert np.isfinite(pwm).all() and (pwm >= 0).all(), (path, name)
        assert np.allclose(pwm.sum(axis=1), 1, atol=2e-5), (path, name)
    return motifs


def audit():
    OUT.mkdir(exist_ok=False)
    frozen = read(BASE / 'manifest.json')
    protected = {relative(BASE / name): digest for name, digest in frozen['output_sha256'].items()}
    assert len(protected) == 210 and frozen['cross_split_leakage_groups'] == 0
    for name, digest in protected.items():
        assert sha(ROOT / name) == digest, name
    old = BRANCH / 'real_bsub_comparison_20261027/summary'
    old_manifest = read(old / 'manifest.json')
    for name, digest in old_manifest['source_sha256'].items():
        assert sha(ROOT / name) == digest, name
        protected[name] = digest
    for name, digest in old_manifest['outputs_sha256'].items():
        p = old / name
        assert sha(p) == digest, name
        protected[relative(p)] = digest
    cnn = validate_lock()
    for name, digest in cnn['source_sha256'].items():
        protected[name] = digest
    for name, digest in cnn['artifact_sha256'].items():
        protected[relative(BRANCH / 'cnn_real_pilot_20261027' / name)] = digest
    streme_manifest_path = STREME / 'public_evidence/run_manifest.json'
    streme_manifest = read(streme_manifest_path)
    protected[relative(streme_manifest_path)] = sha(streme_manifest_path)
    inventory = []
    for cohort in sorted(BASE.glob('tjupan_*/main')):
        species = cohort.parent.name
        master = cohort / 'sequence_master.tsv'
        for control, directory in [('natural_control', 'natural'), ('dinucleotide_null', 'dinucleotide')]:
            path = STREME / species / directory / 'streme.txt'
            native = next(x for x in streme_manifest['runs']
                          if x['species'] == species and x['background'] == directory)
            assert sha(cohort / 'discovery.positive.fasta') == native['positive_sha256']
            assert sha(cohort / f'discovery.{control}.fasta') == native['negative_sha256']
            assert sha(path) == streme_manifest['output_sha256'][f'{species}/{directory}/streme.txt']
            protected[relative(path)] = sha(path)
            inventory.append({'species': species, 'method': 'STREME', 'seed': 'existing',
                              'training_control': control, 'pwm_count': len(checked_pwm(path)),
                              'aligned_to_frozen_data': True, 'status': 'available', 'path': relative(path)})
        for method in ['onehot', 'prokbert']:
            for seed in [20261005, 20261006, 20261007]:
                folder = BRANCH / 'runs' / f'{species}__{method}__L6W10__s{seed}'
                manifest = read(folder / 'manifest.json')
                path = folder / 'motifs.meme'
                assert not manifest['holdout_used'] and manifest['status'] == 'COMPLETED'
                assert sha(path) == manifest['output_sha256']['motifs.meme']
                key = relative(master).replace('/', '\\')
                assert sha(master) == manifest['input_hashes'][key]
                protected[relative(path)] = sha(path)
                protected[relative(folder / 'manifest.json')] = sha(folder / 'manifest.json')
                inventory.append({'species': species, 'method': method.upper(), 'seed': str(seed),
                                  'training_control': 'natural_control_evaluation', 'pwm_count': len(checked_pwm(path)),
                                  'aligned_to_frozen_data': True, 'status': 'available', 'path': relative(path)})
        for control in ['natural_control', 'dinucleotide_null']:
            inventory.append({'species': species, 'method': 'MEME', 'seed': '20261005',
                              'training_control': control, 'pwm_count': 0, 'aligned_to_frozen_data': False,
                              'status': 'pending_rerun' if species == 'tjupan_bacillus_subtilis' else 'not_in_bounded_run',
                              'path': ''})
        for seed in [20261005, 20261006]:
            path = BRANCH / 'cnn_real_pilot_20261027/artifacts' / f'seed_{seed}/SHIFT/motifs.meme'
            available = species == 'tjupan_bacillus_subtilis'
            inventory.append({'species': species, 'method': 'CNN_SHIFT', 'seed': str(seed),
                              'training_control': 'natural_control_supervised',
                              'pwm_count': len(checked_pwm(path)) if available else 0,
                              'aligned_to_frozen_data': available, 'status': 'available' if available else 'not_run',
                              'path': relative(path) if available else ''})
    legacy = ROOT / 'results/motif/tjupan_m3_20261003/summary/matrices'
    for path in legacy.glob('*.meme'):
        protected[relative(path)] = sha(path)
    save(OUT / 'protected_history.json', protected)
    tsv('coverage.tsv', inventory)
    counts = []
    for split in ['discovery', 'development']:
        pos, neg, shuffled, mixed = load_split(split)
        assert not set(pos.leakage_group) & set(load_split('development' if split == 'discovery' else 'discovery')[0].leakage_group)
        counts.append({'split': split, 'positive_sequences': len(pos), 'natural_sequences': len(neg),
                       'paired_shuffle_sequences': len(shuffled), 'mixed_label_groups': len(mixed),
                       'positive_groups': len(set(pos.leakage_group) - mixed),
                       'natural_groups': len(set(neg.leakage_group) - mixed)})
    tsv('data_counts.tsv', counts)
    sources = {relative(PROTOCOL): sha(PROTOCOL), relative(Path(__file__)): sha(Path(__file__)),
               relative(BASE / 'manifest.json'): sha(BASE / 'manifest.json')}
    # Pin the actual shared loader, scorer, ranking and CNN artifacts' provenance.
    for name in ['audit.py', 'discovery.py', 'cnn_real_pilot_20261027/pipeline.py',
                 'refinement_20261026/pipeline.py']:
        path = BRANCH / name
        sources[relative(path)] = sha(path)
    save(OUT / 'audit.json', {'status': 'PASSED', 'frozen_outputs_checked': 210,
                             'protected_files': len(protected), 'source_sha256': sources,
                             'holdout_scored': False, 'new_model_training': False})
    print('AUDITED', len(protected), 'protected files;', len(inventory), 'coverage rows', flush=True)


def unchanged():
    for path, digest in read(OUT / 'protected_history.json').items():
        assert sha(ROOT / path) == digest, f'Historical input/output changed: {path}'
    for path, digest in read(OUT / 'audit.json')['source_sha256'].items():
        assert sha(ROOT / path) == digest, f'Locked implementation changed: {path}'


def linux(path):
    path = path.resolve()
    return f'/mnt/{path.drive[0].lower()}' + path.as_posix()[2:]


def discover():
    unchanged()
    ART.mkdir(exist_ok=False)
    cfg = read(PROTOCOL)['meme']
    binary = '/opt/meme-5.5.9/bin/meme'
    version = subprocess.check_output(['wsl', '-d', DISTRO, '--exec', binary, '-version'], text=True).strip()
    assert version == cfg['version'], version
    native = {'version': version, 'new_run': True, 'historical_output_replacement': False,
              'input_split': 'discovery', 'holdout_scored': False, 'runs': []}
    for control in ['natural_control', 'dinucleotide_null']:
        folder = ART / control
        folder.mkdir()
        positive = DATA / 'discovery.positive.fasta'
        negative = DATA / f'discovery.{control}.fasta'
        # Reproduce the native background writer in run_expanded_m3.py.
        counts = Counter(''.join(seq for _, seq in fasta(negative)))
        size = sum(counts.values())
        bg = np.array([((counts[base] + counts[comp]) / 2 + .25) / (size + 1)
                       for base, comp in zip('ACGT', 'TGCA')])
        bfile = folder / 'discovery_background.txt'
        bfile.write_text('# Discovery control, zero order, strand symmetric\n' +
                         '\n'.join(f'{base} {value:.12f}' for base, value in zip('ACGT', bg)) + '\n', encoding='utf-8')
        command = [binary, linux(positive), '-dna', '-revcomp', '-mod', cfg['model'],
                   '-objfun', cfg['objective'], '-test', cfg['test'], '-neg', linux(negative),
                   '-bfile', linux(bfile), '-markov_order', str(cfg['markov_order']),
                   '-hsfrac', str(cfg['hsfrac']), '-searchsize', str(cfg['searchsize']),
                   '-nmotifs', str(cfg['nmotifs']), '-minw', str(cfg['minw']), '-maxw', str(cfg['maxw']),
                   '-seed', str(cfg['seed']), '-maxsize', str(cfg['maxsize']), '-o', linux(folder / 'native')]
        before = time.monotonic()
        with (folder / 'stdout.log').open('wb') as stdout, (folder / 'stderr.log').open('wb') as stderr:
            result = subprocess.run(['wsl', '-d', DISTRO, '--exec', *command], stdout=stdout, stderr=stderr)
        record = {'control': control, 'exit_code': result.returncode, 'elapsed_seconds': round(time.monotonic() - before, 3),
                  'positive_n': len(fasta(positive)), 'negative_n': len(fasta(negative)),
                  'positive_sha256': sha(positive), 'negative_sha256': sha(negative),
                  'command': [binary, relative(positive), *command[2:]], 'output': relative(folder / 'native')}
        # Keep commands portable and do not publish machine-specific paths.
        record['command'] = [str(x).replace(linux(ROOT) + '/', '') for x in record['command']]
        native['runs'].append(record)
        save(OUT / 'meme_runs.json', native)
        assert result.returncode == 0, f'MEME failed; see {folder}/stderr.log'
        path = folder / 'native/meme.txt'
        motifs = checked_pwm(path)
        record['motifs'] = len(motifs)
        record['native_sha256'] = {p.name: sha(p) for p in (folder / 'native').iterdir() if p.is_file()}
        matrices = OUT / 'matrices'
        matrices.mkdir(exist_ok=True)
        exported = matrices / f'meme_{control}.meme'
        write_meme(exported, motifs, bg)
        assert all(np.allclose(value, checked_pwm(exported)[key], atol=1e-8) for key, value in motifs.items())
        record['matrix_sha256'] = sha(exported)
        print('MEME', control, len(motifs), 'PWM;', record['elapsed_seconds'], 's', flush=True)
        save(OUT / 'meme_runs.json', native)
    native['status'] = 'COMPLETED'
    native['tool_sha256'] = subprocess.check_output(['wsl', '-d', DISTRO, '--exec', 'sha256sum', binary], text=True).split()[0]
    save(OUT / 'meme_runs.json', native)
    unchanged()


def pools():
    items = [('MEME', '20261005', 'primary', OUT / 'matrices/meme_natural_control.meme'),
             ('STREME', 'existing', 'primary', STREME / 'tjupan_bacillus_subtilis/natural/streme.txt')]
    for method in ['onehot', 'prokbert']:
        for seed in read(PROTOCOL)['algorithm_seeds']:
            items.append((method.upper(), str(seed), 'primary', BRANCH / 'runs' /
                          f'tjupan_bacillus_subtilis__{method}__L6W10__s{seed}/motifs.meme'))
    for variant in ['SHIFT', 'RAW']:
        for seed in read(PROTOCOL)['algorithm_seeds']:
            items.append(('CNN_' + variant, str(seed), 'primary' if variant == 'SHIFT' else 'diagnostic',
                          BRANCH / 'cnn_real_pilot_20261027/artifacts' / f'seed_{seed}/{variant}/motifs.meme'))
    items.extend([('MEME_NULL', '20261005', 'robustness', OUT / 'matrices/meme_dinucleotide_null.meme'),
                  ('STREME_NULL', 'existing', 'robustness', STREME / 'tjupan_bacillus_subtilis/dinucleotide/streme.txt')])
    return [(method, seed, purpose, path, checked_pwm(path)) for method, seed, purpose, path in items]


def evaluate(motifs, ranking, pos, neg, shuffled, mixed, bg, positions=False):
    windows, hits, candidate_rows, locations = {}, {}, [], []
    for item in ranking:
        name = item['candidate']
        pwm = motifs[name]
        width = len(pwm)
        if width not in windows:
            windows[width] = [sequence_windows(seq, width) for seq in [pos.sequence.tolist(), neg.sequence.tolist(), shuffled]]
        scanned = [scan_windows(window, pwm, bg) for window in windows[width]]
        hits[name] = [result[0] > item['threshold'] for result in scanned]
        pr, nr, pn, nn = rates(pos, neg, hits[name][0], hits[name][1], mixed)
        candidate_rows.append({**item, 'positive_groups': pn, 'natural_groups': nn,
                               'positive_rate': pr, 'natural_rate': nr,
                               'shuffle_rate': float(groups(pos, hits[name][2]).mean())})
        if positions and item['rank'] <= 3:
            for label, frame, result, flags in zip(['positive', 'natural', 'shuffle'], [pos, neg, pos], scanned, hits[name]):
                table = pd.DataFrame({'group': frame.leakage_group.to_numpy(), 'id': frame.sequence_id.to_numpy(),
                                      'score': result[0], 'start': result[1], 'forward': result[2], 'hit': flags})
                if label != 'shuffle':
                    table = table[~table.group.isin(mixed)]
                table = table[table.hit].sort_values(['score', 'id'], ascending=[False, True]).drop_duplicates('group')
                centers = table.start + (width - 1) / 2 - 60
                for center, count in centers.value_counts().sort_index().items():
                    locations.append({'candidate': name, 'rank': item['rank'], 'class': label,
                                      'center_tss': float(center), 'groups': int(count)})
    families = []
    for view, names in [('TOP5', [x['candidate'] for x in ranking[:5]]), ('ALL', [x['candidate'] for x in ranking])]:
        flags = [np.logical_or.reduce([hits[name][i] for name in names]) for i in range(3)]
        pr, nr, pn, nn = rates(pos, neg, flags[0], flags[1], mixed)
        pg, sg = groups(pos, flags[0]), groups(pos, flags[2])
        families.append({'view': view, 'candidates': len(names), 'positive_groups': pn, 'natural_groups': nn,
                         'positive_rate': pr, 'natural_rate': nr, 'natural_difference': pr - nr,
                         'shuffle_rate': float(sg.mean()), 'paired_positive_rate': float(pg.mean()),
                         'paired_positive_only': int((pg & ~sg).sum()), 'paired_shuffle_only': int((sg & ~pg).sum()),
                         'paired_groups': len(pg)})
    return candidate_rows, families, locations


def score():
    unchanged()
    assert read(OUT / 'meme_runs.json')['status'] == 'COMPLETED'
    assert not (OUT / 'fit_lock.json').exists(), 'Refuse to overwrite a locked comparison'
    pos, neg, _, mixed = load_split('discovery')
    bg = background(neg.sequence.tolist())
    all_pools = pools()
    rankings, sources = [], {}
    for method, seed, purpose, path, motifs in all_pools:
        sources[relative(path)] = sha(path)
        ranked = rank(motifs, pos, neg, bg, .9, mixed)
        assert len(ranked) >= 5
        rankings.extend({'method': method, 'seed': seed, 'purpose': purpose, 'rank': i, **row}
                        for i, row in enumerate(ranked, 1))
    tsv('discovery_rankings.tsv', rankings)
    save(OUT / 'fit_lock.json', {'source_sha256': sources, 'ranking_sha256': sha(OUT / 'discovery_rankings.tsv'),
                                'protocol_sha256': sha(PROTOCOL), 'background': bg.tolist(), 'holdout_scored': False})
    dp, dn, shuffle, dmixed = load_split('development')
    candidates, family, locations, matches = [], [], [], []
    references = [(method, motifs) for method, _, purpose, _, motifs in all_pools
                  if purpose == 'primary' and method in ['MEME', 'STREME']]
    for method, seed, purpose, _, motifs in all_pools:
        ranked = [x for x in rankings if x['method'] == method and x['seed'] == seed]
        cr, fr, loc = evaluate(motifs, ranked, dp, dn, shuffle, dmixed, bg, positions=purpose == 'primary')
        candidates.extend(cr)
        family.extend({'method': method, 'seed': seed, 'purpose': purpose, **x} for x in fr)
        locations.extend({'method': method, 'seed': seed, **x} for x in loc)
        for item in ranked[:5]:
            query = motifs[item['candidate']]
            for refmethod, refs in references:
                if method == refmethod:
                    continue
                choices = []
                for refname, pwm in refs.items():
                    minimum = max(5, math.ceil(.6 * min(len(query), len(pwm))))
                    value, offset, strand = similarity(query, pwm, minimum_overlap=minimum)
                    overlap = min(len(query), offset + len(pwm)) - max(0, offset)
                    choices.append((value, refname, overlap, offset, strand, len(pwm)))
                value, refname, overlap, offset, strand, width = max(choices, key=lambda x: (x[0], x[1]))
                matches.append({'method': method, 'seed': seed, 'candidate': item['candidate'], 'reference_method': refmethod,
                                'reference_candidate': refname, 'cosine': value, 'overlap_bp': overlap,
                                'offset': offset, 'orientation': strand, 'query_width': len(query), 'reference_width': width})
    tsv('candidates.tsv', candidates)
    tsv('family.tsv', family)
    tsv('positions.tsv', locations)
    tsv('similarity.tsv', matches)
    coverage = pd.read_csv(OUT / 'coverage.tsv', sep='\t', keep_default_na=False)
    native = read(OUT / 'meme_runs.json')
    for record in native['runs']:
        mask = (coverage.species == 'tjupan_bacillus_subtilis') & (coverage.method == 'MEME') & (coverage.training_control == record['control'])
        coverage.loc[mask, ['pwm_count', 'aligned_to_frozen_data', 'status', 'path']] = [record['motifs'], True, 'completed_new_run',
                                                                                     relative(OUT / f"matrices/meme_{record['control']}.meme")]
    coverage.to_csv(OUT / 'coverage.tsv', sep='\t', index=False, lineterminator='\n')
    print(pd.DataFrame(family).query("view == 'TOP5' and purpose == 'primary'").to_string(index=False), flush=True)
    unchanged()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['audit', 'discover', 'score'])
    args = parser.parse_args()
    {'audit': audit, 'discover': discover, 'score': score}[args.action]()
