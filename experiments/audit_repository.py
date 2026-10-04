"""Audit local source orientation, historical splits and public artifact hashes."""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import re
import subprocess

import pandas as pd
from experiments.preflight import collect_environment
from preprocessing.source_adapters import parse_gff3_attributes
from experiments.validate_regulondb_sigma_holdout import GROUPS, _split_paired_records, _read_fasta, SEED
from experiments.validate_regulondb_cross_dataset import _candidate_index

ROOT = Path(__file__).resolve().parents[1]


def reverse_complement(s):
    return s.translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def dump(p, data):
    p.write_text(json.dumps(data, ensure_ascii=False, indent=2) + '\n', encoding='utf-8')


def table(p, rows, fields):
    with p.open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def run(output):
    output = output.resolve()
    output.relative_to(ROOT)
    if output.exists() and any(output.iterdir()):
        raise FileExistsError(output)
    output.mkdir(parents=True, exist_ok=True)
    tracked = subprocess.check_output(['git', 'ls-files', '-z'], cwd=ROOT).decode('utf-8').split('\0')
    inventory = [{'path': s, 'bytes': (ROOT/s).stat().st_size, 'sha256': sha(ROOT/s)}
                 for s in tracked if s and (ROOT/s).is_file()]
    table(output/'tracked_inventory.tsv', inventory, ['path', 'bytes', 'sha256'])
    checks = []
    for p in (ROOT/'results').rglob('*manifest.json'):
        d = json.loads(p.read_text(encoding='utf-8'))
        for field in ['output_sha256', 'curated_output_sha256', 'summary_sha256', 'matrix_sha256', 'training_motif_matrix_sha256']:
            mapping = d.get(field)
            if not isinstance(mapping, dict):
                continue
            for name, value in mapping.items():
                expected = value.get('sha256') if isinstance(value, dict) else value
                if not isinstance(expected, str) or not re.fullmatch('[0-9a-fA-F]{64}', expected):
                    continue
                candidates = [p.parent/name, ROOT/name, p.parent/Path(name).name]
                # Some manifests refer to matrices produced by a preceding run,
                # rather than to outputs in their own directory.
                if field == 'matrix_sha256' and p.parent.name == 'regulondb_sigma_fimo_20261004':
                    candidates.append(ROOT/'results/motif/m3_20260930'/name)
                if field == 'matrix_sha256' and p.parent.name == 'cross_species':
                    candidates.append(p.parent.parent/'summary/matrices'/(name+'.meme'))
                existing = next((q for q in candidates if q.is_file()), None)
                actual = sha(existing) if existing else ''
                checks.append({'manifest': p.relative_to(ROOT).as_posix(), 'field': field,
                               'artifact': name, 'expected': expected.lower(), 'actual': actual,
                               'status': ('match' if actual == expected.lower() else 'mismatch') if existing else 'not_available_locally'})
    table(output/'artifact_hash_checks.tsv', checks, ['manifest', 'field', 'artifact', 'expected', 'actual', 'status'])
    genome = ''.join(x.strip() for x in (ROOT/'data/raw/regulondb_e_coli_k12_20260916.fna').read_text(encoding='utf-8').splitlines() if not x.startswith('>')).upper()
    orientation = Counter()
    anomalies = []
    for line in (ROOT/'data/raw/regulondb_promoters_20260916.gff3').read_text(encoding='utf-8').splitlines():
        if not line or line.startswith('#'):
            continue
        f = line.split('\t')
        a = parse_gff3_attributes(f[8])
        s = a.get('Sequence', '')
        idx = [i for i, x in enumerate(s) if x in 'ACGT']
        if not s or len(idx) != 1:
            continue
        t, i = int(f[3])-1, idx[0]
        forward = genome[t-i:t-i+len(s)]
        reverse = reverse_complement(genome[t-(len(s)-1-i):t+i+1])
        expected = forward if f[6] == '+' else reverse
        status = 'matches_transcription_direction' if expected == s.upper() else 'source_genome_disagreement'
        orientation[f[6]+'|'+status] += 1
        if status != 'matches_transcription_direction':
            anomalies.append({'source_record_id': a.get('name'), 'strand': f[6], 'genomic_tss': int(f[3]), 'local_tss': i+1})
    profiles = []
    raw = {}
    for name in ['regulondb', 'dbtbs']:
        d = pd.read_csv(ROOT/f'data/interim/{name}_metadata_20260916.tsv', sep='\t')
        labelled = d.dropna(subset=['sequence', 'sigma_factor_type'])
        multigroup = labelled.groupby('sequence').sigma_factor_type.nunique()
        profiles.append({'source': name, 'rows': len(d), 'nonempty_sequences': int(d.sequence.notna().sum()),
                         'missing_tss': int(d.tss_position.isna().sum()), 'missing_strand': int(d.strand.isna().sum()),
                         'missing_evidence': int(d.evidence_level.isna().sum()), 'missing_sigma': int(d.sigma_factor_type.isna().sum()),
                         'multi_sigma_exact_sequence_groups': int((multigroup > 1).sum()),
                         'multi_sigma_records': int(labelled.sequence.isin(multigroup[multigroup > 1].index).sum()),
                         'evidence_counts': d.evidence_level.value_counts().to_dict(),
                         'sigma_counts': d.sigma_factor_type.fillna('unknown').value_counts().to_dict(),
                         'length_counts': {str(k): int(v) for k,v in d.sequence.str.len().value_counts().items()}})
        if name == 'regulondb':
            raw = {r.sequence_id: r for r in d.itertuples()}
    processed = pd.read_csv(ROOT/'data/processed/m1_regulondb_tss_20260922/promoters_clean.tsv', sep='\t')
    states = Counter()
    for r in processed.itertuples():
        original = raw.get(r.sequence_id)
        if original is None:
            states['missing_source_id'] += 1
            continue
        s = str(original.sequence)
        state = 'unchanged' if r.sequence == s else 'reverse_complemented' if r.sequence == reverse_complement(s) else 'other'
        states[str(r.strand)+'|'+state] += 1
    inp = ROOT/'results/motif_inputs_regulondb_tss_n1fix_20260928'
    im = json.loads((inp/'motif_input_manifest.json').read_text(encoding='utf-8'))
    hm = json.loads((ROOT/'results/motif/regulondb_sigma_holdout_20261004/sigma_fimo_manifest.json').read_text(encoding='utf-8'))
    train, test, reconstruct = [], [], []
    for n, g in enumerate(GROUPS):
        desc = im['outputs'][g]
        part = _split_paired_records(_read_fasta(inp/desc['primary_fasta']), _read_fasta(inp/desc['n1_fasta']), seed=SEED+n)
        hashes = {}
        for side in ['train', 'holdout']:
            payload = ''.join('>'+i+'\n'+s+'\n' for i,s in part[side]['primary']).encode('ascii')
            hashes[side] = hashlib.sha256(payload).hexdigest()
        old = hm['holdout_design']['groups'][g]
        match = hashes['train'] == old['train_primary_sha256'] and hashes['holdout'] == old['holdout_primary_sha256']
        reconstruct.append({'group': g, 'train_count': len(part['train']['primary']),
                            'holdout_count': len(part['holdout']['primary']), 'historical_split_hashes_match': match})
        if match:
            train.extend((g,i,s) for i,s in part['train']['primary'])
            test.extend((g,i,s) for i,s in part['holdout']['primary'])
    leakage = []
    seqs = [s for _,_,s in train]
    index = _candidate_index(seqs)
    for tg,ti,s in test:
        candidates = set()
        for query in [s, reverse_complement(s)]:
            for start in range(0,81,9):
                candidates.update(index.get(query[start:start+9], ()))
        for j in sorted(candidates):
            direct = sum(a == b for a,b in zip(s,seqs[j]))
            rev = sum(a == b for a,b in zip(reverse_complement(s),seqs[j]))
            matches = max(direct,rev)
            if matches >= 73:
                leakage.append({'train_group': train[j][0], 'train_id': train[j][1], 'holdout_group': tg,
                                'holdout_id': ti, 'matching_bases': matches, 'length': 81,
                                'orientation': 'direct' if direct >= rev else 'reverse_complement'})
    table(output/'legacy_holdout_near_duplicates.tsv', leakage,
          ['train_group','train_id','holdout_group','holdout_id','matching_bases','length','orientation'])
    dump(output/'data_profiles.json', profiles)
    summary = {'audit_type': 'source/evidence audit; no new motif experiment',
               'git_head': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
               'tracked_files': len(inventory), 'hash_check_counts': dict(Counter(r['status'] for r in checks)),
               'source_orientation': dict(orientation), 'source_genome_anomalies': anomalies,
               'historical_processed_orientation': dict(states), 'legacy_split_reconstruction': reconstruct,
               'legacy_holdout_full_length_90pct_pairs': len(leakage),
               'affected_holdout_records': len(set(r['holdout_id'] for r in leakage)),
               'environment': collect_environment(),
               'scope_limits': ['Local historical genome only; flagged disagreement needs provenance reconciliation.',
                                'Split check covers full-length ungapped 81bp, <=8 substitutions in either orientation; not shifted/indel alignments.',
                                'Missing native XML/site TSV/logs cannot be certified from public summary hashes.']}
    dump(output/'audit_manifest.json', summary)
    print(json.dumps({k:v for k,v in summary.items() if k != 'environment'},ensure_ascii=False,indent=2))
    return summary


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output-dir', type=Path, required=True)
    run(p.parse_args().output_dir)
