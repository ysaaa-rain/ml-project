"""Read-only audit; never rebuild or split the frozen input."""
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
from collections import Counter

import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
BRANCH = Path(__file__).resolve().parent
DATA = ROOT / 'data/processed/pr01_02_data_v2'


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding='utf-8')


def fasta(path):
    records = []
    for line in Path(path).read_text().splitlines():
        if line.startswith('>'):
            records.append([line[1:].split()[0], ''])
        elif line.strip():
            records[-1][1] += line.strip()
    return records


def load(source, split, label=1):
    if split not in ('discovery', 'development'):
        raise ValueError('Holdout is disabled until a separately locked final protocol.')
    frame = pd.read_csv(DATA / source / 'main/sequence_master.tsv', sep='\t')
    return frame[(frame.split == split) & (frame.label == label)].sort_values('sequence_id').reset_index(drop=True)


def run():
    out = BRANCH / 'audit'
    out.mkdir(exist_ok=True)
    manifest = json.loads((DATA / 'manifest.json').read_text(encoding='utf-8'))
    checks = []
    for section, base in [('output_sha256', DATA), ('raw_sha256', ROOT)]:
        for name, expected in manifest[section].items():
            path = base / name
            actual = sha(path) if path.exists() else None
            checks.append({'section': section, 'path': name, 'expected': expected,
                           'actual': actual, 'match': actual == expected})
    pd.DataFrame(checks).to_csv(out / 'hash_checks.tsv', sep='\t', index=False)
    frames, counts, bindings = [], [], []
    for path in sorted(DATA.glob('*/**/sequence_master.tsv')):
        d = pd.read_csv(path, sep='\t')
        frames.append(d)
        counts.extend({'source_tier': str(path.parent.relative_to(DATA)), 'split': s,
                       'label': int(l), 'n': len(g), 'gc_mean': g.sequence.str.count('[GC]').div(81).mean()}
                      for (s, l), g in d.groupby(['split', 'label']))
        assert d.sequence.str.fullmatch('[ACGT]{81}').all()
        assert (d.local_tss_index == 61).all()
        for split in ['discovery', 'development', 'holdout']:
            pos = dict(fasta(path.parent / f'{split}.positive.fasta'))
            expected = dict(zip(d[(d.split == split) & (d.label == 1)].sequence_id,
                                d[(d.split == split) & (d.label == 1)].sequence))
            assert pos == expected
            null = dict(fasta(path.parent / f'{split}.dinucleotide_null.fasta'))
            # Binding is accepted only by exact frozen identifier relation.
            for key, seq in null.items():
                parent = key.removesuffix('__dinucleotide_null')
                if parent not in pos:
                    parent = key.split('|')[0].removesuffix('_dinuc').removesuffix('_shuffle')
                bindings.append({'source_tier': str(path.parent.relative_to(DATA)),
                                 'split': split, 'null_id': key, 'parent_id': parent,
                                 'bound': parent in pos,
                                 'dinucleotide_preserved': parent in pos and
                                 Counter(a+b for a,b in zip(seq,seq[1:])) ==
                                 Counter(a+b for a,b in zip(pos[parent],pos[parent][1:]))})
            natural_path = path.parent/f'{split}.natural_control.fasta'
            if natural_path.exists():
                expected_negative = dict(zip(d[(d.split == split) & (d.label == 0)].sequence_id,
                                             d[(d.split == split) & (d.label == 0)].sequence))
                assert dict(fasta(natural_path)) == expected_negative
    all_data = pd.concat(frames)
    cross = int((all_data.groupby('leakage_group').split.nunique() > 1).sum())
    tju = all_data[all_data.source.str.startswith('tjupan_')]
    assert int((tju.label == 1).sum()) == 9708
    assert int(((tju.label == 1) & (tju.split == 'discovery')).sum()) == 6745
    assert int((tju.label == 0).sum()) == 8641
    aux = {}
    for source,tier in [('regulondb_ecoli','core'),('regulondb_ecoli','extended'),('dbtbs_bsub','core')]:
        folder = DATA/source/tier
        d = pd.read_csv(folder/'sequence_master.tsv',sep='\t')
        sig = pd.read_csv(folder/'sigma_associations.tsv',sep='\t')
        external_ids = [i for i,s in fasta(folder/'external_validation.positive.fasta')]
        ext = d[d.sequence_id.isin(external_ids)]
        assert len(ext) == len(external_ids)
        assert not set(ext.leakage_group)&set(tju.leakage_group)
        assert (ext.split == 'holdout').all()
        aux[f'{source}/{tier}'] = {'positive_n':int((d.label==1).sum()),'external_candidates':len(ext),
              'external_group_overlap_all_tju_splits':0,
              'multi_sigma_sequences':int((sig.groupby('sequence_id').sigma_id.nunique()>1).sum()),
              'sigma_evidence_fields':sig.columns.tolist(),
              'independence_basis':'frozen exhaustive <=8-edit both-strand leakage graph; not newly recomputed alignments'}
    traditional = json.loads((ROOT/'results/motif/tju_streme_v2_20261008/public_evidence/run_manifest.json').read_text(encoding='utf-8'))
    assert traditional['input_manifest_sha256'] == sha(DATA/'manifest.json')
    save(out/'auxiliary_and_traditional.json', {'auxiliary':aux,'streme_status':traditional['status'],
         'streme_main_candidates':traditional['main_candidates'], 'streme_runs':len(traditional['runs']),
         'expanded_manifest_present':(ROOT/'tmp/m3_expanded_20261009/manifest.json').exists(),
         'm1m2_missing_stats':'preserved original pending status'})
    pd.DataFrame(counts).to_csv(out / 'counts_gc.tsv', sep='\t', index=False)
    pd.DataFrame(bindings).to_csv(out / 'shuffle_bindings.tsv', sep='\t', index=False)
    import torch
    env = {'python': platform.python_version(), 'platform': platform.platform(),
           'cuda_available': torch.cuda.is_available(),
           'packages': {p: importlib.metadata.version(p) for p in
                        ['torch', 'transformers', 'numpy', 'pandas', 'scipy', 'scikit-learn', 'matplotlib', 'huggingface-hub']}}
    save(out / 'environment.json', env)
    references = {}
    from pypdf import PdfReader
    for p in ROOT.glob('docs/**/motif-embedding*.pdf'):
        text = '\n'.join(page.extract_text() or '' for page in PdfReader(p).pages)
        (out / (p.stem + '.txt')).write_text(text, encoding='utf-8')
        references[str(p.relative_to(ROOT))] = {'sha256': sha(p), 'pages': len(PdfReader(p).pages),
                                              'unresolved_citation_markers': text.count('TRAE_REF')}
    result = {'status': 'COMPLETED' if all(c['match'] for c in checks) and cross == 0 else 'FAILED',
              'data_manifest_sha256': sha(DATA / 'manifest.json'), 'hash_files': len(checks),
              'hash_failures': [c for c in checks if not c['match']], 'cross_split_groups': cross,
              'shuffle_unbound': sum(not b['bound'] for b in bindings),
              'shuffle_dinucleotide_failures':sum(not b['dinucleotide_preserved'] for b in bindings),
              'references': references, 'git_head': subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
              'limitations': ['Upstream TJU TSS alignment accepted, not independently reconstructed.',
                             'Reference PDFs contain unresolved citations; hypotheses only.',
                             'DBTBS external holdout n=1 is not statistical validation.',
                             'Existing traditional holdout has been viewed; future EMB evaluation is exploratory.']}
    if result['shuffle_unbound'] or result['shuffle_dinucleotide_failures']:
        result['status'] = 'FAILED'
    save(out / 'audit.json', result)
    print(json.dumps(result, ensure_ascii=True))
    return result


if __name__ == '__main__':
    run()
