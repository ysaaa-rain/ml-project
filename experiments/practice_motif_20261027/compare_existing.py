"""Summarize existing partial-evaluation results only within matched batches."""
import hashlib
import json
from pathlib import Path

import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SOURCE = ROOT / 'reports/practice_motif_20261027/synthetic_partial_per_data_seed.tsv'
OUT = ROOT / 'reports/practice_comparison_20261027'
CONFIG = {
    'budget_validation_20261026': [('BG_ONLY', 'EMB128'), ('BG_ONLY', 'OH24'), ('BG_ONLY', 'OH48')],
    'cnn_seed_validation_20261027': [('standard', 'OH48_BG'), ('standard', 'TIED_RAW'), ('standard', 'TIED_SHIFT')],
}

def main():
    data = pd.read_csv(SOURCE, sep='\t')
    print('columns:', data.columns.tolist())
    groups = []
    for version, conditions in CONFIG.items():
        batch = data[(data.version == version) & (data.view == 'CAP5') & (data.split == 'test')]
        for variant, representation in conditions:
            subset = batch[(batch.variant == variant) & (batch.representation == representation)].copy()
            if subset.empty:
                raise ValueError(f'missing comparison: {version} {variant} {representation}')
            groups.append(subset)
    selected = pd.concat(groups, ignore_index=True)
    for (version, scenario), group in selected.groupby(['version', 'scenario']):
        seeds = [set(g.data_seed) for _, g in group.groupby(['variant', 'representation'])]
        if len(seeds) != len(CONFIG[version]) or any(s != seeds[0] for s in seeds):
            raise ValueError(f'unmatched data seeds: {version}, {scenario}')
    OUT.mkdir(parents=True, exist_ok=False)
    selected.to_csv(OUT / 'per_data_seed.tsv', sep='\t', index=False, lineterminator='\n')
    summary = selected.groupby(['version','scenario','variant','representation'], as_index=False).agg(
        f1_mean=('f1','mean'), f1_min=('f1','min'), f1_max=('f1','max'),
        blank_fpr=('family_blank_fpr','mean'), n_data=('data_seed','nunique'))
    summary.to_csv(OUT / 'summary.tsv', sep='\t', index=False, lineterminator='\n')
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.5), sharey=True)
    for ax, (version, conditions) in zip(axes, CONFIG.items()):
        part = summary[summary.version == version]
        scenarios = sorted(part.scenario.unique())
        for i, (variant, representation) in enumerate(conditions):
            x = [j + (i-1)*0.23 for j in range(len(scenarios))]
            vals = [part[(part.scenario==s)&(part.variant==variant)&(part.representation==representation)].f1_mean.iloc[0] for s in scenarios]
            ax.bar(x, vals, width=.22, label=f'{representation} / {variant}')
        ax.set_xticks(range(len(scenarios)), [s.replace('_','\n') for s in scenarios], fontsize=8)
        ax.set_title('BERT batch' if 'budget' in version else 'CNN batch')
        ax.set_ylim(0,1); ax.grid(axis='y',alpha=.2); ax.legend(fontsize=8)
    axes[0].set_ylabel('Reference-assisted partial localization F1')
    fig.suptitle('Same evaluation rule; comparisons within each batch only')
    fig.tight_layout(); fig.savefig(OUT / 'matched_batches.png', dpi=160); plt.close(fig)
    manifest = {'source': str(SOURCE.relative_to(ROOT)), 'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
                'source_protocol': 'reports/practice_motif_20261027/manifest.json', 'split': 'test', 'view': 'CAP5',
                'selection': CONFIG, 'rows': len(selected), 'summary_rows':len(summary),
                'limitations': 'post-hoc reference-assisted partial matching; batches use different data; CNN has supervised training'}
    (OUT / 'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n', encoding='utf8')
    print(summary.to_string(index=False))

if __name__ == '__main__':
    main()
