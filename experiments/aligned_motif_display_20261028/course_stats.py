"""Small read-only M1/M2 supplement; no motif model discovery or split changes."""
from collections import Counter
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from .pipeline import ROOT, BASE, HERE, read, save, sha, fasta, unchanged

OUT = ROOT / 'results/eda/m1m2_course_20261028'
PATTERNS = {'template_10': 'TATAAT', 'template_35': 'TTGACA'}


def hit_positions(sequence, pattern, mismatch):
    return [i for i in range(len(sequence) - len(pattern) + 1)
            if sum(a != b for a, b in zip(sequence[i:i + len(pattern)], pattern)) <= mismatch]


def main():
    unchanged()
    OUT.mkdir(exist_ok=False)
    inventory, qc, sigma_rows, hit_rows, position_rows, pair_rows = [], [], [], [], [], []
    inputs = {}
    for cohort in sorted(BASE.glob('*/**/sequence_master.tsv')):
        folder = cohort.parent
        # Only canonical cohort tables, not auxiliary per-sigma FASTA directories.
        source = folder.parent.name
        tier = folder.name
        if tier not in ['main', 'core', 'extended']:
            continue
        key = folder.relative_to(BASE).as_posix()
        inputs[cohort.relative_to(ROOT).as_posix()] = sha(cohort)
        data = pd.read_csv(cohort, sep='\t', keep_default_na=False)
        assert data.sequence.str.fullmatch('[ACGT]{81}').all()
        assert data.local_tss_index.eq(61).all()
        assert data.groupby('leakage_group').split.nunique().max() == 1
        sigma_file = folder / 'sigma_associations.tsv'
        mapping = {}
        if sigma_file.exists():
            inputs[sigma_file.relative_to(ROOT).as_posix()] = sha(sigma_file)
            assoc = pd.read_csv(sigma_file, sep='\t', keep_default_na=False)
            positive_ids = set(data.loc[data.label.eq(1), 'sequence_id'])
            assoc = assoc[assoc.sequence_id.isin(positive_ids) &
                          ~assoc.sigma_id.str.lower().str.endswith(':unknown') & assoc.sigma_id.ne('')]
            mapping = assoc.groupby('sequence_id').sigma_id.apply(lambda x: sorted(set(x))).to_dict()
            for sigma in sorted(set(assoc.sigma_id)):
                ids = {k for k, values in mapping.items() if sigma in values}
                single = {k for k in ids if len(mapping[k]) == 1}
                sigma_rows.append({'cohort': key, 'sigma': sigma, 'any_association_sequences': len(ids),
                                   'single_sigma_sequences': len(single), 'multi_sigma_sequences': len(ids - single)})
        counts = Counter(data.label.astype(int))
        inventory.append({'cohort': key, 'sequences': len(data), 'positive': counts[1], 'natural': counts[0],
                          'leakage_groups': data.leakage_group.nunique(),
                          'sigma_labeled_sequences': len(mapping),
                          'multi_sigma_sequences': sum(len(v) > 1 for v in mapping.values()),
                          'missing_sigma_sequences': int(data.label.eq(1).sum()) - len(mapping),
                          'tss_missing': int(data.local_tss_index.ne(61).sum())})
        groups = [('overall', data)]
        for sigma in sorted({v[0] for v in mapping.values() if len(v) == 1}):
            ids = {k for k, values in mapping.items() if values == [sigma]}
            groups.append(('single:' + sigma, data[data.sequence_id.isin(ids)]))
        for group, frame in groups:
            for (split, label), subset in frame.groupby(['split', 'label']):
                gc = subset.sequence.str.count('[GC]').div(81)
                common = {'cohort': key, 'group': group, 'split': split, 'label': int(label)}
                qc.append({**common, 'sequences': len(subset), 'groups': subset.leakage_group.nunique(),
                           'length_min': int(subset.sequence.str.len().min()), 'length_max': int(subset.sequence.str.len().max()),
                           'gc_mean': float(gc.mean()), 'gc_std': float(gc.std(ddof=0))})
            # Template scans use Discovery and Development. Holdout contributes only QC counts.
            for split in ['discovery', 'development']:
                subset = frame[frame.split == split]
                for label, chunk in subset.groupby('label'):
                    views = [('positive' if label == 1 else 'natural', chunk.sequence.tolist())]
                    if label == 1:
                        null_path = folder / f'{split}.dinucleotide_null.fasta'
                        inputs[null_path.relative_to(ROOT).as_posix()] = sha(null_path)
                        paired = {k.removesuffix('__dinucleotide_null'): seq for k, seq in fasta(null_path)}
                        views.append(('paired_shuffle', [paired[k] for k in chunk.sequence_id]))
                    for view, sequences in views:
                        for mismatch in [0, 1]:
                            meta = {'cohort': key, 'group': group, 'split': split, 'view': view, 'mismatches_le': mismatch}
                            by_template = {}
                            for name, pattern in PATTERNS.items():
                                hits = [hit_positions(seq, pattern, mismatch) for seq in sequences]
                                by_template[name] = hits
                                total = Counter(start - 60 for row in hits for start in row)
                                expected = (-15, -6) if name == 'template_10' else (-40, -29)
                                hit_rows.append({**meta, 'template': name, 'sequences': len(sequences),
                                                 'sequences_with_hit': sum(bool(row) for row in hits), 'all_sites': sum(map(len, hits)),
                                                 'expected_start_low': expected[0], 'expected_start_high': expected[1],
                                                 'sequences_with_expected_start': sum(any(expected[0] <= x - 60 <= expected[1] for x in row) for row in hits)})
                                position_rows.extend({**meta, 'template': name, 'start_tss': start, 'sites': count}
                                                     for start, count in sorted(total.items()))
                            gaps, cooccur, reversed_or_overlap, all_pairs = [], 0, 0, 0
                            for left, right in zip(by_template['template_35'], by_template['template_10']):
                                if left and right:
                                    cooccur += 1
                                pairs = [(a, b) for a in left for b in right]
                                all_pairs += len(pairs)
                                eligible = [(a, b) for a, b in pairs if b >= a + 6]
                                reversed_or_overlap += len(pairs) - len(eligible)
                                if eligible:
                                    # Fixed rule: expected starts, closeness, then leftmost; no outcome tuning.
                                    a, b = min(eligible, key=lambda p: (
                                        -(int(-40 <= p[0] - 60 <= -29) + int(-15 <= p[1] - 60 <= -6)),
                                        abs(p[0] - 60 + 35) + abs(p[1] - 60 + 12), p[0], p[1]))
                                    gaps.append(b - a - 6)
                            pair_rows.append({**meta, 'sequences': len(sequences), 'sequences_with_both_templates': cooccur,
                                              'all_unfiltered_pairs': all_pairs, 'overlap_or_wrong_order_pairs': reversed_or_overlap,
                                              'sequences_with_selected_pair': len(gaps),
                                              'gap_median': float(np.median(gaps)) if gaps else None,
                                              'gap_min': min(gaps) if gaps else None, 'gap_max': max(gaps) if gaps else None})
    outputs = [('m1_source_inventory.tsv', inventory), ('m1_sigma_coverage.tsv', sigma_rows), ('m2_gc_length.tsv', qc),
               ('m2_known_element_hits.tsv', hit_rows), ('m2_tss_position.tsv', position_rows), ('m2_pair_spacing.tsv', pair_rows)]
    for name, rows in outputs:
        pd.DataFrame(rows).to_csv(OUT / name, sep='\t', index=False, lineterminator='\n')
    # Independently recompute exact occurrences with str.startswith for every aggregate input sequence.
    exact_checked = 0
    for path in sorted(BASE.glob('*/**/sequence_master.tsv')):
        if path.parent.name not in ['main', 'core', 'extended']:
            continue
        frame = pd.read_csv(path, sep='\t')
        for sequence in frame[frame.split.isin(['discovery', 'development'])].sequence:
            for pattern in PATTERNS.values():
                assert hit_positions(sequence, pattern, 0) == [i for i in range(76) if sequence.startswith(pattern, i)]
                exact_checked += 1
    table = pd.DataFrame(inventory).set_index('cohort')
    assert table.missing_sigma_sequences.ge(0).all()
    assert sum(table.loc[k, 'positive'] for k in table.index if k.startswith('tjupan_')) == 9708
    assert sum(table.loc[k, 'natural'] for k in table.index if k.startswith('tjupan_')) == 8641
    assert table.loc['regulondb_ecoli/core', 'positive'] == 1793
    assert table.loc['dbtbs_bsub/core', 'positive'] == 432
    overview = pd.DataFrame(qc)
    chart = overview[(overview.group == 'overall') & (overview.split == 'discovery')]
    fig, ax = plt.subplots(figsize=(11, 4.5))
    for label, color, offset in [(1, '#467397', -.18), (0, '#c59643', .18)]:
        rows = chart[chart.label == label].set_index('cohort').reindex(table.index)
        ax.bar(np.arange(len(table)) + offset, rows.gc_mean, width=.35, label='positive' if label else 'natural', color=color)
    ax.set_xticks(range(len(table)), [x.replace('tjupan_', '').replace('/main', '') for x in table.index], rotation=30, ha='right', fontsize=8)
    ax.set_ylabel('Mean GC fraction (Discovery)')
    ax.set_ylim(0, 1)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(OUT / 'data_overview.png', dpi=170)
    plt.close(fig)
    save(OUT / 'input_sha256.json', inputs)
    save(OUT / 'audit.json', {'status': 'PASSED', 'source_code_sha256': sha(Path(__file__)),
                             'frozen_outputs_checked': 210, 'exact_scan_independent_checks': exact_checked,
                             'new_model_training': False, 'template_splits': ['discovery', 'development'],
                             'holdout_use': 'counts/length/GC QC only; no template scan, candidate generation or model evaluation'})
    text = ['# M1/M2 课程数据统计补充', '',
            '统一读取冻结v2，补充三来源、六物种及来源标签分组的数量、GC/长度、两个固定字符串的出现位置与间距。仅为课程数据分析，无模型训练或新候选发现。', '',
            '| 来源/层级 | 正类 | 天然负类 | 同源组 |', '| --- | ---: | ---: | ---: |']
    for row in inventory:
        text.append(f"| {row['cohort']} | {row['positive']} | {row['natural']} | {row['leakage_groups']} |")
    text += ['', '![Discovery数据组成](data_overview.png)', '',
             '所有纳入序列均为81 bp，TSS局部第61位。来源与层级不合并计独立样本；RegulonDB extended单列，DBTBS与TJU的跨来源重叠保留。天然正负类组成差异按实际值展示。', '',
             '模板仅在转录正向输入直接匹配TATAAT与TTGACA，分别报告exact与≤1错配；位点起点相对TSS为start0−60。预先固定的宽描述区间为TATAAT起点[−15,−6]、TTGACA起点[−40,−29]，仅用作统一描述字段；TJU缺逐条sigma标签，其他sigma组不据此解释为适用结构。', '',
             '间距表同时保留未过滤共现、全部组合数及重叠/反序组合数；每序列只选一个不重叠、顺序正确的组合。优先两个描述窗口内的起点，其次距离−35/−12起点的绝对差总和，再按左端坐标。gap=start10−start35−6，不临时按gap大小删值。', '',
             '全部split只做数量、长度与GC质检；模板位置/间距仅用Discovery和Development，没有新增Holdout模型评价。每来源标签多对多关系保留，single分组只纳入恰有一个不同sigma标签的序列。旧组员表的过滤定义不完全相同，本补充单独保存且不覆盖旧表，最终展示采用这里统一定义的表格。', '',
             '完整表：[数量](m1_source_inventory.tsv)、[标签覆盖](m1_sigma_coverage.tsv)、[GC/长度](m2_gc_length.tsv)、[固定模板命中](m2_known_element_hits.tsv)、[位置](m2_tss_position.tsv)、[间距](m2_pair_spacing.tsv)、[核验](audit.json)。']
    (OUT / 'findings.md').write_text('\n'.join(text) + '\n', encoding='utf-8')
    save(OUT / 'manifest.json', {'source_sha256': inputs, 'output_sha256': {p.name: sha(p) for p in OUT.iterdir() if p.is_file()}})
    unchanged()
    print('COURSE_STATS', len(inventory), 'cohorts;', exact_checked, 'independent exact scans', flush=True)


if __name__ == '__main__':
    main()
