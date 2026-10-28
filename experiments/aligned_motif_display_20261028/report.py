"""Build report-facing tables, genuine PWM logos and same-scanner position plots."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D
import numpy as np
import pandas as pd

from .pipeline import HERE, ROOT, OUT, DATA, pools, read, save, sha, relative, unchanged

NAMES = ['MEME', 'STREME', 'ONEHOT', 'PROKBERT', 'CNN_SHIFT']
COLORS = ['#467397', '#69a6ad', '#c59643', '#8d79a4', '#58966d']


def logo(ax, pwm):
    colors = dict(A='#40916c', C='#4d80b0', G='#d29b38', T='#c96460')
    font = FontProperties(family='DejaVu Sans', weight='bold')
    with np.errstate(divide='ignore', invalid='ignore'):
        bits = 2 + np.sum(np.where(pwm > 0, pwm * np.log2(pwm), 0), axis=1)
    for i, (row, information) in enumerate(zip(pwm, bits)):
        bottom = 0.
        for index in np.argsort(row):
            height = float(row[index] * information)
            if height < .005:
                continue
            char = 'ACGT'[index]
            path = TextPath((0, 0), char, size=1, prop=font)
            bounds = path.get_extents()
            transform = (Affine2D().translate(-bounds.x0, -bounds.y0)
                         .scale(.88 / bounds.width, height / bounds.height)
                         .translate(i + .06, bottom))
            ax.add_patch(PathPatch(path, transform=transform + ax.transData, color=colors[char], lw=0))
            bottom += height
    ax.set_xlim(0, len(pwm))
    ax.set_ylim(0, 2)
    ax.set_xticks(np.arange(len(pwm)) + .5, np.arange(1, len(pwm) + 1), fontsize=7)
    ax.set_yticks([0, 1, 2], ['0', '1', '2'], fontsize=7)
    ax.spines[['top', 'right']].set_visible(False)


def main():
    unchanged()
    assert not (OUT / 'REPORT.md').exists(), 'Refuse to overwrite the report'
    family = pd.read_csv(OUT / 'family.tsv', sep='\t', dtype={'seed': str})
    ranking = pd.read_csv(OUT / 'discovery_rankings.tsv', sep='\t', dtype={'seed': str})
    positions = pd.read_csv(OUT / 'positions.tsv', sep='\t', dtype={'seed': str})
    counts = pd.read_csv(OUT / 'data_counts.tsv', sep='\t')
    primary = family[(family.purpose == 'primary') & (family.view == 'TOP5')]
    conditions = {'MEME': 'positive vs natural-control differential discovery',
                  'STREME': 'positive vs natural-control differential discovery',
                  'ONEHOT': 'Discovery-positive unsupervised clustering',
                  'PROKBERT': 'Discovery-positive frozen features and unsupervised clustering',
                  'CNN_SHIFT': 'Discovery-positive/natural sequence-label supervised training'}
    shared = {name + '_sha256': sha(DATA / (name + '.fasta')) for name in
              ['discovery.positive', 'discovery.natural_control', 'development.positive',
               'development.natural_control', 'development.dinucleotide_null']}
    alignment = [{'method': row.method, 'seed': row.seed, 'fit_split': 'discovery',
                  'evaluation_split': 'development', 'learning_condition': conditions[row.method],
                  'scan': 'same dual-strand zero-order PWM scanner', 'threshold': 'Discovery natural q0.9 higher strict >',
                  'candidate_budget': 5, 'positive_groups': row.positive_groups, 'natural_groups': row.natural_groups,
                  **shared} for row in primary.itertuples()]
    pd.DataFrame(alignment).to_csv(OUT / 'alignment.tsv', sep='\t', index=False, lineterminator='\n')
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.1), sharey=True)
    for ax, column, title in zip(axes, ['positive_rate', 'natural_rate', 'shuffle_rate'],
                                 ['Positive groups', 'Natural negative groups', 'Paired shuffled groups']):
        for i, (name, color) in enumerate(zip(NAMES, COLORS)):
            values = primary[primary.method == name][column].to_numpy()
            ax.bar(i, values.mean(), color=color, width=.65)
            ax.scatter(i + np.linspace(-.06, .06, len(values)), values, s=26, color='#27313b', zorder=3)
            ax.text(i, min(1.015, max(values) + .035), f'{values.mean():.2f}', ha='center', fontsize=9)
        ax.set_title(title, fontsize=11)
        ax.set_ylim(0, 1.1)
        ax.set_xticks(range(len(NAMES)), NAMES, rotation=30, ha='right', fontsize=9)
        ax.grid(axis='y', alpha=.18)
        ax.spines[['top', 'right']].set_visible(False)
    axes[0].set_ylabel('Group hit fraction, common TOP5 budget')
    fig.suptitle('Five methods on the same B. subtilis Development data', fontsize=13)
    fig.text(.5, .015, 'Dots: fixed algorithm seeds on one dataset. CNN uses labels. Descriptive course comparison.', ha='center', fontsize=8)
    fig.tight_layout(rect=(0, .04, 1, .95))
    fig.savefig(OUT / 'comparison.png', dpi=180)
    plt.close(fig)

    primary_pools = [x for x in pools() if x[2] == 'primary']
    fig, axes = plt.subplots(len(primary_pools), 3, figsize=(12, len(primary_pools) * 1.8), squeeze=False)
    selected = []
    for row, (method, seed, _, _, motifs) in enumerate(primary_pools):
        top = ranking[(ranking.method == method) & (ranking.seed == seed)].sort_values('rank').head(3)
        for col, item in enumerate(top.itertuples()):
            ax = axes[row, col]
            logo(ax, motifs[item.candidate])
            ax.set_title(f'#{item.rank}  {item.candidate}', fontsize=9, loc='left')
            if col == 0:
                ax.set_ylabel(f'{method}\n{seed}\nBits', fontsize=8)
            selected.append({'method': method, 'seed': seed, 'rank': item.rank, 'candidate': item.candidate,
                             'selection': 'Discovery ranking only'})
    fig.suptitle('Actual PWM logos: first three Discovery-ranked candidates from every run', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, .975), h_pad=1.3)
    fig.savefig(OUT / 'candidate_logos.png', dpi=180)
    plt.close(fig)
    pd.DataFrame(selected).to_csv(OUT / 'figure_candidates.tsv', sep='\t', index=False, lineterminator='\n')

    fig, axes = plt.subplots(len(primary_pools), 3, figsize=(12, len(primary_pools) * 1.75), sharex=True, sharey=True, squeeze=False)
    dev = counts[counts.split == 'development'].iloc[0]
    denominators = {'positive': dev.positive_groups, 'natural': dev.natural_groups, 'shuffle': dev.positive_groups}
    bins = np.arange(-60.5, 22, 2)
    for row, (method, seed, _, _, _) in enumerate(primary_pools):
        top = ranking[(ranking.method == method) & (ranking.seed == seed)].sort_values('rank').head(3)
        for col, item in enumerate(top.itertuples()):
            ax = axes[row, col]
            for label, color in [('positive', '#467397'), ('natural', '#9b9fa6'), ('shuffle', '#c59643')]:
                frame = positions[(positions.method == method) & (positions.seed == seed) &
                                  (positions.candidate == item.candidate) & (positions['class'] == label)]
                values, _ = np.histogram(frame.center_tss, bins=bins, weights=frame.groups / denominators[label])
                ax.stairs(values, bins, label=label, color=color, lw=1.2)
            ax.axvline(0, color='#aaa', lw=.6, ls='--')
            ax.set_xlim(-60, 20)
            ax.set_title(f'#{item.rank}  {item.candidate}', fontsize=9, loc='left')
            ax.grid(axis='y', alpha=.15)
            if col == 0:
                ax.set_ylabel(f'{method}\n{seed}', fontsize=8)
            if row == len(primary_pools) - 1:
                ax.set_xlabel('Best-site center relative to TSS (bp)', fontsize=8)
    fig.legend(*axes[0, 0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(.5, .976), ncol=3, frameon=False)
    fig.suptitle('Common Development scanner: group hit fraction per 2 bp bin', fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, .953), h_pad=1.2)
    fig.savefig(OUT / 'candidate_positions.png', dpi=180)
    plt.close(fig)

    mean = primary.groupby('method')[['positive_rate', 'natural_rate', 'natural_difference', 'shuffle_rate']].mean().reindex(NAMES)
    mean.to_csv(OUT / 'method_means.tsv', sep='\t', lineterminator='\n')
    lines = ['# 机器学习课程实验：同数据五方法候选展示', '',
             '本轮以枯草芽孢杆菌为完整案例，补齐冻结 Discovery 上的 MEME 后，对 MEME、STREME、One-hot、ProkBERT 和 CNN 使用同一套评价程序。已有模型及历史结果保留，没有新增 BERT/CNN/One-hot 训练。', '',
             '## 数据与横向对齐', '',
             '逐方法输入哈希、学习条件、评价分区和候选预算见 [横向对齐验收表](alignment.tsv)。输入/评价口径一致，方法自身的训练方式保留。', '',
             '| 分区 | 正类序列 | 天然负类序列 | 配对Shuffle序列 | 正类有效组 | 天然负类有效组 |',
             '| --- | ---: | ---: | ---: | ---: | ---: |']
    for row in counts.itertuples():
        lines.append(f'| {row.split} | {row.positive_sequences} | {row.natural_sequences} | {row.paired_shuffle_sequences} | {row.positive_groups} | {row.natural_groups} |')
    lines += ['', '输入为既有 v2、81 bp、TSS局部第61位。候选产生和共同阈值设定仅用 Discovery；Development 用于课程探索展示，未新增 Holdout 评价。各候选统一采用双链扫描、同一零阶背景和 Discovery 天然负类最大分数的90%分位阈值（higher，严格大于）。按 Discovery 正类减天然负类的同源组命中率排序，信息量和候选ID作为并列顺序。天然比较排除混标签组，Shuffle始终与正类配对。', '',
              'TOP5用于等候选预算比较；ALL用于展示实际产出规模。单候选10%附近的校准比例不等于候选集合10%误报率，多候选合并可能增加对照命中。共同扫描分数与阈值不称作FIMO p/q，也不更改历史严格统计。', '',
              '## 五方法比较', '', '![同数据TOP5比较](comparison.png)', '',
              '柱为既定运行均值，点为每次运行。MEME/STREME各1次固定运行，One-hot、BERT、CNN各2个既定种子；种子共享数据。CNN使用序列级标签监督，方法生成候选的过程和监督条件不同。', '',
              '| 方法 | 种子 | 正类组命中 | 天然负类组命中 | 差值 | 配对Shuffle组命中 |',
              '| --- | --- | ---: | ---: | ---: | ---: |']
    for row in primary.itertuples():
        lines.append(f'| {row.method} | {row.seed} | {row.positive_rate:.3f} | {row.natural_rate:.3f} | {row.natural_difference:+.3f} | {row.shuffle_rate:.3f} |')
    lines += ['', '表中命中率描述候选对不同数据的覆盖，不能只根据正类覆盖率给方法排优劣；同时观察天然负类和配对打乱对照。两种子差异反映已有算法运行的波动，不作为独立数据重复。', '',
              '## 候选与位置', '', '![真实PWM Logo](candidate_logos.png)', '',
              '每次运行固定展示 Discovery 前3候选，字母高度来自实际PWM与信息量；未按Development表现挑选。候选ID、排序及阈值见 [完整候选表](candidates.tsv) 和 [图中候选来源](figure_candidates.tsv)。', '',
              '![统一位置展示](candidate_positions.png)', '',
              '所有方法都使用同一个 Development 扫描器；每候选每序列取最高分位点，再取每同源组最高分位点。各2 bp区间以该类全部有效组为分母，展示命中的分布。位置图不是已知真值的定位准确率。', '',
              '允许偏移和正反向的局部PWM相似结果见 [similarity.tsv](similarity.tsv)，保留重叠长度及方向；它用于展示不同方法候选是否相似。', '',
              '## 补跑、范围与收束', '',
              '本轮仅新增 MEME 5.5.9 的两套固定参数发现：Discovery正类分别对天然负类与既有dinucleotide-null。天然背景结果进入主比较，打乱背景结果与STREME打乱背景另列稳健性视图；不根据Development表现切换主候选池。CNN RAW与SHIFT另作已有结果诊断，见 [全部候选集合评价](family.tsv)。', '',
              '10月3日旧MEME六物种60个PWM继续保留；该批使用全量数据，包含当前Development和Holdout，因此未混入本轮同划分比较。完整五方法对齐范围为枯草芽孢杆菌；其余五物种已有STREME、One-hot和BERT候选的覆盖情况见 [coverage.tsv](coverage.tsv)，CNN未运行、新版MEME未恢复的格子如实留缺。', '',
              '本轮以课程展示完整性作为停止条件：数据和规则对齐、真实候选可追溯、图表与表格一致即可收束。BERT保持阶段收束；CNN现有批次作为课程方法案例，不继续为了提高分数扩大调参。', '',
              '[协议](../protocol.json)｜[数据计数](data_counts.tsv)｜[输入审计](audit.json)｜[MEME运行记录](meme_runs.json)｜[锁定排序](discovery_rankings.tsv)｜[结果核验](verification.json)。']
    native = read(OUT / 'meme_runs.json')
    recovered = [r for r in native['runs'] if r.get('postprocess_exit_code') == 0]
    if recovered:
        lines += ['', '运行说明：天然背景首次搜索已达到10个候选并生成完整XML，但HTML导出因缺少XML::Parser依赖而返回1。安装依赖后仅重做HTML导出并成功，原退出码和日志保留；没有重跑搜索、选择参数或更换PWM。']
    (OUT / 'REPORT.md').write_text('\n'.join(lines) + '\n', encoding='utf-8')
    output_hashes = {p.relative_to(OUT).as_posix(): sha(p) for p in OUT.rglob('*')
                     if p.is_file() and p.name not in ['manifest.json', 'verification.json']}
    save(OUT / 'manifest.json', {'status': 'COMPLETED', 'protocol_sha256': sha(HERE / 'protocol.json'),
                                 'report_source_sha256': sha(Path(__file__)), 'outputs_sha256': output_hashes,
                                 'primary_methods': NAMES, 'primary_runs': len(primary), 'new_meme_runs': 2,
                                 'new_cnn_bert_onehot_training': 0, 'holdout_scored': False,
                                 'scope': 'B. subtilis complete case; other species inventory only'})
    print('REPORT_READY', OUT / 'REPORT.md', flush=True)


if __name__ == '__main__':
    main()
