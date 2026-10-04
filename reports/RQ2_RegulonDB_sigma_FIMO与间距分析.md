# RQ2：RegulonDB 六个 σ 组 FIMO 与间距分析

> 2026-10-04审计说明：旧留出检出4对跨训练/测试近重复，且输入负链二次翻转、multi_sigma丢失；不能作新规同源隔离或严格位置/语法证据。 [审计证据](../docs/08_全方位核查审计_20261004.md)。

更新日期：2026-10-04

## 问题与数据

本实验补充回答“不同 σ 因子组的 motif 命中频率和组合间距是否不同”。TJU Pan 六物种没有 σ 标签，因此本实验使用辅助数据 RegulonDB *E. coli* K-12 六组：Sigma70 1,942 条、Sigma24 513 条、Sigma32 300 条、Sigma38 176 条、Sigma28 140 条、Sigma54 96 条。每条正序列配一条保留二核苷酸频数的 N1 对照，共 3,167 对。组间未发现同向或反向互补的完全相同正序列。

## 方法

将六组各自 MEME 得到的 10 个矩阵（共 60 个）分别扫描全部六组正序列和 N1 对照；每个来源组运行一次 FIMO 5.5.9。统一使用合并 N1 对照估计的零阶 A/C/G/T 背景、双链扫描、motif pseudocount 0.1，并输出 raw site p≤0.05 的位点及 FIMO q 值。

- **组内 motif 命中：** 对每个 motif、目标 σ 组和 raw site p≤0.05/0.01/0.001 阈值，比较正序列与其一一对应 N1 对照是否命中，采用精确 McNemar 检验；每档对 360 项“来源 motif×目标组”检验分别做 BH。
- **σ 组特异性：** 对每个 motif 比较目标组正序列与其余五组正序列的命中率，采用双侧 Fisher 检验；每档对 360 项检验分别做 BH。
- **motif 组合与间距：** 对每组自身的 10 个 motif 共 45 对，按同样阈值比较正序列与配对 N1 的共同命中，精确 McNemar 检验；每档对 270 项“组×motif 对”检验做 BH。另对共同命中的序列报告最近位点对的有符号间隔；重叠为负值，间隔是未匹配碱基数。

矩阵在原发现组数据上发现，随后也扫描了该组数据；因此来源组命中和间距属于探索性结果，需留出或独立数据复核。

## 结果

| 位点 raw p 阈值 | 来源组正序列高于配对 N1 的 motif（全部 / E-value≤0.05） | 来源组命中率高于其他五组的候选 motif（E-value≤0.05） | 正序列共现高于 N1 且 pair BH q≤0.05 |
| ---: | ---: | ---: | ---: |
| 0.05 | 0 / 0 | 2 / 14 | 0 / 270 |
| 0.01 | 22 / 13 | 10 / 14 | 129 / 270 |
| 0.001 | 35 / 14 | 11 / 14 | 134 / 270 |

“来源组命中率高于其他五组”是按 motif 来源组与其余组比较，再从 14 个 MEME E-value≤0.05 候选中筛出。按组拆分的候选数量为：

| σ 组 | 来源组高于 N1，p≤0.01 | 来源组高于其他组，p≤0.01 | 来源组高于 N1，p≤0.001 | 来源组高于其他组，p≤0.001 |
| --- | ---: | ---: | ---: | ---: |
| Sigma70 | 5 | 2 | 5 | 2 |
| Sigma24 | 3 | 3 | 4 | 4 |
| Sigma32 | 2 | 2 | 2 | 2 |
| Sigma38 | 1 | 1 | 1 | 1 |
| Sigma28 | 1 | 1 | 1 | 1 |
| Sigma54 | 1 | 1 | 1 | 1 |

在 p≤0.01 时，129 个通过校正且正序列共现较多的 motif 对中，正序列相对 N1 的共现率差中位数为 12.5 个百分点；p≤0.001 时，134 对的差值中位数为 4.51 个百分点。p≤0.05 时共现命中接近饱和，没有 pair 通过校正。

各组 45 个 pair 的“最近间隔中位数”再取组内中位数如下；该统计只描述位点间距，不是组间差异检验：

| σ 组 | p≤0.01 | p≤0.001 |
| --- | ---: | ---: |
| Sigma70 | 0 bp | 8 bp |
| Sigma24 | 1 bp | 11 bp |
| Sigma32 | 1 bp | 10 bp |
| Sigma38 | 0 bp | 6 bp |
| Sigma28 | 0 bp | 13 bp |
| Sigma54 | 1 bp | 14.5 bp |

严格 FIMO site q≤0.05 的命中较少；上述组间与组合结果主要来自 raw site p 阈值分析。14 个候选 motif 中，来源组对 N1 的富集有 13 个在 p≤0.01、14 个在 p≤0.001 通过校正；来源组相对其他 σ 组的命中率差异分别有 10 个和 11 个通过校正。结果显示候选 motif 的命中频率具有 σ 组相关性，但其强度和组合间距随 site 阈值变化。

## 对 RQ2 的阶段回答

现有结果支持“不同 σ 组的候选 motif 命中频率和组合背景存在差异”这一探索性结论。尚不能确认稳定的 σ 特异性调控语法：motif 与扫描序列有重叠；raw p≤0.05/0.01/0.001 得出的显著 pair 数不同；各组间距只作描述，未对组间间距分布做独立推断。下一步用留出序列复核 motif 频率和 spacing/order，并与 TJU Pan 主分析的结论分开报告。

## 同来源留出复核

为避免在同一批序列上发现并检验，六组分别按配对启动子/N1 对照划分训练集和留出集：训练集约 80%，留出集为 `ceil(0.20×N)`。训练集 MEME 每组输出 10 个 motif；FIMO 只扫描未参与发现的留出序列，共 3,167 对中的 636 对。FIMO 汇总表含 360 行 motif 覆盖、1,080 行组特异性及 810 行 motif-pair 共现/间距统计。拆分种子、逐组样本数、软件参数和文件哈希见留出复核 manifest。

15 个训练集 motif 的 MEME E-value≤0.05。对 15 个 motif 的来源组留出序列，raw site p≤0.01 时正类命中率合计为 79.7%，配对 N1 为 74.6%，配对检验经每档 360 项“来源 motif×目标 σ 组”BH 校正后没有一项显著；p≤0.001 时两者分别为 23.2% 和 15.6%，3/15 个来源组比较通过校正。p≤0.05 时两类命中率均约 99%，接近饱和，不能区分正类与背景。严格 FIMO site q≤0.05 下，仅 Sigma70 MEME-10（正类 4 条、N1 2 条）和 Sigma54 MEME-1（正类 1 条、N1 0 条）在正类留出序列出现命中。

来源组与其余五组比较时（每档对 360 项检验做 BH），Sigma70 的 MEME-10 和 Sigma24 的 MEME-3 在 p≤0.01、0.001 两档均通过；p≤0.01 命中率分别为 65.6% 对 51.4%、94.2% 对 76.2%。其余候选未通过。270 个组内 motif pair 中（每档对 270 项检验做 BH），p≤0.01 有 13 个共现率高于 N1 且通过校正，均来自 Sigma70（9 个）和 Sigma24（4 个）；正类与 N1 共现率差中位数为 13.4 个百分点，最近有符号间隔中位数为 0 bp。p≤0.05 和 p≤0.001 均无 pair 通过校正。

留出复核表明，少数 motif 的命中频率和组合在同一数据来源的新序列中可重现；多数候选未通过校正，pair 结果只在中间阈值成立。RQ2 目前只能回答为“存在部分 σ 组相关候选”，不足以确认普遍、稳定的 σ 特异 motif 或调控语法。该实验只验证 RegulonDB 内部留出，不等同于独立数据集或物种验证。

结果文件位于 `results/motif/regulondb_sigma_holdout_20261004/`：覆盖、组特异性和 pair 间距表，六组训练 PWM 矩阵及 `sigma_fimo_manifest.json`。完整输入拆分和逐位点 FIMO 记录保存在项目内忽略目录 `tmp/regulondb_sigma_holdout_20261004/`。

## 结果文件与复现

- `results/motif/regulondb_sigma_fimo_20261004/sigma_motif_coverage.tsv`：60 个来源 motif 在六组正序列/N1 的命中率与配对检验。
- `results/motif/regulondb_sigma_fimo_20261004/sigma_motif_group_specificity.tsv`：目标 σ 组与其他五组的命中率比较。
- `results/motif/regulondb_sigma_fimo_20261004/sigma_motif_pair_spacing.tsv`：270 组内 motif 对的共现、配对背景检验和最近间距。
- `results/motif/regulondb_sigma_fimo_20261004/sigma_fimo_manifest.json`：软件、输入与矩阵哈希、参数、运行退出码、行数及输出哈希。
- 脚本：`experiments/analyze_regulondb_sigma_fimo.py`。

复现命令：

```bash
.venv/bin/python -m experiments.analyze_regulondb_sigma_fimo --work-dir tmp/regulondb_sigma_fimo_reproduction --result-dir results/motif/regulondb_sigma_fimo_reproduction
.venv/bin/python -m experiments.validate_regulondb_sigma_holdout --work-dir tmp/regulondb_sigma_holdout_reproduction --result-dir results/motif/regulondb_sigma_holdout_reproduction
```

逐位点 FIMO 输出、FASTA 和序列 ID 留在本机 `tmp/`；Git 报告表不含逐序列 ID 或 DNA 序列。
