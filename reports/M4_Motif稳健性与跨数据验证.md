# M4 motif 稳健性与跨数据验证

更新日期：2026-10-04

范围：TJU Pan 六物种 PWM 阈值与位置复核、六物种矩阵比较、E. coli 去重跨数据集复核。

## 1. 参考 PWM 阈值与位置

四个理想化 E. coli 参考 PWM 使用六物种各自的 label=0 零阶背景，报告单个位点原始 p≤0.05、0.01、0.001。总原始命中数依次为 585,721、131,123、15,112。按独特序列比较正类与对照、每档 24 项 Fisher 检验分别做 BH 后，q≤0.05 项依次为 18、23、22 项；这是阈值化命中比例差异，不是 FIMO 位点显著性。

按上游所述 −60 至 +20 窗口假设，σ70 −10 参考矩阵的预设起点为 49–54，−35 为 24–29。对 E. coli 两个窗口，在每个原始 p 阈值、两个输入链上按 81 bp GC 数分层置换标签 10,000 次，再对 12 项检验做 BH：

| 参考矩阵 | 位点 p 阈值 | 链 | 正类窗口命中 | 对照窗口命中 | GC 分层置换 BH q |
| --- | ---: | --- | ---: | ---: | ---: |
| −10 | 0.01 | + | 252/1,637 | 21/1,704 | 0.0003 |
| −10 | 0.01 | − | 154/1,637 | 27/1,704 | 0.0003 |
| −35 | 0.01 | + | 157/1,637 | 79/1,704 | 0.2925 |
| −35 | 0.01 | − | 114/1,637 | 57/1,704 | 0.0511 |
| −10 | 0.001 | + | 19/1,637 | 1/1,704 | 0.1483 |
| −10 | 0.001 | − | 6/1,637 | 1/1,704 | 0.6280 |
| −35 | 0.001 | + | 14/1,637 | 1/1,704 | 0.2178 |
| −35 | 0.001 | − | 4/1,637 | 2/1,704 | 0.6280 |

−10 窗口在原始 p≤0.01 时高于 GC 分层随机背景；收紧到 p≤0.001 后不再通过校正。−35 窗口在两档均未通过 12 项校正。两条链都出现命中；TJU Pan 输入方向仍未确认，因此不能把窗口命中表述为转录方向上的元件定位。四个参考 PWM 在原有 FIMO site q≤0.05 扫描中仍均无显著位点。

结果：`results/motif/tjupan_m3_20261003/reference_elements/pvalue_sensitivity.tsv`、`pvalue_position_distribution.tsv`、`ecoli_sigma70_expected_window_tests.tsv`、`pvalue_sensitivity_manifest.json`。复现：

```bash
.venv/bin/python -m experiments.analyze_reference_pwm_sensitivity
```

## 2. 六物种 motif 矩阵一致性

将六物种 60 个 de novo PWM 合并后，用 Tomtom 计算 3,600 个有方向的矩阵配对，其中 3,000 项为跨物种配对。每个 query 自身的 60 个目标比较中，有 18 个跨物种配对达到 q≤0.05，涉及 15 个无序物种对中的 7 对。对全部 3,000 项跨物种 p 值统一做 BH 后，没有配对达到 q≤0.05，最小全局 q=0.0644。

候选中，*B. subtilis* MEME-4 与 *C. diphtheriae* MEME-2/3 有 13 bp 相似重叠区；全局校正后仍不显著。其余候选多为 6 bp 最小重叠。当前结果不支持六物种间存在统计稳健的共同 motif；矩阵相似不等于功能相同。

结果：`results/motif/tjupan_m3_20261003/cross_species/` 下的 `species_pair_similarity.tsv`、`significant_cross_species_motif_matches.tsv`、`cross_species_manifest.json`。显著配对表按全局 BH q≤0.05 定义，本轮为空；每-query 候选和全局校正结果都保存在汇总表与 manifest。复现：

```bash
.venv/bin/python -m experiments.compare_tjupan_species_motifs
```

比较采用 pooled label=0 单核苷酸背景、Ed 距离、最小重叠 6 bp，并启用目标矩阵反向互补比较。Tomtom q 值用于每个 query 的目标库筛选；跨 query 的结论以 3,000 项全局 BH 为准。

## 3. 去重后的 RegulonDB 独立复核

RegulonDB 清洗集含 3,807 条独特的 81 bp、TSS 位于第 61 位的 E. coli 序列。与 TJU Pan E. coli 的 1,637 条独特正类发现序列比较后，排除完全同向重复 375 条、完全反向互补 369 条、无缺口且身份率至少 90%（73/81）的近同源 2 条；保留 3,061 条独立序列。

每条保留序列生成 5 份独立的二核苷酸保持随机序列，随机种子为 20261003–20261007。每轮随机序列与其原始序列一一配对，且每轮的 3,061 条正序列和 3,061 条随机序列均无重复或跨标签共享。FIMO 扫描 E. coli 的 10 个 M3 motif；主要检验是配对精确 McNemar 检验，每轮、每个命中阈值下对 10 个 motif 做 BH；另输出非配对 Fisher 结果供对照。结果：

- FIMO site q≤0.05：5 轮中 10 个 motif 均无显著位点。
- 原始 site p≤0.05：真实序列与随机序列覆盖率接近，5 轮均无 motif 通过 paired McNemar BH。
- 原始 site p≤0.01：MEME-1、2、4、5、6、7、9 在 5/5 轮通过 paired McNemar BH；MEME-3 为 4/5 轮，MEME-8 为 3/5 轮。前 7 项真实序列覆盖率中位数比随机序列高 3.2–9.1 个百分点。
- 原始 site p≤0.001：MEME-1、4、5、7、9 在 5/5 轮通过 paired McNemar BH，覆盖率中位数高 2.7–11.8 个百分点；MEME-6、8 各为 3/5 轮。

独立数据中可重复观察到若干 PWM 分数高于二核苷酸随机背景的模式，但这些位点没有通过 FIMO site q≤0.05。当前证据支持阈值依赖的序列模式重现，不确认严格位点显著性、调控功能或 σ 特异性。

过滤计数与序列级汇总：`results/motif/tjupan_m3_20261003/cross_dataset/dataset_overlap_summary.tsv`、`regulondb_external_motif_validation.tsv`、`regulondb_validation_manifest.json`。逐条序列、随机序列和 FIMO 原始位点留在本机 `tmp/`。复现：

```bash
.venv/bin/python -m experiments.validate_regulondb_cross_dataset
```

## 4. 研究问题的阶段性回答

- **RQ1，发现的 motif 是否对应已知元件？** 严格 FIMO site q≤0.05 下没有参考 PWM 位点通过。按预设窗口和 GC 分层置换，E. coli −10 窗口只在原始 site p≤0.01 时高于随机背景，收紧到 p≤0.001 后不再显著；输入方向未确认，因此目前是位置假设下的候选信号，不能确认元件身份或转录方向位置。
- **RQ2，不同 σ 因子的 motif、频率和间距是否不同？** RegulonDB 六组矩阵交叉扫描后，14 个 E-value≤0.05 候选中，13/14 个在 raw site p≤0.01、14/14 个在 p≤0.001 时高于配对 N1 背景；10/14 与 11/14 个在来源组内的命中率高于其他五组。组合检验在两档阈值分别有 129/270 和 134/270 个 pair 高于匹配背景，但 p≤0.05 时为 0。候选命中频率显示 σ 组相关性，spacing 结论依赖阈值，尚不能确认稳定 σ 特异语法。详见[RQ2 六 σ 组 FIMO 与间距分析](RQ2_RegulonDB_sigma_FIMO与间距分析.md)。
- **RQ3，规律能否跨物种或数据集重现？** 六物种矩阵比较的 3,000 项跨物种检验经全局 BH 后没有显著配对。去重后的 RegulonDB E. coli 序列在较宽松原始位点阈值下复现了部分 motif 覆盖差异，但严格 site q≤0.05 下无显著位点。当前是阈值依赖的局部重现，不支持普遍规律。
- **实验 5，motif 与强度标签是否相关？** 去除与 TJU Pan E. coli 发现/对照序列片段重叠后，11,137 条序列中有 5 个 E-value≤0.05 的 motif 在 GC 校正后与 dRNA-seq 派生数值标签呈弱负相关（|ρ|=0.028–0.049，BH q<0.05）。数据没有条件、单位和重复信息，且测量独立性未确认；这不是因果或强度预测结论。
- **实验 6，motif 组合是否高于序列背景？** 在原始 site p≤0.01/0.001 下，分别有 176/80 个 motif 对高于二核苷酸保持背景；没有 motif 对在三档阈值下都稳健，不能据此确认稳定调控语法。

## 5. 下一步

追溯 dRNA-seq 标签的条件、数值定义和数据来源关系；补做 de novo motif 留出集、背景阶数和低复杂度敏感性；对 RegulonDB σ 组候选 motif 频率和 spacing/order 做留出复核；取得云盘序列方向证据后再检验 TJU Pan TSS 坐标下的位置语法。

实验 6 结果、逐项统计口径和复现方法见[实验 6：motif 组合与背景对照](实验6_Motif组合与背景对照.md)。
