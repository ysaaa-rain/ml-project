# TJU Pan 六物种 M2 阶段报告：描述性数据探索

## 当前采用规则（2026-10-05）

更新：2026-10-05。当前数据版本：`data/processed/pr01_02_data_v2/`。

TJU 六物种为主数据；RegulonDB / DBTBS 主要负责 sigma annotation、sigma-specific motif 与隔离后的不同来源验证。

本研究采用 PromLoop 发布的标准化 promoter 序列，并按其数据集定义视为已经完成 TSS 对齐及必要的方向标准化。81 bp 对应相对 TSS 的 [-60,+20]，TSS 为局部第61位；TJU 不重新截取或反向互补。genomic coordinate、逐条 strand 和原始数据库映射不是 TJU 主实验必要输入。

本轮仅修复和构建数据，没有运行正式 MEME/STREME/DREME/FIMO。旧 RegulonDB 3,807 条处理输入及其实验结果已从项目移除，不能继续引用旧覆盖率、σ显著性或位置结论。

以下数值如属于旧TJU全量探索，保留为历史结果，不作为v2冻结验证结论。旧RegulonDB相关结果已撤回，链接若指向旧输入不再有效。当前构建/状态见《数据修复与输入验收_20261005》。


更新日期：2026-10-03
分析范围：TJU Pan 云盘 reg_and_gen/Datasets 的六物种完整数据。本报告不作 motif 发现或功能显著性结论。

## 1. 目的和口径

本阶段检查每个物种正负样本规模、序列长度、GC 含量、重复情况、拆分间完全相同的序列，以及正类各序列位置的碱基组成。由于 CSV 没有 strand/TSS 逐条字段，位置统计只按 1–81 序列索引解释，不能写成相对 TSS 的上游/下游位置。

每个物种的 Dataset.csv 是完整表；train/dev/test 是该表的拆分。positive_samples.csv 的正类与 Dataset.csv 中 label=1 已按记录 ID 和序列核对。详细输入文件哈希见 results/eda/tjupan_m2_20261003/run_manifest.json。

## 2. 样本数、长度与重复

六物种共 18,370 条记录：正类 9,724、负类 8,646。所有正类均为 81 bp A/C/G/T。负类各有一条 80 bp 记录的物种共四个；MEME 的等长负类对照中将这四条排除，剩余 8,642 条。

| 物种 | 正类记录 | 负类记录 | 正类唯一序列 | 负类长度范围 |
| --- | ---: | ---: | ---: | --- |
| B. subtilis | 691 | 809 | 691 | 81 |
| A. baumannii | 1,540 | 1,360 | 1,540 | 81 |
| Bradyrhizobium | 1,987 | 1,953 | 1,983 | 80–81 |
| C. diphtheriae | 1,654 | 1,656 | 1,654 | 80–81 |
| E. coli | 1,645 | 1,705 | 1,637 | 80–81 |
| Staphylococcus | 2,207 | 1,163 | 2,206 | 80–81 |

表中负类数为原始计数；四条长度 80 bp 的记录不进入等长 MEME 对照。正类重复记录合计 13 条，保留在原始输入；后续 motif 唯一序列覆盖率不能直接拿命中数作分子。拆分间精确相同序列在 Bradyrhizobium 发现 2 条、E. coli 发现 3 条，其余物种未检出；这项检查针对拆分泄漏，不意味着各拆分是独立数据集。

## 3. GC 组成观察

| 物种 | 正类平均 GC | 负类平均 GC | 正类减负类 |
| --- | ---: | ---: | ---: |
| B. subtilis | 0.3287 | 0.4606 | −0.1319 |
| A. baumannii | 0.3314 | 0.4129 | −0.0815 |
| Bradyrhizobium | 0.6057 | 0.6501 | −0.0444 |
| C. diphtheriae | 0.4663 | 0.5511 | −0.0848 |
| E. coli | 0.4416 | 0.5297 | −0.0881 |
| Staphylococcus | 0.2563 | 0.3425 | −0.0862 |

六个物种的正类平均 GC 都低于对应负类，差值约为 4.4–13.2 个百分点。这是重要的背景组成差异：即使 MEME 以真实负类作对照，候选序列模式也要检查是否主要由 GC/A-T 组成差异推动；不能只凭 E-value 或 logo 形状解释成调控元件。后续可用组成匹配/二核苷酸背景作稳健性对照。

## 4. 序列位置组成

每个物种正类均按 81 个序列索引计算 A/C/G/T 频率、Shannon 熵和相对均匀四碱基背景的 2−H 信息量。该值是未校正有限样本偏差的描述性统计，不是 MEME 显著性、进化保守性或功能证据。最大逐位点 2−H 在各物种中约为：B. subtilis 位置 50、0.405 bits；A. baumannii 位置 32、0.328；Bradyrhizobium 位置 71、0.094；C. diphtheriae 位置 61、0.605；E. coli 位置 61、0.151；Staphylococcus 位置 32、0.439。峰值位置和高度只是待复核的组成描述，不代表已发现 motif。逐物种结果和图见结果目录。

由于六物种的链方向尚未从 CSV 得到独立证明，图表横轴为序列位置 1–81。方向审计通过后，才考虑将这些位置换算到共同的 TSS 坐标。

## 5. 复现与产物

运行命令：

~~~bash
python -m experiments.eda_tjupan
~~~

脚本默认读取 data/raw/tju_pan_promoter/reg_and_gen/Datasets/，输出到 results/eda/tjupan_m2_20261003/。结果包括：

- species_class_summary.tsv：按物种和正负类的数量、唯一序列、长度及 GC 摘要；
- sequence_type_counts.tsv：记录类型计数；
- positive_position_composition.tsv：正类逐位置碱基比例和信息量；
- split_overlap_summary.tsv：train/dev/test 完全相同序列重叠数；
- 三张汇总图；
- run_manifest.json：输入文件及产物 SHA256、样本数、坐标限制和统计解释。

当前结果包不含逐条序列，也不含样本 ID。云盘数据及其衍生结果的公开许可仍需核实，确认前不推送可能受限制的数据产物。

## 6. M2 判定及下一步

六物种的样本规模、GC/长度、重复和序列内逐位点组成等描述性 EDA 已完成，可作为 M2 主数据基线。仍未完成的是链方向证据、TSS 相对位置解释以及以已知元素/背景模型作适当稳健性比较。接下来先完成六物种 MEME 并验收结果，再把高 GC/A-T 组成差异纳入 motif 解释和对照设计。
