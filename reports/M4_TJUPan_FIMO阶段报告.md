# M4 阶段报告：TJU Pan 六物种 FIMO 首轮定位

更新日期：2026-10-03
项目：PR01-02 启动子调控模式发现
阶段状态：六物种 de novo FIMO、理想化已知参考矩阵扫描、motif pair 共现/最近间距首轮统计均已完成。de novo 有 353 个 site q≤0.05 位点；参考矩阵在相同 q 阈值下没有位点通过；540 项 motif pair 共现检验经 BH 后无显著项。结果仍属探索性，不是独立验证，也没有可确认的 TSS 有符号位置、稳定 spacing/order 或跨数据集结论。

## 1. 本轮做了什么

六物种 MEME 阶段各发现 10 个概率矩阵，本轮把 60 个矩阵分别用于对应物种，扫描该物种的正类启动子和长度合格的负类对照。合计扫描 9,724 条正类记录与 8,642 条对照记录。扫描前没有把序列反向互补或重排，所有输出位置均保留为输入 FASTA 内的 1–81 bp 坐标。

本轮目的是先建立每条序列上的 motif 命中、正负类唯一序列覆盖率和位置分布，供后续检验位置、间距和组合规律使用。FIMO 负责回答“这个矩阵在这些序列里哪些位置达到指定统计阈值”；它本身不能证明 motif 有调控功能。

## 2. 固定方法与统计口径

- 软件：MEME Suite FIMO 5.5.9。
- 输入：每物种自己的 10 个 MEME motif；合并扫描该物种正类与等长负类对照，FASTA 标题保留 label=1/0。
- 背景：每物种用其有效负类序列计算 A/C/G/T 零阶频率，作为 FIMO 背景；背景文件、计数和哈希见 `fimo_run_manifest.json`。
- 链：motif 文件声明 `strands: + -`，FIMO 扫描正向及反向互补链。`+/-` 表示相对于当前输入字符串的链，不证明字符串已按转录方向统一。
- 显著位点门槛：FIMO `--qv-thresh --thresh 0.05`，最多保存 1,000,000 个 score，实际未触及截断上限。
- q 值范围：核对本机 FIMO 5.5.9 的日志和实现后，确认它对**每个 motif 分别**计算 q 值，不是把同一物种的 10 个 motif 合并成一个 q 值集合；每个 motif 的运行日志显示它用 10,000 个随机抽取的 p 值估计 π₀。每个 motif 的 π₀ 记录在 run manifest。FIMO 的阈值和输出字段说明见 [FIMO 官方文档](https://meme-suite.org/meme/doc/fimo.html)。
- 正负类覆盖比较：把同一 motif 的 FIMO q≤0.05 命中按序列去重，构成正类命中/未命中与对照命中/未命中 2×2 表，计算双侧 Fisher 精确检验；再对 60 个“物种×motif”比较统一做 Benjamini–Hochberg 校正，报告 `fisher_bh_q_value_60_tests`。这项 BH 是本项目汇总脚本执行的，不是 FIMO 输出的 site q 值。
- 重复：命中覆盖率按独特 DNA 序列计算，不让正类内完全重复序列重复计数。六个物种的正负标签间没有完全相同的序列，故 Fisher 检验没有因跨标签重复而排除序列。MEME/FIMO site-level 扫描本身仍沿用源文件记录；少量正类重复存在，故覆盖率/组间比较以去重后的独特序列为准。
- 范围边界：MEME motif 是在同一批正类中发现，再在这批正类上进行 FIMO 定位。因此这不是独立测试集验证；FIMO q 值也不校正“先用同一数据挑选 motif”带来的选择偏差。所有正负类差异均先称为探索性关联。

## 3. 六物种结果概览

“显著位点”指 FIMO 按上述 per-motif q≤0.05 输出的 site；“有命中 motif”只表示至少有一个 site 通过，不等于每个 motif 都可信。

| 物种 | 正类记录 | 对照记录 | FIMO 显著位点数 | 至少有一个显著位点的 motif 数 | 正负覆盖差异经 60 项 BH 后 q≤0.05 的 motif 数 |
| --- | ---: | ---: | ---: | ---: | ---: |
| *Bacillus subtilis* | 691 | 809 | 126 | 2 | 2 |
| *Acinetobacter baumannii* | 1,540 | 1,360 | 0 | 0 | 0 |
| *Bradyrhizobium* sp. | 1,987 | 1,952 | 0 | 0 | 0 |
| *Corynebacterium diphtheriae* | 1,654 | 1,655 | 227 | 4 | 3 |
| *Escherichia coli* | 1,645 | 1,704 | 0 | 0 | 0 |
| *Staphylococcus* sp. | 2,207 | 1,162 | 0 | 0 | 0 |
| **合计** | **9,724** | **8,642** | **353** | — | **5 个物种×motif 项** |

四个物种没有任何位点通过当前 q≤0.05 门槛。这只表示在本轮背景模型、矩阵和显著性门槛下未检出过线位点；不等于序列里没有相似片段，也不等于这些物种没有启动子 motif。

## 4. 当前值得复核的探索性信号

| 物种 | motif / consensus | FIMO 显著位点 | 正类唯一序列覆盖 | 对照唯一序列覆盖 | 覆盖率差 | Fisher 双侧 p | 60 项 BH q |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| *B. subtilis* | MEME-2 `TTTTTYMWBWTDWWH` | 118 | 85/691（12.30%） | 19/809（2.35%） | +9.95 个百分点 | 1.64×10⁻¹⁴ | 3.28×10⁻¹³ |
| *B. subtilis* | MEME-9 `ATKTTHMTTYAAAA` | 8 | 8/691（1.16%） | 0/809（0%） | +1.16 个百分点 | 0.00198 | 0.0238 |
| *C. diphtheriae* | MEME-2 `KYCATGBBBNCCABBMTA` | 109 | 105/1,654（6.35%） | 4/1,655（0.24%） | +6.11 个百分点 | 1.84×10⁻²⁷ | 5.53×10⁻²⁶ |
| *C. diphtheriae* | MEME-3 `TCATRRYYAHHHABYSTA` | 103 | 102/1,654（6.17%） | 1/1,655（0.06%） | +6.11 个百分点 | 2.08×10⁻³⁰ | 1.25×10⁻²⁸ |
| *C. diphtheriae* | MEME-9 `ATKAGTWTSACAWT` | 10 | 10/1,654（0.60%） | 0/1,655（0%） | +0.60 个百分点 | 0.000960 | 0.0144 |

解释这些表时要特别留意三点：

1. MEME-9 在两个出现它的物种中，MEME 阶段 E-value 分别为 2.4 和 2.0，**并未达到**前文用于“候选观察”的 E-value≤0.05 线。因此即使这两个项的正负覆盖 Fisher q≤0.05，也只能列为待复核线索，不能说成已确认 motif。
2. FIMO site q 值和 Fisher BH q 值回答不同问题：前者判断单个 motif 位点是否越过 FIMO 阈值；后者比较该 motif 在正类/对照中的独特序列命中比例。它们都没有解决同一批数据先发现再扫描的选择偏差。
3. *C. diphtheriae* 的 MEME-2 与 MEME-3 很多命中起点紧邻或相同：在正类命中中，MEME-2 的起点 48 有 89 个 site，MEME-3 的起点 47 有 71 个 site。这提示两矩阵可能标记相邻/重叠的同一局部模式，不能把它们直接当成两个独立调控元件。位置仍是 FASTA 序列内坐标，不转成 TSS 坐标。

## 5. 结果文件与复现

报告级文件不包含逐条 DNA 序列或样本 ID：

- `results/motif/tjupan_m3_20261003/fimo/motif_hit_summary.tsv`：60 个物种×motif 项的显著位点、唯一序列覆盖和正负类差异统计。
- `results/motif/tjupan_m3_20261003/fimo/position_distribution.tsv`：按物种、motif、标签、输入序列起点和链汇总的位点数。
- `results/motif/tjupan_m3_20261003/fimo/fimo_summary_manifest.json`：统计口径、来源哈希、输出哈希与解释边界。
- `results/motif/tjupan_m3_20261003/fimo/fimo_run_manifest.json`：六次运行命令、软件版本、参数、样本数、负类背景频率、用时、返回码、pi0 估计和输出哈希。
- FIMO 原始 TSV/XML、合并 FASTA、背景文件及逐条日志仍位于 Git 忽略的 `tmp/`；原始 `fimo.tsv` 含序列名称和匹配序列，未作公开上传。

完整重跑命令：

```bash
.venv/bin/python -m experiments.run_tjupan_fimo
.venv/bin/python -m experiments.summarize_tjupan_fimo
```

第一次 Bacillus FIMO 尝试的自定义背景文件格式不符合工具要求，程序立即返回错误；修正为每行一项的背景频率格式后重跑成功。两次状态都留在 FIMO run manifest 的 `attempt_history`，正式结果仅使用成功运行。

## 6. 尚未完成与下一步

- 方向审计：云盘表无 strand 字段。当前坐标只写输入序列 1–81，不解释为 TSS 上下游。
- 已知元件真实位置对照第一轮已跑，site q≤0.05 下无参考位点；原始 p 值和位置敏感性已完成，GC 分层后 −10 窗口信号只在较宽松阈值下稳定，不能作为确认结果。详见[M4 motif 稳健性与跨数据验证](M4_Motif稳健性与跨数据验证.md)。
- motif 间关系：已完成第一轮共现、最近位点间距和输入字符串顺序描述；540 项共现 Fisher 检验均未通过 BH。只观察到两种小样本共现组合，尚无位置随机化，因此没有可靠的调控语法结论。
- 稳健性：参考 PWM 阈值敏感性与 GC 分层检验、二核苷酸随机背景外部复核已完成；MEME/FIMO 的背景阶数、低复杂度和训练/验证拆分敏感性仍待分析。
- 跨物种/数据集：六物种 60 个 PWM 比较和去重后的 TJU Pan *E. coli*/RegulonDB 复核已完成；跨物种全局 BH 无显著配对，独立数据的 FIMO site q≤0.05 也无显著位点。
- M3 中期汇报 PPT 暂停，当前优先推进实验；实验 6 第一轮组合分析和实验 5 第一轮强度关联已完成，后续推进两项的留出/敏感性复核、调控语法和自定义 RQ4。

本轮结果应作为 M4 第一轮定位，不改变数据分工：TJU Pan 六物种仍是主数据，RegulonDB 只承担 σ 分组和独立数据集辅助验证，DBTBS 仍是可选增强。

## 7. 已知参考矩阵 FIMO 扫描与坐标重叠（2026-10-03）

### 目的与方法

Tomtom 比较的是 motif 矩阵形状是否相似，不会告诉我们 motif 在每条 DNA 序列的哪个位置。因此另用 FIMO 扫描 `baselines/known_promoter_elements.meme` 中四个参考矩阵：大肠杆菌 σ70 −10、−35、UP proximal、UP distal。六物种分别使用与 de novo FIMO 完全相同的合并 FASTA 和本物种 label=0 零阶背景，扫描两条链，site q≤0.05；参考 motif 各自独立计算 q 值。坐标仍为输入序列 1–81。

随后把每个显著 de novo site 与同一条序列中的显著参考 site 做区间相交：坐标采用 FIMO 的 1-based inclusive 约定，共享至少一个碱基即记为 overlap；同时输出任意链重叠与输入链相同的重叠数。聚合表只含物种、motif 名、计数和比例，不含序列 ID、命中序列或完整逐位点输出。

### 结果

- 六个物种 × 四个参考 motif 共 24 个组合；FIMO site q≤0.05 的显著参考位点总数为 0。
- 24 项参考 motif 正类/负类唯一序列覆盖率 Fisher 检验，经 24 项 BH 校正后 q≤0.05 的比较数为 0。
- De novo FIMO 原有 353 个显著位点；由于参考位点集合为空，本阈值下坐标 overlap 数为 0。这个零值首先是“没有显著参考位点可供相交”的结果，不能写成“de novo motif 与已知元件没有任何碱基重合”，更不能写成“启动子没有 −10/−35/UP 元件”。
- 12 份 FIMO/CisML XML 均通过 `xmllint` 解析。

结果文件：

- `results/motif/tjupan_m3_20261003/reference_elements/known_reference_coverage_association.tsv`
- `results/motif/tjupan_m3_20261003/reference_elements/denovo_known_reference_overlap.tsv`
- `results/motif/tjupan_m3_20261003/reference_elements/reference_fimo_run_manifest.json`
- `results/motif/tjupan_m3_20261003/reference_elements/reference_overlap_summary_manifest.json`

运行脚本为 `experiments/run_tjupan_reference_fimo.py` 和 `experiments/summarize_tjupan_reference_overlap.py`；原始逐位点结果和 FIMO 日志只保存在本机忽略目录 `tmp/fimo_runs/tjupan_promloop_20261003/reference_elements/`。

### 解释边界

四个参考是少数理想化的大肠杆菌 σ70 模板，并非每个物种、每个 σ 因子的完整 PWM 文库；在非大肠杆菌上的命中只能作为探索性参考。严格 site-q 扫描没有参考位点通过，不能据此推断真实调控元件不存在。新增原始 p 阈值、序列坐标和 GC 分层敏感性结果见[M4 motif 稳健性与跨数据验证](M4_Motif稳健性与跨数据验证.md)。云盘数据方向仍未核实，输入坐标不能直接解释为转录方向上的“−10”“−35”。

方向证据审计检查了上游公开说明：[PromLoop README](https://github.com/KevinHZS/PromLoop#-datasets) 说明 81 bp 窗口相对 TSS 为 −60 至 +20，但没有说明反向链是否已经统一反向互补；当前云盘 CSV 字段也没有 strand 或基因组 TSS 坐标。因此方向结论是“未确认”，不是已证实统一或随机。后续若要把序列内坐标转换成上游/下游坐标，必须取得能追溯到每条序列的链方向或上游生成流程证据。

## 8. Motif 共现、最近位点间距与顺序首轮统计（2026-10-03）

### 统计方法

本轮只使用六物种 de novo FIMO q≤0.05 位点。先在每个物种、正类/对照内按 DNA 完全相同序列去重，再为每个 motif pair 建立“是否命中 A × 是否命中 B”的 2×2 表，做双侧 Fisher 精确检验。六物种、两类标签、每物种 45 个 pair，共 540 项检验，统一做 BH 校正。

对同时命中两个 motif 的唯一序列，分别取两组位点中间隔最小的一对，报告它们在输入字符串上的左/右次序及 signed gap：不重叠时为两区间间未匹配碱基数；重叠时为共享碱基数的负值。每条唯一序列只贡献一个最近位点对，避免多重 hit 序列把同一 pair 过度计数。此间距是描述统计，本轮没有做位置随机化检验。

### 结果

- 540 项 motif-pair 共现检验经 BH 校正后，q≤0.05 项数为 0。
- 非零正类共现只有两对：

| 物种 | motif pair | 共现唯一序列 | 最近位点结果 | 统计解释 |
| --- | --- | ---: | --- | --- |
| *B. subtilis* | MEME-2 / MEME-9 | 4 | 中位间隔 18 bp；1 条序列为 A 在 B 左侧、3 条相反 | Fisher 原始 p≈0.0101，但 540 项 BH q=1；探索线索 |
| *C. diphtheriae* | MEME-3 / MEME-9 | 2 | 两条最近位点对均重叠 4 bp | Fisher 原始 p≈0.1227，BH q=1；样本极少 |

- 其余 species×label×motif-pair 组合没有同一条唯一序列同时命中两者，故没有可计算的最近间距。
- 这些结果支持“当前显著位点集合里几乎观察不到稳健共现”，不支持“细菌启动子不存在 motif 组合”。非显著可能来自 hit 数少、q 阈值严格、正/负标签构造、方向未知或候选 motif 尚未独立验证。

结果文件：

- `results/motif/tjupan_m3_20261003/pairs/motif_pair_cooccurrence.tsv`：540 项物种×标签×motif pair 汇总，含 BH 校正值。
- `results/motif/tjupan_m3_20261003/pairs/motif_pair_spacing.tsv`：最近位点间距、重叠和输入字符串顺序描述。
- `results/motif/tjupan_m3_20261003/pairs/motif_pair_manifest.json`：样本去重口径、检验范围、哈希和限制。
- 复现命令：`.venv/bin/python -m experiments.analyze_tjupan_motif_pairs`

脚本为 `experiments/analyze_tjupan_motif_pairs.py`。新增参考 PWM、跨物种和跨数据集复核结果见[M4 motif 稳健性与跨数据验证](M4_Motif稳健性与跨数据验证.md)。

## 9. 实验 6：匹配随机背景下的 motif 组合

前述 540 项检验只纳入 FIMO q≤0.05 的严格位点，因共现极少而检验能力有限。后续实验 6 对六物种 MEME 矩阵分别使用原始 site p≤0.05、0.01、0.001，并为每条正类序列生成三份保持二核苷酸频数的配对随机背景，比较 motif 对共现、最近间距和输入字符串次序。

三种子最差 p 值经 270 项 BH 后，稳健 pair 数为：p≤0.05 时 0/270，p≤0.01 时 176/270，p≤0.001 时 80/270；没有 pair 在三档阈值下都稳健。部分共现高于二核苷酸随机背景，但结果依赖位点阈值，且 motif 来自同批正类序列，尚不能解释为已验证的调控组合。方向和独立验证问题仍未解决。

分析方法、物种分布、间距描述与限制见[实验 6 报告](实验6_Motif组合与背景对照.md)；机器可读结果位于 `results/motif/tjupan_m3_20261003/combinations/`。
