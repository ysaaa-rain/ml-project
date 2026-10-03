# M3 阶段报告：TJU Pan 六物种 motif 发现

更新日期：2026-10-03
项目：PR01-02 启动子调控模式发现
阶段状态：六物种全量 MEME 搜索、已知元件初步矩阵比较和报告级 motif 可视化均已完成；FIMO 定位及 M4 后续实验未完成。

说明：本文记录 M3 阶段收口时的结果快照；之后的 FIMO 首轮扫描已另行完成，当前状态和结果见[六物种 M4 FIMO 阶段报告](M4_TJUPan_FIMO阶段报告.md)。

## 1. 本阶段目标与边界

M3 要从启动子 DNA 序列中找出反复出现的短序列模式（motif），分别在六个物种中进行相同设置的全体发现，然后把发现的 motif 矩阵与已知启动子元件参考矩阵作第一轮相似性比较。本阶段回答的是“算法在这些正类序列和对照下提出了哪些候选模式”，不是“这些模式已被证实有调控功能”。

本报告只记录 TJU Pan / PromLoop 云盘 `reg_and_gen/Datasets` 六物种主数据。另有 RegulonDB 辅助 σ 分组结果，记录在[RegulonDB 辅助 M3 报告](M3阶段报告.md)；它不计入本文六物种 motif 数量。这里的 Tomtom 比较是矩阵相似性检查，不等于 FIMO 在启动子上的定位。

## 2. 数据与输入核验

每个物种使用 `positive_samples.csv` 作为 MEME primary 序列；它已逐 ID、逐序列核对为 `Dataset.csv` 中 `label=1` 的完整子集。`Dataset.csv` 里的 `label=0` 记录作 differential-enrichment 对照，只保留长度为 81 bp 的记录。输入准备脚本会检查字符集、长度、ID 唯一性、拆分交叠和文件 SHA256；原始序列 FASTA 保存在 Git 忽略的本地 `tmp/`，manifest 不保存 DNA 字符串。

| 物种 | 正类 primary 数 | 等长负类对照数 | 发现 motif 数 | MEME 报告 E-value ≤ 0.05 的候选数 | 最小 E-value 候选（consensus；宽度；sites；E-value） | 成功运行耗时 |
| --- | ---: | ---: | ---: | ---: | --- | ---: |
| *Bacillus subtilis* | 691 | 809 | 10 | 5 | `AGGAGG`；6 bp；161；2.7e-21 | 175.32 秒 |
| *Acinetobacter baumannii* | 1,540 | 1,360 | 10 | 10 | `TGCCVY`；6 bp；297；7.4e-34 | 868.53 秒 |
| *Bradyrhizobium* sp. | 1,987 | 1,952 | 10 | 7 | `GGCGCG`；6 bp；451；3.5e-60 | 1,512.84 秒 |
| *Corynebacterium diphtheriae* | 1,654 | 1,655 | 10 | 6 | `GGGGTG`；6 bp；264；5.6e-30 | 1,041.50 秒 |
| *Escherichia coli* | 1,645 | 1,704 | 10 | 9 | `CRSCAR`；6 bp；599；7.5e-84 | 1,065.75 秒 |
| *Staphylococcus* sp. | 2,207 | 1,162 | 10 | 10 | `KTYRMC`；6 bp；721；1.9e-64 | 1,417.72 秒 |
| **合计** | **9,724** | **8,642** | **60** | **47** | — | **6,081.66 秒（约 101.4 分钟）** |

对照比正类少不是另行抽样造成：原始数据本来各物种正负比例不同；脚本只剔除了四条不符合 81 bp 窗口的负类。输入数量、来源和哈希见 [`input_manifest.json`](../results/motif/tjupan_m3_20261003/input_manifest.json)。

## 3. MEME 方法与参数

软件为 MEME Suite 5.5.9。六物种使用同一参数：DNA 字母表；`zoops`（每条 primary 序列最多一个零次或一次 motif 位点的模型）；`de` differential enrichment 目标；各物种自己的 81 bp 负类对照；最多 10 个 motif；宽度 6–20 bp；搜索正反向互补链；随机种子 20261003；`-searchsize 0`；`-brief 50`。详细固定配置见 [`m3_tjupan_meme_20261003.yaml`](../configs/m3_tjupan_meme_20261003.yaml)，实际命令、版本、时间、样本数、返回码、输入输出 SHA256 见 [`run_manifest.json`](../results/motif/tjupan_m3_20261003/run_manifest.json)。参数语义以 [MEME 官方说明](https://meme-suite.org/meme/doc/meme.html) 为准：`searchsize=0` 禁用序列抽样并使用完整 primary 序列集合；`de` 目标自己的留出机制仍适用。

Staphylococcus 第一次搜索在 2026-10-03 运行时意外中断，只留下 3 个 motif 的部分文本且没有 XML，未被算作成功。部分输出和日志被保留在本地 Git 忽略目录 `tmp/meme_runs/tjupan_promloop_20261003/interrupted_attempts/staphylococcus_20261003_155433/`；同一输入和参数重跑后以退出码 0 完成。两次尝试都在 run manifest 的 `attempt_history` 中可追踪，正式结果只采用第二次完整运行。

## 4. motif 结果怎样解读

全体 60 个发现矩阵均可从 XML 解析，六组各有 10 个 motif。按各次 MEME 输出中的 E-value≤0.05 作为“候选观察线”，六组分别得到 5、10、7、6、9、10 个候选，共 47 个。这个筛线**没有对六个物种的 60 个 motif 再做跨组多重比较校正**，所以不能把 47 个候选说成 47 个已验证的统计发现。

候选模式在不同物种间不完全相同；例如 *E. coli* 中 `CRSCAR`（E=7.5e-84），*Bradyrhizobium* 中 `GGCGCG`（E=3.5e-60），以及 Staphylococcus 中 `KTYRMC`（E=1.9e-64）在本次各自数据/搜索条件下排位靠前。60 条 consensus 字符串没有完全相同的重复项，但这只是字符串完全相等的描述统计：IUPAC 模糊码、反向互补和概率矩阵间的近似关系都还没有做正式跨物种检验。

表里的 `sites` 是 MEME 估计/报告的位点数，不等于“有多少条不同启动子命中”；同一序列贡献、序列覆盖率和准确坐标需在 FIMO 阶段另算。`-brief 50` 限制逐序列贡献位置的报告细节，故当前 manifest 和摘要不能代替完整逐序列定位结果。

特别需要人工复核的例子是 *B. subtilis* 的 `AGGAGG`。它常见于细菌翻译起始附近的 Shine–Dalgarno 样序列；由于当前 81 bp 窗口延伸到 TSS 下游 +20 bp，这个发现可能来自转录起始之后、与翻译相关的局部序列，而不一定是启动子调控元件。这里只把它标为需要结合位置验证的候选，不作功能定论。

## 5. 与 −10/−35/UP 参考矩阵初步比较

使用四个理想化参考矩阵：σ70 `−10 TATAAT`、σ70 `−35 TTGACA`、UP 近端 `AAAAAARNR`、UP 远端 `AWWWWWTTTTT`；Tomtom 参数为 Euclidean distance、最小重叠 6 nt、允许参考反向互补方向、阈值 1（保留所有配对）。参考文件 [`known_promoter_elements.meme`](../baselines/known_promoter_elements.meme) 新增自动共识断言，确保 −10/−35 实际矩阵分别仍为 `TATAAT` / `TTGACA`。

六物种共 60 个 de novo motif × 4 个参考，产生 240 条 query-reference 配对。Tomtom 明确警告 4 个 target motif 太少、无法准确估计 p-value；每个 query 只有 8 个方向感知比较，q-value 也不稳定。因此表中的最佳参考仅是“四个模板里排在最前者”，p/q 值不作为显著性检验或确认性证据。就 47 个 E-value 候选而言，最靠前的参考模板计数为：

| 物种 | 最佳模板为 σ70 −10 | 最佳模板为 σ70 −35 | 最佳模板为 UP 近端 | 最佳模板为 UP 远端 |
| --- | ---: | ---: | ---: | ---: |
| *B. subtilis* | 0 | 0 | 5 | 0 |
| *A. baumannii* | 3 | 3 | 4 | 0 |
| *Bradyrhizobium* sp. | 0 | 1 | 5 | 1 |
| *C. diphtheriae* | 2 | 2 | 2 | 0 |
| *E. coli* | 2 | 5 | 2 | 0 |
| *Staphylococcus* sp. | 1 | 2 | 7 | 0 |

这里许多配对只重叠最小允许的 6 nt；模板来源和长度本身也不均衡（两个 σ70 核心框各 6 nt，UP 模板较长），所以最佳类别数不能直接解释成“哪个元件在物种中更多”。Tomtom 不检查 motif 是否落在各条 promoter 的 −10、−35 或 UP 基因组坐标上，也不检查共现、转录活性、功能或 σ 因子特异性。真实位置重叠率必须等 FIMO 命中、方向/坐标口径审计通过后再计算。

## 6. 可复核的报告级结果

正式整理结果目前包含概率矩阵、summary 表、已知元件配对和可读 logo，不含逐条 DNA 序列：

- [`summary/motif_summary.tsv`](../results/motif/tjupan_m3_20261003/summary/motif_summary.tsv)：60 个 motif 的宽度、consensus、E-value、sites、信息量及 GC 概要。
- [`summary/matrices/`](../results/motif/tjupan_m3_20261003/summary/matrices/)：六个物种的 MEME 概率矩阵文件，可作为 FIMO/Tomtom 输入。
- [`logos/six_species_motif_logos.pdf`](../results/motif/tjupan_m3_20261003/logos/six_species_motif_logos.pdf) 及各物种 PNG：60 个 motif 的统一信息量 logo；字母高度按 `2 − Shannon entropy` 计算，未作有限样本修正。
- [`comparison/all_species_reference_matches.tsv`](../results/motif/tjupan_m3_20261003/comparison/all_species_reference_matches.tsv)：240 条完整 Tomtom 配对。
- [`comparison/best_reference_match_per_motif.tsv`](../results/motif/tjupan_m3_20261003/comparison/best_reference_match_per_motif.tsv)：60 个 motif 的最佳模板排序摘要。
- 三份 `*_manifest.json`：阶段完整性、软件/参数、输入输出 SHA256 与统计/生物学限制。
- 原生 XML/HTML、FASTA、Tomtom HTML/XML、运行日志和中断尝试归档仍位于 Git 忽略的 `tmp/`，便于本机复核，不直接提交。

## 7. M3 验收状态与下一步

| 检查项 | 状态 | 证据/说明 |
| --- | --- | --- |
| 六物种 primary/control 输入验证 | 已完成 | input manifest、81 bp 长度和输入 hash |
| 六物种全量 MEME 搜索 | 已完成 | 每组退出码 0、10 个 motif、XML 可解析 |
| 矩阵/summary/logo 整理 | 已完成 | 60 个 motif、6 张物种 PNG + 1 份 PDF |
| 与 −10/−35/UP 模板初步比较 | 已完成（描述性） | 240 对；模板小导致 p/q 不可靠 |
| FIMO 定位及唯一启动子覆盖率 | 未完成 | 需要锁定 motif、阈值及多重校正策略 |
| TSS 坐标解释 | 未完成 | CSV 无链方向字段，需先做方向审计 |
| motif 间距、顺序、组合与跨数据集复核 | 未完成 | 属于后续 M4/实验 5、6 等工作 |

下一步先做 FIMO 方案锁定及 81 bp 序列方向核查，随后扫描并分别报告 motif 的命中条数、唯一启动子数、序列内位置分布和多重比较校正值；方向证据不足时，位置只写成输入 FASTA 内的 1–81 bp 坐标，不换算成 TSS 上下游。之后再做 spacing/order、跨物种 PWM 比较及 RegulonDB 去重后的独立数据集验证。云盘 CSV 没有 σ 标签，不能仅凭这六物种回答 σ 组差异。

M3 搜索产物复现命令见[运行说明](../docs/运行说明.md)；当前仓库测试为 86 项全部通过（16 条 matplotlib/pyparsing 第三方弃用警告）。数据及其衍生结果公开再分发许可仍需最终核对；在许可确认前，不把原始序列或完整序列可恢复产物推送到 GitHub。
