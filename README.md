# PR01-02 多物种启动子调控模式发现

研究细菌启动子的 motif 序列、相对 TSS 的位置及间距/排列，回答已知元件恢复、sigma 特异性和跨数据稳定性（RQ1–RQ3）。

## 2026-10-05 当前最终数据方案

当前路线以用户 2026-10-05 的明确决定为最高优先级：

> **TJU Pan / PromLoop 六物种为主数据；RegulonDB / DBTBS 主要负责 sigma annotation、机制解释和独立来源验证。**

在数据处理方案验收完成前，不运行新的正式 MEME/FIMO。embedding、Nucleotide Transformer、clustering、occlusion、SAE 等全部后置为探索项。

| 来源 | 当前角色 | 当前处理原则 |
|---|---|---|
| TJU Pan / PromLoop 六物种 | **RQ1/RQ3 主序列数据** | 云盘 `reg_and_gen/Datasets` 已与 PromLoop 固定提交逐文件核验一致；采用上游发布的 81 bp、相对 TSS `[-60,+20]` 标准化序列，不重新截 TSS、不自行按 strand 再反向互补 |
| RegulonDB E. coli K-12 | **RQ2 主 sigma 数据 + E. coli 跨数据集验证** | 官方 GFF3 `Sequence` 已按转录方向给出；不得再次按 genomic strand reverse-complement；保留完整 promoter↔sigma 多对多关系 |
| DBTBS 4.1 B. subtilis | **RQ2 补充 + B. subtilis 跨数据集验证** | 从原始网页重新解析 `Direction`、`Location`、`Absolute position`；只有证据完整且坐标一致的记录进入 TSS-aligned 子集，`Location=ND` 不进入依赖精确 TSS 的分析 |
| `strenth` 等课程扩展数据 | 可选后续实验 | 与主 motif 输入分开；测量定义、单位和来源未验收前不升级为核心证据 |

### TJU/PromLoop 的已确认边界

- 六物种 CSV 来自课程 TJU Pan，且与 `KevinHZS/PromLoop` 固定提交对应文件 30/30 SHA256 一致。
- PromLoop 发布说明定义序列为 81 bp、相对 TSS `[-60,+20]`，并说明做过 CD-HIT 90% 去冗余及 train/dev/test 划分。本项目把这些发布序列作为已经标准化的主输入。
- 不声称已经逐条独立验证 PromLoop 对每个原始负链记录的转换过程；项目层面采用其发布的数据定义。
- TJU CSV 没有逐条 sigma 标签，因此不作为 RQ2 sigma-specific 主证据。
- `Dataset.csv` / `positive_samples.csv` / train/dev/test 仍需按 ID、序列、标签、exact/RC 重复和近同源关系重新验收。原始文件不静默修改。
- label=0 只是一种背景候选。由于六物种正负样本 GC 存在系统差异，正式敏感性分析还需 GC-matched 与 dinucleotide-shuffled 背景。

## 已确认的数据问题与修复状态

### RegulonDB

旧处理曾对 **1,985 条已经按转录方向提供的负链序列再次 reverse-complement**，因此旧 3,807 条输入及直接依赖它的 motif/FIMO/position/spacing 结果不能作为正式证据。

当前 `preprocessing/source_adapters.py` 已按正确语义处理：GFF3 `Sequence` 直接作为转录方向序列，genomic strand 只保留为 provenance；TSS 由原始 GFF3 的显式标记读取，缺失/歧义不猜。新的 G0 重建同时保留 multi-sigma 信息并使用同源簇隔离划分。

旧输入和旧错误派生结果正在从正式结果目录清理；**原始 RegulonDB 下载快照不删除**。

### DBTBS 4.1

旧 parser 把 `tss_position` 和 `strand` 留空，但原始 DBTBS 页面实际含有：

- Genes 表的 `Direction`；
- promoter 的 `Location`；
- `Absolute position`；
- cis-element sequence 与实验 evidence。

DBTBS 官方 terms 明确：`Location` 是相对 transcription start site 的位置，`Absolute position` 按 NCBI accession `NC_000964` 重算，粗体标记 TSS。因此这些字段可以用于恢复大量记录的 strand/TSS，而不是把缺失当作数据库本身没有。

新增 `preprocessing/recover_dbtbs_coordinates.py`：只重解析冻结的原始 4.1 页面；只有 relative span、absolute span、序列长度和 gene direction 全部一致时才标记 `recovered_from_dbtbs`。`Location=ND`、方向混合/缺失、坐标不一致等记录保留 provenance，但不进入 TSS-dependent 核心分析。

在使用当前 NCBI `NC_000964.3` 重新截取 `[-60,+20]` 前，还必须用 DBTBS 页面序列/绝对坐标验证 accession version 的坐标兼容性，不能静默把无版本号的历史坐标当作当前版本已经验证。

## RQ 数据安排

- **RQ1：** TJU 六物种主 motif discovery；RegulonDB/DBTBS 的可靠注释用于已知元件和机制对照。
- **RQ2：** RegulonDB E. coli sigma-specific 为主；DBTBS 可恢复且样本量足够的 sigma 组作补充。TJU 不承担逐条 sigma 主证据。
- **RQ3：** 分开报告三类验证：
  1. TJU 内不同 species subset 的跨物种比较；
  2. 同来源按同源簇隔离的 discovery/holdout；
  3. 同物种不同来源：TJU E. coli ↔ RegulonDB、TJU B. subtilis ↔ DBTBS。跨来源前必须排除 exact、reverse-complement 和预定义阈值的近同源重叠。

“跨物种”和“跨数据集”不得混称；来源不同也不自动代表样本独立。

## 当前数据流水线

```text
raw snapshot
→ provenance / schema audit
→ source-specific normalization
→ TSS/orientation handling（TJU采用PromLoop发布标准化；外部库按原始证据恢复）
→ sequence/label QC
→ exact + reverse-complement dedup audit
→ near-homology clustering
→ sigma association normalization
→ discovery / holdout / cross-dataset isolation
→ background construction
→ FASTA + metadata + exclusion table + SHA256 manifest
→ 数据方案验收
→ 才开始正式 MEME/STREME/DREME/FIMO
```

所有正式派生数据记录来源版本、source record ID、排除原因、TSS/strand 证据、sigma 标签、cluster/split、输入输出 SHA256、处理命令和日志。

## 当前可用代码入口

- `preprocessing/source_adapters.py`：RegulonDB GFF3 正确转录方向适配。
- `preprocessing/recover_dbtbs_coordinates.py`：从冻结 DBTBS 4.1 页面恢复可验证的 strand/TSS/绝对坐标。
- `tests/test_source_adapters.py`：RegulonDB/旧 DBTBS 适配测试。
- `tests/test_dbtbs_coordinate_recovery.py`：DBTBS 坐标恢复边界测试。
- `experiments/prepare_tjupan_motif_inputs.py`：TJU/PromLoop 输入准备历史入口；正式版本仍需完成新的标签、重复、近同源和背景验收。

## 文档优先级

1. `AGENTS.md`：当前强制执行规则。
2. 本 README：当前数据角色和项目状态。
3. 2026-10-05 后续新增的数据审计/冻结说明。
4. `docs/07_后续实验执行规范_20261004.md`、`docs/08_全方位核查审计_20261004.md` 等保留其历史方法和审计证据；其中与 2026-10-05 数据角色冲突的旧表述已被当前方案取代，不能继续作为全局数据角色依据。

历史报告可保留用于追溯，但凡依赖已确认错误 RegulonDB 输入的结果均不是正式科学证据。原始下载数据始终保留，不因派生结果淘汰而删除。
