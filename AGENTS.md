# Repository instructions

## 后续实验路线（用户指定，2026-10-08）

现行文档入口为 `docs/README.md`。M1/M2 是课程要求的数据采集、按 σ 分组和只读统计阶段；不运行 de novo motif 发现。本轮甲乙丙旧稿已收敛为 `docs/M1M2_阶段调查与描述性发现_20261008.md`，并用 `docs/M1M2_待补只读分析与验收_20261008.md` 标注未完成的三来源统一统计。

在启动 PR01-02 新 STREME/FIMO/已知元件对照之前，必须阅读 `docs/PR01-02_motif_reference_and_grammar_protocol_20261008.md`。该文件是 2026-10-08 最新预注册规则：TJU 主 STREME 使用发布的 `label=0` 天然对照；另用已生成的 `discovery.dinucleotide_null.fasta` 进行必须执行的 STREME 稳健性分析。两套背景独立运行、禁止混合；不新增 GC-matched 背景；背景组成差异作为解释限制保留。DBTBS 不再主张有效规模的独立外部验证，保留 sigma/TSS/位点注释及明确标记的探索工作。真实参考 PWM 库须按该文件整理验收后才可宣称“库已冻结”，不能把理想 consensus 当实测 PWM。


开始任何实验设计、数据处理、运行、结果解释或报告前，阅读 `docs/README.md`、`docs/07_后续实验执行规范_20261004.md`、`docs/PR01-02_motif_reference_and_grammar_protocol_20261008.md` 和当前数据审计/交接文档。后续实验以用户最新决定为准；与旧路线冲突时，以本文件及最新数据规范为准。

PR01-02 当前数据架构固定为：**TJU Pan / PromLoop 六物种是主序列数据源；RegulonDB 用于 E. coli σ 专项与可隔离的外部验证；DBTBS 用于 B. subtilis σ/TSS/位置和文献注释，不作独立外部统计验证。** 不再采用“RegulonDB/DBTBS 为主发现、TJU 仅为辅助探索”的旧全局路线。

TJU `reg_and_gen/Datasets` 已与 PromLoop 固定提交逐文件核验一致；项目采用 PromLoop 已发布的 81 bp、相对 TSS `[-60,+20]` 标准化序列。对 TJU 不重新执行 TSS 截取或按 strand 反向互补，不因 CSV 未逐条提供 genomic coordinate/strand 而阻塞主 motif 分析。不得声称已经逐条独立验证上游每个负链样本的方向处理；报告中应表述为“采用 PromLoop 发布的 TSS-aligned 标准化序列”。TJU 缺逐条 sigma 标签，因此不作为 RQ2 的主 sigma 证据。

核心实验顺序为：先完成数据审计与清洗 → TJU 六物种主 motif discovery（STREME 双背景分开运行，MEME 代表性对照，DREME 暂不跑）→ FIMO/position/spacing/arrangement → 已知元件对比 → RQ3 分层验证；RegulonDB 用于 E. coli sigma-specific RQ2 和 TJU E. coli 的同物种跨数据集验证；DBTBS 用于 B. subtilis σ/TSS/位点注释；因外部严格留出仅 1 条，不做正式独立统计验证。跨物种、同来源留出、同物种跨数据集三种验证必须分开报告。

embedding / Nucleotide Transformer / clustering / occlusion / SAE 等路线全部后置为探索项。不能把 pooled embedding、attention 或分类性能直接称为 motif/因果发现。

在当前数据方案验收完成前，不运行新的正式 MEME/FIMO。可以构建、核验和冻结处理后的 FASTA/metadata/background/manifest，但不得把历史错误输入上的结果继续当作正式证据。

## 已确认的数据纠错要求

RegulonDB：官方 GFF3 `Sequence` 已按转录方向给出，genomic strand 只作为 provenance，**不得再次按负链 reverse-complement**。从原始 GFF3 重新适配/构建，保留可靠 TSS/坐标证据；序列表与 promoter↔sigma 多对多关联分开保存，不能因序列去重丢失多 sigma 标签。旧二次翻转数据以及直接由其产生的 motif/FIMO/spacing 结果不得作为正式结果；确定仅由错误输入产生且无继续用途的派生产物应删除，原始下载快照不得删除。

DBTBS：原始 release 4.1 页面中的 Gene `Direction`、promoter `Location`、`Absolute position` 和 cis-element sequence 必须重新解析。能可靠恢复 TSS/strand/统一窗口的记录进入 TSS-aligned 核心子集；`Location=ND` 或坐标证据不足的记录保留原始 provenance，但不进入依赖精确 TSS 的位置/间距主分析。不得通过 motif 命中或模型预测反推缺失 TSS。对应参考基因组版本必须记录并校验。

TJU：`positive_samples.csv` 与 `Dataset.csv label=1` 的关系必须逐 ID/序列验收后冻结主正样本入口；原 train/dev/test 只视为同一数据集发布拆分，不自动视作独立验证。已发现的拆分标签冲突、exact/RC 重复及近同源泄漏必须在派生数据中显式隔离，原始 CSV 不做静默修改。TJU `label=0` 固定为 STREME 主背景；已生成的 `dinucleotide-shuffled`（`discovery.dinucleotide_null.fasta`）固定为第二套稳健性背景。明确不新增 GC-matched；两套对照必须分别运行、记录输入哈希与参数，不能根据结果切换主结论。

所有正式派生数据必须保存：来源版本、原始记录 ID、处理规则、排除原因、序列长度、TSS 相对定义（如适用）、strand/orientation 来源（如适用）、sigma 规范化结果、去重/同源簇信息、split、输入输出 SHA256 和处理日志。原始资料与用户已有修改优先保留；删除仅限已确认错误或明确淘汰的派生产物。

## 通用实验约束

按同源簇隔离发现与验证，motif 发现和阈值/特征选择仅用 discovery/train。验证集不重发现、不调阈值。严格 site q 与 raw p 候选分开；阴性结果须保留。课程原编号与新规范编号分开标注，课程必做消融不得由 ML 消融替代。

同步远端前检查工作区；只提交明确授权路径，不自动覆盖原始数据。每次新运行单独保存配置、结果、失败/限制和 manifest。

## PowerPoint default workflow

For any task that creates, edits, redoes, polishes, lays out, reproduces, or visually reviews PowerPoint, PPT, PPTX, slides, or decks, default to loading and following `$ppt-quality-loop`.

Unless the user explicitly requests file-only output or explicitly declines visual review, when PowerPoint MCP is available do not declare PowerPoint work complete without real PowerPoint screenshot visual verification.
