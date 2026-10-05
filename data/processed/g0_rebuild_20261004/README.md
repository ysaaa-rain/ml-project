# 已核查的G0协作输入

发布日期：2026-10-05。用户明确要求共享最新方案下已重建的数据。本目录随仓库拉取，队友无需重新下载或清洗即可继续M2。完整方法见[数据重建报告](../../../docs/09_G0重建与冻结划分_20261004.md)。

## 从哪里开始

- `standard_metadata.tsv`：三个窗口的全部有效位点与审计注释。筛选`upstream == 80`得到101bp主窗口；做Gold代表EDA再筛选`gold_representative == True`。`gold_sigma_labels`为JSON列表，多sigma不能取首项当单标签。
- `window_80_20/train/all_gold.fasta`：751条主窗口训练代表；`holdout/all_gold.fasta`：349条固定留出代表。主窗口TSS序列内索引81，所有序列已转录定向，勿再次根据strand翻转。
- `window_60_20/`、`window_100_20/`：81/121bp敏感性窗口，同代表、同划分。不能把三个窗口合并当三倍样本。
- `cluster_members.tsv`、`excluded_source_records.tsv`、`sequence_hash_and_selection.tsv`：簇、排除和选择依据。不得重新随机划分。
- `single_sigma_ml.tsv`：961条严格ML候选，不表示已训练分类器。
- 各训练目录下`*_N1.fasta`为二核苷酸保持背景；`*_other_matched_primary.fasta`与`*_other_matched.fasta`为逐条匹配输入/对照。Sigma70背景B只有159对；Sigma28训练仅8条，不作正式发现。
- `dbtbs_audited_candidates.tsv`：698条DBTBS候选，证据待审阅且TSS/strand未确认；不是已经通过的Gold数据。

总量：3,822有效位点、3,296联合同源簇、1,100 Gold代表；发现集与留出集按簇冻结。全部位点EDA与Gold代表EDA分别报告，按sigma展开后的组数量之和不等于独特代表数。发现参数和背景设计以train为依据，holdout不参与调参。

## 证据与文件完整性

`share_manifest.json`提供已发布文件的字节数和SHA256、原构建参数、输入来源hash、独立核查及补充产物hash。原始本机`run_manifest.json`含机器环境信息，保留本地；其原字节hash在共享manifest登记，共享文件不冒充该原文件。科学配置、结果状态和构建限制原样保留。

发布时逐项校验冻结构建与补充manifest，并核对Git暂存blob的字节hash；FASTA/TSV采用LF，避免跨系统换行改变DNA文件hash。分享仅增加传递文件，不重跑或改变数据重建结果。

当前完成的是来源/数据准备及内部同源隔离核查，正式新motif发现、FIMO扫描、独立验证和ML实验尚未运行。raw来源下载与其他历史大数据没有随本目录上传；需要从原始来源完整重建时，须按冻结配置另外取得raw输入。
