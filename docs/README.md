# PR01-02 文档导航（2026-10-09）

2026-10-28 **课程展示入口**：[同数据五方法结果](../experiments/aligned_motif_display_20261028/summary/REPORT.md)｜[逐方法对齐验收](../experiments/aligned_motif_display_20261028/summary/alignment.tsv)｜[六物种覆盖清单](../experiments/aligned_motif_display_20261028/summary/coverage.tsv)｜[M1/M2统一统计](../results/eda/m1m2_course_20261028/findings.md)。枯草芽孢杆菌完整比较已补齐MEME；按机器学习课程目标收束，其他物种覆盖如实标注，不追加BERT/CNN调参。

2026-10-27 [真实数据同口径比较与CNN位置图](../experiments/pr01_02_embedding_v2/real_bsub_comparison_20261027/summary/REPORT.md)：同一枯草芽孢杆菌数据，完整STREME候选池与CNN、One-hot、BERT按统一规则扫描；MEME PWM待补齐。

2026-10-27 [CNN真实数据小试](../experiments/pr01_02_embedding_v2/cnn_real_pilot_20261027/summary/REPORT.md)：沿用B. subtilis真实数据、两次CNN训练与统一候选扫描；保留天然及Shuffle对照，尚未完成全部方法的横向比较。

2026-10-27 [真实数据统一比较重做计划](PR01-02_真实数据统一比较重做计划_20261027.md)｜[候选展示规范](PR01-02_课程实践候选与展示规范_20261027.md)：采用 MEME/STREME 的基础数据，所有方法统一使用适合课程实践的候选和局部匹配评价；CNN 待收尾后再对齐。

2026-10-27 [CNN边界处理多算法种子独立复核](../experiments/pr01_02_embedding_v2/cnn_seed_validation_20261027/REPORT.md)：5组新数据×6算法种子×三条件，90次CNN训练、90次OH48聚类和900组评价。30%退化匹配候选下RAW/SHIFT F1约0.506/0.539，5组数据2胜1负2平；4个新算法种子仅1个平均差为正，未通过进一步研究门槛。同批OH48约0.510；CNN使用序列级监督，差距受单组数据影响。算法种子同时控制初始化与训练顺序，真实Holdout未用。

2026-10-27 [CNN参考无关窗口边界与方向消融](../experiments/pr01_02_embedding_v2/cnn_boundary_20261027/REPORT.md)已完成：60次CNN训练、30次OH48聚类、150组PWM池、420组评价。30%退化匹配候选下，CNN平移前/后F1约0.350/0.564，5数据种子4胜1平，空白误报约0.186/0.170；同批OH48约0.579。10次算法运行的4次收益均来自同一初始化，需新数据与更多初始化复核；不能声称稳定超越OH48。25项相关测试及产物核验通过，真实Holdout未用。

2026-10-26 当前Embedding状态：[BERT收束与清理材料](../experiments/pr01_02_embedding_v2/closeout_20261026/CLOSEOUT.md)｜[老师推荐方法核对](../experiments/pr01_02_embedding_v2/closeout_20261026/teacher_methods.md)｜[CNN局部特征首轮结果](../experiments/pr01_02_embedding_v2/cnn_features_20261026/REPORT.md)。已腾退84份旧编码缓存约9.94GiB，保留复现实验依据。CNN新路线完成36训练+18随机模型+18次聚类、72组PWM池与180组评价；弱模式F1约0.277/One-hot48约0.547，混合链约0.894/0.871，尚无稳定总体优势。探索性诊断提示卷积窗口边界问题，后续应另开协议检验；22项相关测试及独立产物核验通过，真实Holdout未用。

2026-10-26 Embedding第四轮：[背景拒绝独立复核与候选预算匹配](../experiments/pr01_02_embedding_v2/budget_validation_20261026/REPORT.md)。135次聚类、270组PWM池、900组评价已完成；30%退化双方各5个候选时，BERT F1约0.295，完整One-hot48约0.543，5/5新数据种子均较低且候选集误报更高。背景拒绝有自身改进，但稳定优势未复现；建议收束当前流程的常规调参，保留补充研究定位。37项相关测试及产物核验通过，真实Holdout未用。

2026-10-26 Embedding第三轮：[候选片段对齐与背景拒绝报告](../experiments/pr01_02_embedding_v2/refinement_20261026/REPORT.md)。48次聚类、240组输出已完成；背景拒绝改善BERT自身的30%退化恢复（0→约0.482），但只在1/3数据种子超过同处理One-hot且误报更高，仍无稳定优势。源码、冻结协议、图表、候选集误报及实例/PWM核验见同目录。

2026-10-25 Embedding第二轮：[多种子复核与PCA/实例瓶颈报告](../experiments/pr01_02_embedding_v2/validation_20261025/REPORT.md)。已完成324次实际聚类、414组PWM及25个监督探针；原流程15%退化约0.695，128维约0.790，One-hot约0.900。探针说明局部信息仍可读出，未证明无监督发现优势；真实数据门槛未满足。源码、两份冻结协议、逐数据种子统计、图表和核验见同目录。

2026-10-25 Embedding探索更新：[定位瓶颈与mini-c阶段报告](../experiments/pr01_02_embedding_v2/bottleneck_20261025/REPORT.md)。已完成6条件、42个合成拟合；局部重编码+token拼接在15%退化中约0.900，One-hot约0.919；30%退化与mini-c均值方案仍失败。两数据种子仅作描述性诊断，未进入真实数据或使用真实Holdout。现行主线规范和冻结结果不改写。

2026-10-09 实测完成：新增 42 组 MEME/STREME（27 MEME、15 STREME）与原 12 组 STREME 合计 54 套；540 个报告 PWM 均接 FIMO，344 次扫描成功且无截断。157 项测试通过。8 次辅助 SEA 因 development 仅 1 条正样本失败，保留日志；原 17 主候选和旧结果哈希不变。扩展 holdout 出现率 3 项 BH 显著、σ 特异性 0 项，全部按探索性报告。 [完整结果](../reports/M3_扩展双方法与sigma实验_20261009.md)。

最新执行入口：[扩展 M3 规范](PR01-02_expanded_M3_protocol_20261009.md)（2026-10-09）。MEME/STREME 双方法与 σ 整体/分组、每套 FIMO，覆盖旧文档的方法范围；DBTBS 不恢复独立外部验证。

本目录的**现行规范仅以下列文档为准**，不要从文件日期或旧组员报告推断当前实验状态。数据冻结目录为 `data/processed/pr01_02_data_v2/`。本次文档整编不更改任何数据、脚本参数、冻结 FASTA 或已有实验结果；本地 agent 正在启动的 STREME 不受本次整编影响。

| 阅读顺序 | 文件 | 定位 |
| --- | --- | --- |
| 1 | [00_项目约束与任务定义](00_项目约束与任务定义.md) | 老师的 M1/M2 阶段任务、RQ 与实验边界 |
| 2 | [PR01-02_数据源选择与处理规范_20261005](PR01-02_数据源选择与处理规范_20261005.md) | v2 来源、数量、证据、排除与泄漏规则 |
| 3 | [07_后续实验执行规范_20261004](07_后续实验执行规范_20261004.md) | 按 10-08 最终决定更新的主流水线 |
| 4 | [PR01-02_motif_reference_and_grammar_protocol_20261008](PR01-02_motif_reference_and_grammar_protocol_20261008.md) | 双 STREME 背景、参考库构建、FIMO 与 BH/grammar 预注册 |
| 5 | [M1M2_阶段调查与描述性发现_20261008](M1M2_阶段调查与描述性发现_20261008.md) | 已有证据可支持的 M1/M2 结果，不冒充 M3 实验 |
| 6 | [M1M2_待补只读分析与验收_20261008](M1M2_待补只读分析与验收_20261008.md) | 现有报告尚未覆盖的三来源分析及本地补算方案 |
| 7 | [运行说明](运行说明.md)、[当前进度与协作交接](当前进度与协作交接.md) | 入口与交接；实况以最新 commit、manifest、任务日志为准 |

## 证据和历史

- [M1/M2 精选表格证据](evidence/M1M2_20261007/README.md)：保留从乙组统计、丙组归档中筛选的可追溯 TSV 原字节，不代表全部发现已被独立复算。
- [数据审计机器证据](data_audit_evidence_20261005/)：manifest、来源哈希、源数据异常与清理记录，不删改。
- [PromLoop 来源哈希核验](PromLoop数据来源哈希核验_20261003.md)：六物种原始 CSV 的上游来源与 SHA256。
- [组员原稿与旧方案溯源](archive/README.md)：甲、乙、丙原稿通过 Git 固定提交追溯；旧 G0 审计/运行日志与旧实验记录在 archive 下，只供历史查询，不作新方案依据。
- [课程参考资料](参考资料/) 和 [会议纪要](会议纪要/) 为原始资料，仍保留。

## 阶段与问题不要混为一谈

**M1（问题定义/数据构建）和 M2（探索性数据分析）是课程阶段，不是 RQ1/RQ2，也不是运行 MEME/STREME/FIMO 的阶段。** M1/M2 只读检查、计数、分析已知 −10/−35 候选位置、GC/长度和 σ 标签分布，并报告观察和限制。**M3 起**才是 de novo motif discovery，六物种先 STREME，两套背景分开跑。

本轮已废止 GC-matched 分组/背景要求。采用真实 `label=0` 作为 STREME 主对照、已冻结 `dinucleotide_null` 作为稳健性对照。DBTBS 不做独立外部统计验证；仍供 B. subtilis σ/TSS/文献注释。TJU 正式流水线仅使用 Dataset.csv；positive_samples.csv 原件保留 raw 归档。

> **目录整编状态（2026-10-08）：已实施。** 甲乙丙 84 个历史交付文件退出当前工作树，15 份重要 TSV 按原 blob SHA 转存至 `evidence/`，6 份独立历史审计记录转入 `archive/history/`；过时的根级重复说明已清除。组员完整原稿可从 `archive/README.md` 的固定 Git 提交还原。本次未执行 M1/M2 尚欠缺的三来源只读复算，亦未参与本地 agent 正在进行的 STREME。

## M3 中期汇报（2026-10-10）

- [M3 中期 RQ1 综合分析与汇报初稿](../reports/M3_中期RQ1综合分析与汇报初稿_20261010.md)：MEME/STREME 双方法、已知元件、FIMO/SEA、探索性 holdout 和现有局限。
- [本地 agent 的 RQ1 已知位点重合补证任务](M3_RQ1_本地agent只读补证任务_20261010.md)：指定本地原生位点/参考注释的只读汇总口径，禁止虚构实际位点重合率。

上述两份是阶段性分析与工作交接，不覆盖冻结的原始 17 个 STREME 主候选或既有实验清单。
