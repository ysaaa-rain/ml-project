# PR01-02：启动子调控模式发现

2026-10-28 **机器学习课程五方法对齐展示已完成**：[结果报告](experiments/aligned_motif_display_20261028/summary/REPORT.md)｜[对齐验收表](experiments/aligned_motif_display_20261028/summary/alignment.tsv)｜[比较图](experiments/aligned_motif_display_20261028/summary/comparison.png)。枯草芽孢杆菌同一冻结Discovery/Development，MEME、STREME、One-hot、BERT、CNN统一阈值规则、TOP5/ALL和天然/Shuffle评价；仅补跑MEME两背景共20个PWM，已有模型不重训。24组候选集合、384次独立标量扫描和24个图中候选核验通过，339个历史文件不变；三张展示图已目视检查。其他物种覆盖单列，不声称CNN全物种完成。[M1/M2统一课程统计](results/eda/m1m2_course_20261028/findings.md)也已补充。以下日期段保留当时的历史进展。

2026-10-27 [枯草芽孢杆菌同数据候选比较](experiments/pr01_02_embedding_v2/real_bsub_comparison_20261027/summary/REPORT.md)：完整10个STREME候选池及两种子CNN、One-hot、BERT使用同一Discovery排序及90%分位阈值、Development正类和两类对照，另展示CNN训练实例位置。结果为探索，MEME完整PWM尚缺；CNN使用序列级监督。

2026-10-27 [CNN真实数据小试](experiments/pr01_02_embedding_v2/cnn_real_pilot_20261027/summary/REPORT.md)：使用传统方法相同的枯草芽孢杆菌数据和划分，训练2个CNN并生成候选，统一扫描正类、天然负类和配对Shuffle。前3候选正类命中0.848/0.818，天然负类0.325/0.247，Shuffle 0.667/0.652；为监督训练和已观察Development的探索，不是完整方法排名。

2026-10-27：以 MEME/STREME 使用的真实数据与划分为共同基础，所有线路统一采用适合课程实践的较宽松候选及局部匹配评价；CNN 尚待收尾，未来横向比较时同样对齐。[重做计划](docs/PR01-02_真实数据统一比较重做计划_20261027.md)｜[候选规范](docs/PR01-02_课程实践候选与展示规范_20261027.md)。

2026-10-27 CNN边界处理多算法种子独立复核：三条件×5组新合成数据×6算法种子（2旧+4新），90次CNN训练、90次OH48聚类、270组PWM池及900组预算评价。30%退化、匹配5候选时，RAW/SHIFT F1约0.506/0.539；5组数据2胜1负2平，4个新算法种子仅1个平均差为正，未通过进一步合成研究门槛。同批OH48约0.510，CNN采用序列级监督，不能据此声称总体优势。算法种子同时控制初始权重与批次顺序；真实Holdout未使用。[完整报告](experiments/pr01_02_embedding_v2/cnn_seed_validation_20261027/REPORT.md)。

2026-10-27 CNN参考无关窗口边界与方向消融完成：60次CNN训练、30次OH48聚类、150组PWM池和420组评价。30%退化、匹配5候选时，同一模型的窗口平移前F1约0.350、平移后约0.564，5组新数据4胜1平、空白误报约0.186/0.170；通过预设进一步合成复核门槛。同批OH48约0.579，CNN仍未稳定超过基线。4次算法运行改善均来自同一模型初始化，下一步优先扩大初始化种子并用新数据复核。25项相关测试、实例/PWM及历史哈希核验通过，真实Holdout未用。[完整报告](experiments/pr01_02_embedding_v2/cnn_boundary_20261027/REPORT.md)。

2026-10-26 已按用户决定收束当前BERT聚类/筛选线路，整理[四轮材料与清理记录](experiments/pr01_02_embedding_v2/closeout_20261026/CLOSEOUT.md)，腾退84份可再生成旧编码缓存（约9.94GiB），保留原始数据、权重、实例/PWM与负面结果。优先老师文档的CNN卷积特征路线，已完成36次CNN训练、18随机模型和18次OH48聚类；30%退化匹配预算F1约0.277/0.547，未通过小试门槛，混合链约0.894/0.871。评分后诊断提示窗口边界错位，下一步优先研究CNN的参考无关边界处理，不恢复BERT常规调参。22项相关测试及独立产物核验通过，真实Holdout未用。[CNN首轮报告](experiments/pr01_02_embedding_v2/cnn_features_20261026/REPORT.md)｜[老师方法核对](experiments/pr01_02_embedding_v2/closeout_20261026/teacher_methods.md)。

2026-10-26 Embedding第四轮独立复核完成：5个新数据种子×3聚类种子，135次聚类、270组PWM池、900组预算评价。30%退化、背景拒绝、双方各5个候选时，BERT F1约0.295，完整One-hot48约0.543；BERT在5/5组数据均较低，候选集误报约0.163/0.146。上一轮均值优势未复现，建议收束当前流程的常规调参，保留为补充探索。37项相关测试及实例/PWM/排名核验通过，真实Holdout未使用。[完整报告](experiments/pr01_02_embedding_v2/budget_validation_20261026/REPORT.md)。

2026-10-26 Embedding第三轮小试完成：48次聚类、240组方向/偏移/背景拒绝对照。30%退化中BERT BG_ONLY的F1从本轮RAW的0提高到约0.482；同处理One-hot约0.418，但BERT仅1/3数据种子更好，候选集空白误报更高，不能宣称稳定优于One-hot。对齐组合未形成一致收益。28项相关测试通过，旧证据和真实Holdout不改。[完整报告](experiments/pr01_02_embedding_v2/refinement_20261026/REPORT.md)。

2026-10-25 Embedding第二轮完成：5数据种子×3聚类种子复核及PCA/筛选消融，含压力测试，合计324次聚类拟合、414组PWM输出，另25个监督定位探针。原流程15%退化F1约0.695；128维约0.790，仍低于One-hot约0.900。30%退化监督定位率约0.730表明信息仍可读出，但不等于无监督发现成功。20项相关测试通过；当前优先研究聚类/实例对齐，未进入真实数据或使用真实Holdout。[完整报告](experiments/pr01_02_embedding_v2/validation_20261025/REPORT.md)。

2026-10-25 Embedding 独立分支推进：完成定位瓶颈与 mini-c 对照小试，6个新合成条件、42个拟合。15%退化中局部重编码+token拼接F1约0.900，仍低于One-hot的0.919；30%退化仍未恢复；mini-c均值方案未显示优势。13项相关测试通过，冻结输入和历史证据不覆盖，未使用真实Holdout。[阶段报告](experiments/pr01_02_embedding_v2/bottleneck_20261025/REPORT.md)｜[协议与复现](experiments/pr01_02_embedding_v2/bottleneck_20261025/README.md)。

2026-10-09 实测完成：新增 42 组 MEME/STREME（27 MEME、15 STREME）与原 12 组 STREME 合计 54 套；540 个报告 PWM 均接 FIMO，344 次扫描成功且无截断。157 项测试通过。8 次辅助 SEA 因 development 仅 1 条正样本失败，保留日志；原 17 主候选和旧结果哈希不变。扩展 holdout 出现率 3 项 BH 显著、σ 特异性 0 项，全部按探索性报告。 [完整结果](reports/M3_扩展双方法与sigma实验_20261009.md)。

2026-10-09 新执行要求：[扩展 M3 规范](docs/PR01-02_expanded_M3_protocol_20261009.md)；MEME/STREME 均覆盖 TJU 六物种双背景及 RegulonDB/DBTBS 整体和可用 σ 组，每套 PWM 接 FIMO。[首次 FIMO 诊断](reports/M3_首次FIMO诊断_20261009.md)未发现明显输入/截断错误，原严格阈值与 17 主候选不变。

后续实测（2026-10-08）：参考库4个短元件PWM本地构建、17主候选不变，已完成18组FIMO及四组BH统计。确认富集/共现均无BH显著，位置/间距因样本不足只描述。详见[参考/FIMO报告](reports/M3_参考库_FIMO与语法检验_20261008.md)。

本地实测更新：2026-10-08 已完成六物种双背景12组STREME，冻结natural主候选17个；全套142项测试通过。正式报告：reports/TJU_STREME_主发现与稳健性_20261008.md。命令及哈希：results/motif/tju_streme_v2_20261008/public_evidence/run_manifest.json。诊断初次输出仅保留本地，诊断参数/哈希已归档。此段记录初次 STREME 完成时的状态；参考库/FIMO 后续实测见上方报告，候选不会按匹配或holdout表现筛除。

**最新文档入口：[docs/README.md](docs/README.md)**；M1/M2 阶段报告：[M1M2_阶段调查与描述性发现](docs/M1M2_阶段调查与描述性发现_20261008.md)。历史组员原稿仅从 Git 固定提交追溯。

更新：2026-10-08。当前数据版本：`data/processed/pr01_02_data_v2/`。

TJU 六物种为主数据；RegulonDB 负责 E. coli sigma 注释/专项与独立来源候选，DBTBS 负责 B. subtilis sigma、TSS、文献注释（因严格外部 holdout 仅 1 条，不做正式独立统计验证）。

本研究采用 PromLoop 发布的标准化 promoter 序列，并按其数据集定义视为已经完成 TSS 对齐及必要的方向标准化。81 bp 对应相对 TSS 的 [-60,+20]，TSS 为局部第61位；TJU 不重新截取或反向互补。genomic coordinate、逐条 strand 和原始数据库映射不是 TJU 主实验必要输入。

2026-10-05 的数据重建提交未运行正式 MEME/STREME/DREME/FIMO；后续 STREME 已由本地 agent 接手，实时进展应以新的结果 manifest 和本地执行日志为准。旧 RegulonDB 3,807 条处理输入及其实验结果已从项目移除，不能继续引用旧覆盖率、σ显著性或位置结论。

## 数据链条与实验主线

TJU 六物种 positive promoter → 清洗/去重/泄漏隔离 → STREME + MEME（natural 与 dinucleotide 两背景）→ FIMO → motif position/spacing/cross-species。

RegulonDB E. coli → σ关系表/σ-specific motif → 与 TJU E. coli 做去重叠 cross-dataset validation。

DBTBS B. subtilis → 经 NC_000964.2 核验的 σ/TSS/位点与文献注释；不作为独立外部统计验证集。

同来源holdout与不同来源验证分别报告；共享序列用于注释，不重复计作独立验证。不同发布来源不保证原始研究独立。

## 当前数据数量

| 数据 | discovery 正类 | development 正类 | holdout 正类 |
| --- | ---: | ---: | ---: |
| dbtbs_bsub/core | 296 | 39 | 97 |
| regulondb_ecoli/core | 1284 | 160 | 349 |
| regulondb_ecoli/extended | 1421 | 191 | 408 |
| tjupan_bacillus_subtilis/main | 473 | 67 | 151 |
| tjupan_baumannii/main | 1099 | 162 | 279 |
| tjupan_bradyrhizobium/main | 1364 | 171 | 448 |
| tjupan_diphtheria/main | 1186 | 156 | 311 |
| tjupan_escherichia_coli/main | 1143 | 163 | 329 |
| tjupan_staphylococcus/main | 1480 | 248 | 478 |

跨来源可用样本以v2生成清单为准。RegulonDB核心中完全不关联任何TJU正负输入且属于冻结留出的外部候选250条；DBTBS核心相同条件只有1条，无法支撑稳健的独立来源结论。不能退回重复样本或降低隔离规则来增加显著性；DBTBS当前以σ注释和可靠坐标数据为主要补充，若要充分验证需要新独立采集来源。

Promoter的Strong/Confirmed不等于σ关联已逐条实验确认，σ证据状态单独保留。小组需总数≥50、发现≥30、留出≥20才作为确认候选；20–49探索，<20描述。详细实际分组见sigma_group_summary.tsv。

E.coli 50nt强度原文件保留为Exp5候选，单位/条件/聚合定义未确认；kucao_test保留待核，不进入正式输入。确定不用的酵母/蓝细菌/派生物化特征/图片资产已移出项目。历史TJU实验可保留为探索，不能替代v2冻结检验。

## 构建与复现

```bash
python -m pip install -r requirements.txt
python -m preprocessing.rebuild_audited_data
python -m preprocessing.finalize_audited_inputs
python -m pytest -q
```

重建命令拒绝覆盖已有v2目录；要复现到新目录，可在Python调用run(output=Path("新的项目内目录"))并将同一目录传给finalize的run(base=...)。固定原始快照必须存在，按manifest原始SHA256核对；不会自动重新下载或覆盖现有数据。

首份修复后的辅助FASTA：`data/processed/pr01_02_data_v2/regulondb_ecoli/core/discovery.positive.fasta`（1284条）。正式主实验从各`tjupan_*/main/discovery.positive.fasta`开始；每个物种使用两套已生成、按manifest验收的 STREME 对照：`discovery.natural_control.fasta` 为正式主背景，`discovery.dinucleotide_null.fasta` 为组成控制的稳健性背景。两次 discovery STREME 分开运行；不新增 GC-matched。详见[2026-10-08 主分析与语法统计预注册](docs/PR01-02_motif_reference_and_grammar_protocol_20261008.md)。

当前不得使用旧实验脚本默认的20260930/20261003输入路径直接启动实验。2026-10-09 扩展实验入口已锁定 v2 与参数；新增运行、完成状态与限制以新 manifest/报告为准。

详细方案：[数据源选择与处理规范](docs/PR01-02_数据源选择与处理规范_20261005.md)。修复/清理记录：[数据修复与输入验收](reports/数据修复与输入验收_20261005.md)。
