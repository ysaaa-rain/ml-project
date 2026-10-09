# PR01-02：启动子调控模式发现

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
