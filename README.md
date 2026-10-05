# PR01-02：启动子调控模式发现

更新：2026-10-05。当前数据版本：`data/processed/pr01_02_data_v2/`。

TJU 六物种为主数据；RegulonDB / DBTBS 主要负责 sigma annotation、sigma-specific motif 与隔离后的不同来源验证。

本研究采用 PromLoop 发布的标准化 promoter 序列，并按其数据集定义视为已经完成 TSS 对齐及必要的方向标准化。81 bp 对应相对 TSS 的 [-60,+20]，TSS 为局部第61位；TJU 不重新截取或反向互补。genomic coordinate、逐条 strand 和原始数据库映射不是 TJU 主实验必要输入。

本轮仅修复和构建数据，没有运行正式 MEME/STREME/DREME/FIMO。旧 RegulonDB 3,807 条处理输入及其实验结果已从项目移除，不能继续引用旧覆盖率、σ显著性或位置结论。

## 数据链条与实验主线

TJU 六物种 positive promoter → 清洗/去重/泄漏隔离 → MEME/STREME → FIMO → motif position/spacing/cross-species。

RegulonDB E. coli → σ关系表/σ-specific motif → 与 TJU E. coli 做去重叠 cross-dataset validation。

DBTBS B. subtilis → 可靠坐标恢复/σ-specific motif → 与 TJU B. subtilis 做去重叠 cross-dataset validation。

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

首份修复后的辅助FASTA：`data/processed/pr01_02_data_v2/regulondb_ecoli/core/discovery.positive.fasta`（1284条），配对`discovery.dinucleotide_null.fasta`。正式主实验从各`tjupan_*/main/discovery.positive.fasta`开始；STREME背景选择需在正式实验配置中固定，不能边看结果边换。

当前不得使用旧实验脚本默认的20260930/20261003输入路径直接启动实验。正式实验启动前需要将MEME/STREME/FIMO入口更新到v2并锁定参数；本轮没有执行它们。

详细方案：[数据源选择与处理规范](docs/PR01-02_数据源选择与处理规范_20261005.md)。修复/清理记录：[数据修复与输入验收](reports/数据修复与输入验收_20261005.md)。
