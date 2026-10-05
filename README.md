# PR01-02 多物种启动子调控模式发现

研究细菌启动子的motif序列、相对TSS的位置及间距/排列，回答已知元件恢复、sigma特异性和跨数据稳定性（RQ1–RQ3）。

## 当前依据

用户指定的[实验执行规范](docs/07_后续实验执行规范_20261004.md)与[原始指导](docs/参考资料/PR01-02_后续实验指导_20261004.txt)为后续路线依据。两份长期[参考资料](docs/参考资料/README.md)已保存。当前证据与差距见[全方位审计](docs/08_全方位核查审计_20261004.md)、[交接](docs/当前进度与协作交接.md)。

## 后续路线

数据/evidence/方向/同源簇审计 → 新1发现 → 新2定位 → 新3known recovery → 新4sigma grammar → 新5独立验证 → 新6sigma预测消融 → 新7CNN/归因/TF-MoDISco独立motif发现。

| 来源 | 后续角色 | 已确认边界 |
|---|---|---|
| RegulonDB E. coli K-12 | 经审计的主发现输入 | 真基因组窗口与统一同源簇划分已重建；Gold限Strong/Confirmed；multi_sigma保留；旧结果不覆盖 |
| DBTBS B. subtilis | 经标注审计的主发现候选 | 本地698记录缺TSS/strand，暂限非位置发现；不能填猜测坐标 |
| TJU Pan/PromLoop六物种 | 其他来源/物种的验证候选 | CSV无sigma及可确认方向；其旧60PWM发现结果保留原角色，不能以同批样本验证自身发现 |
| strenth E_coli.txt | 课程原选做5历史扩展 | 弱关联，标签条件/单位仍需追溯；不是新5验证 |

课程原实验编号与新编号分别标注，课程必做消融不能由ML A–D替代。embedding是有明确motif实例提取、PWM重构和验证链的扩展。

## 2026-10-04审计状态

- 来源与代码确认旧RegulonDB1,985条负链二次反向互补，旧TSS位置/方向解释暂停引用；已修复方向/窗口代码，未重跑旧实验。
- 原RegulonDB123组、DBTBS21组相同序列具有多个sigma标签。修复keep-first标签丢失，保留完整重复注释；严格ML排除multi_sigma。
- 重建历史80/20六组留出、全部12份primary hash匹配；检出4对跨split近重复，涉及2条留出记录。旧留出不满足新同源簇隔离要求。
- 80项已定位的产物hash均匹配；通过LF属性修正Windows CRLF转换问题。哈希匹配不等于全部统计过程或生物结论通过。
- 旧六物种FIMO353个位点与跨物种比较阴性等均保留原参数和证据级别。新1–5尚未全部按新标准验收，新6/7未做。

## 本地运行

2026-10-05依用户要求，已将[重建数据目录](data/processed/g0_rebuild_20261004/)纳入共享，拉取即可继续M2，无需重新下载/清洗。范围包括标准表、三个窗口FASTA、固定划分与两套训练背景；原机器环境manifest仍本地保留。源文件至共享输入的逐项核查与M1未通过项见[M1验收说明](docs/11_M1验收与数据共享_20261005.md)。RegulonDB内部M1通过，不表示DBTBS/外部来源和整个多来源M1已全部完成。

2026-10-04已完成[RegulonDB G0重建与冻结划分](docs/09_G0重建与冻结划分_20261004.md)：3,822个位点、3,296个联合同源簇、1,100个Gold代表，三个窗口穷举跨split核查均无符合冻结标准的近同源对。两套训练背景已准备；Sigma28仅8条训练样本，Sigma70其他sigma背景仅159匹配对。DBTBS证据/坐标与外部来源隔离仍待审计，新motif实验尚未运行。

当前Windows尚无可运行MEME Suite；远端历史构建记录不能当作本机环境。当前Python依赖版本也不同于requirements锁定值。原始下载数据在被Git忽略的 `data/raw/course_share/`、RegulonDB和DBTBS本地目录中；本轮仅明确共享上述G0派生输入，其他原始下载与大型缓存不提交。

```powershell
python -B -m experiments.preflight --output tmp/preflight_audit.json
python -B -m experiments.audit_repository --output-dir results/data_audit/<new_run>
```

审计输出必须新目录；不运行新motif实验。原参数复现与前置文件见[运行说明](docs/运行说明.md)。旧留出复现必须显式 `--legacy-record-split`，不能用其作为新规范的默认划分。

## 目录与协作

队友M1/M2原稿与当前实现的衔接、统一数据口径及建议分工见[协作核查](docs/10_仓库协作核查与报告衔接_20261004.md)。原DOCX为问题定义/分析计划稿，正式M2结果需引用新输入EDA；旧暂存修改不能整批恢复。

`preprocessing/` 数据处理；`experiments/` 审计与统计入口；`results/` 审阅过的摘要/图/manifest；`docs/` 指导、验收、审计、运行和交接；`reports/` 历史阶段报告（受影响内容已加审计说明）。

保留历史结果与会议原稿；每次运行记录Git状态、软件、命令、参数、输入输出hash、随机种子、split/cluster和限制。每次推送更新交接。完整实验报告、课程必做项与答辩材料仍须按证据交付。
