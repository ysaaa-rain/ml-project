# M3 扩展实验执行规范（2026-10-09）

本规范落实用户新要求，覆盖旧规范中“MEME 仅挑代表物种”和“DBTBS 仅注释”的方法范围。TJU 六物种仍是主序列来源；原 natural STREME 的 17 个主候选与原严格 FIMO/holdout 结果保留，不用新增结果替换。新增分析均标为扩展探索，不能把已经查看过的 holdout 再称作全新确认性检验。

## 实验矩阵

- TJU：六物种分别运行 MEME 与 STREME；同一 discovery positive 分别配 natural_control 和 dinucleotide_null。已有 12 组 STREME 复用完整输出，补跑 12 组 MEME。
- RegulonDB、DBTBS：core 整体及每个 discovery≥10 的单 σ 组，均运行 MEME 与 STREME，使用各自已冻结的 dinucleotide_null。RegulonDB 六组；DBTBS 七组；DBTBS 另外六组 n<10 登记跳过两种方法。10–19 条组保留低功效及工具失败记录。
- 共 42 组新增发现、12 组已有发现；每套成功输出的全部 PWM 均接独立 FIMO，包括不显著候选。没有 PWM 的运行记录原因，不虚构扫描。
- DBTBS 自身整体/σ 探索不等于独立外部验证；与 TJU 的高重叠限制仍有效，外部 holdout 1 条不做统计实验。

## 固定参数及方法差异

工具 MEME Suite 5.5.9；seed=20261005；宽度 5–15；最多 10 个 motif。仅 discovery 发现。

MEME：`-dna -revcomp -mod zoops -objfun de -test mhg -neg <control> -hsfrac 0.5 -searchsize 100000 -nmotifs 10 -minw 5 -maxw 15 -seed 20261005 -maxsize 10000000`；显式 discovery control 零阶、互补对称背景与 `-markov_order 0`。STREME：`--dna --p <positive> --n <control> --order 2 --minw 5 --maxw 15 --thresh 0.05 --nmotifs 10 --seed 20261005`。5.5.9 的参数解析中 `--thresh` 重置 motif 数量，因此实际命令先 thresh 后 nmotifs。

MEME 内部评估比例 0.5，STREME 0.1；内部拆分由工具实现，不保证同源簇/正样本与其 shuffle 成对隔离，原生显著性只描述，簇级下游检验另算。背景模型、优化目标及显著性不同，不能将两者的 p/E 直接等同。MEME 按原生 E≤0.05 描述，STREME 按原生内部测试 p≤0.05，并要求内部测试正负样本各≥5。全部报告 PWM 仍扫描，不按下游好坏筛除。

## FIMO 与诊断

每套 PWM 扫描对应 discovery/development/holdout 正样本及背景；固定 discovery control 零阶背景，双链，pseudo=0.1，max-stored-scores=5000000，正式 `--qv-thresh --thresh 0.05`。每个 motif 的原生 site q 家族随 case/split 扫描确定；不同扫描家族的 q 不直接横比。空分区登记跳过。

原 17 主候选诊断仅 discovery/development：全 p<1 位点表、分布、每序列最佳分数、SEA 序列层面富集、背景与截断日志。新增诊断 raw p<1e-4 也仅 discovery/development。SEA 原生序列富集为辅助诊断，不代替原严格 site q 或独立簇统计。shuffle 与原序列成对，不能把两者当独立生物样本。

σ 分组候选另扫描同来源各单 σ 组真实正样本池，以目标 σ 与其他 σ 作对照；不能凭候选出自该组或相对 shuffle 富集宣称 σ 特异。

## 扩展汇总统计

只使用严格位点。每序列/每 motif 选 q 最小、score 最大、start 最小、stop 最小、strand 字典序第一的主位点，TSS 相对中心为 `(start+stop)/2-61`。位置只报告描述性分布；不据此改候选。

出现率以同源簇为单位，规则不依赖命中结果：天然对照/σ 对照的混标签簇排除，纯标签簇取字典序最小 ID；shuffle 对照按正样本簇选一个 ID 与其 shuffle 成对，采用 discordant pairs 的精确二项检验。天然正负与目标 σ/其他 σ 用双侧 Fisher。每项报告有效簇数、出现率差与原始 p。汇总阶段统一补充所有效应的 95% 保守区间：两个边际各取 97.5% Clopper–Pearson，再用 Bonferroni 差值界；成对 shuffle 使用正/负 discordant 事件边际，不按独立两组估计。区间规则不改变任何检验或候选。

BH 分别按 split 和比较目的跨所有扩展 case/PWM 校正：`natural_presence`、`paired_null_presence`、`sigma_specificity`。空组无法检验记 NA，不补显著性。小组的推断限制始终保留。新增 holdout 检验全部探索性，不优化参数。原 17 的既有四组检验家族不被替换。

MEME/STREME 对应 source/group/control 的 PWM 由 Tomtom Pearson、min-overlap=5 比较；保留全部输出匹配，不能按结果重选主候选。参考匹配仅 E. coli overall/Sigma70 与 B. subtilis overall/SigA 使用现有适用参考；其他 σ/物种无合格参考时记未注释，不套用 σ70/SigA。

## 保存与发布

本地保留完整原生命令、版本、输入输出 SHA256、失败日志、矩阵与全位点表。公开提交代码、规范、聚合表和运行证据；受限原始序列、完整含序列的 MEME/XML/FIMO、RegulonDB/DBTBS 派生 PWM/参考库只留本地 tmp。原始数据不删。
