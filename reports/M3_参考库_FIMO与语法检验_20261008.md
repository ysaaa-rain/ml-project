# M3 后续：真实短元件参考、固定 PWM FIMO 与 holdout 语法统计

日期：2026-10-08。先完成并推送M3双背景STREME（远端6ad477e），再按既定顺序构建参考、比较冻结候选、核验discovery/development、最后评估holdout。natural主候选17个身份及PWM SHA256保持不变；所有未匹配/零命中候选均保留。

## 1. 真实参考库 v1

本地文件：`tmp/reference_library_v1_20261008/reference_library_v1.meme`及同名TSV，manifest含来源/产出SHA256与构建规则。数据库数据再分发许可未确认，因此逐条短序列、完整原始API响应和数据库派生参考矩阵不推送；公开引用、规则、计数和哈希。

| reference | distinct verified genomic sites |
| --- | ---: |
| E. coli Sigma70 minus10 | 26 |
| E. coli Sigma70 minus35 | 20 |
| B. subtilis SigA minus10 | 218 |
| B. subtilis SigA minus35 | 198 |

RegulonDB官方GraphQL下载2605/2605 operon关联对象，2286个唯一关联promoter；不认证所有孤立promoter覆盖。接口getDatabaseInfo最新版本为14.5.0、releaseDate2026-01-28、NC_000913.3。只取S/C、实验PMID、无COMP混合证据的Sigma70显式minus10/minus35六碱基boxes，并核对基因组短序列及方向。不是对完整81bp promoter固定窗口切片。

DBTBS4.1官方[terms](https://dbtbs.hgc.jp/ver4/terms.html)说明红色对应文献中的精确结合片段。只从已可靠恢复坐标的记录取连续6bp红色片段、有PMID和实验方法、无HM。SigA类别以相对TSS范围事前规则区分：minus10起点−15..−7，minus35起点−38..−30。该分类是坐标规则，不是原页面另有显式−10/−35字段；不截取未标红碱基来补齐6bp。唯一基因组区间作为去重复单位，相同短序列位于不同区间不直接删除；不同区间也不保证研究或同源独立。每列总伪计数1（四碱基各0.25），无额外序列对齐调整。

证据边界：这些是实验支持promoter中的数据库短元件注释，不认证每个位点每个碱基都经突变实验证实。参考来源与TJU可能共享原始研究或序列，Tomtom是注释相似比较，不是独立验证。UP元件和其余四物种均未获得本轮合格参考；不引入E.coli通用共识替代。

## 2. 固定17候选的已知元件比较

仅E.coli/B.subtilis用对应物种参考各2个PWM。Tomtom Pearson、minimum overlap5、threshold1，保存全部匹配与offset/overlap/orientation/p/E/q。小库警告实际出现（每query仅4个含方向p值，pi0不足100时设为1），只作描述。

B.subtilis `1-TATAATAAAWA` 对 SigA minus10最佳匹配 overlap6，p=0.00260839、q=0.0104336；不能解释为已验证功能。E.coli各候选最佳q均>0.05。其余四物种共11个候选记录无合格参考。`all_17_candidate_annotations.tsv`记录全部17个，不按匹配程度筛除；未匹配始终为未注释候选。

## 3. 固定PWM FIMO扫描

FIMO5.5.9，只用17个冻结PWM。同物种同split将positive与natural_control合并扫描，保留标签映射本地，分别discovery/development/holdout共18组。背景仅用discovery natural负类估计zero-order、方向对称、总伪计数1，此背景对该物种各split固定；不使用holdout估计背景。显式motif pseudo0.1、双链、`--qv-thresh --thresh 0.05 --max-stored-scores 5000000`。最大潜在搜索位点小于此上限，日志无丢弃/截断警告。先完成全部discovery/development核验（各组非零/零位点照常保存）才开始holdout，未改阈值。

FIMO native q在每个motif的同物种同split正负合并搜索集合中计算。其内部使用BH/原生pi0估计实现；不从截断表另算q、不把site q称为后续检验BH q。没有采用`--text`（它禁用q计算）。完整native命中表含真实短匹配序列，留在tmp不公开。参见[官方FIMO说明](https://meme-suite.org/meme/doc/fimo.html)。

| species subset | holdout strict sites q≤0.05 |
| --- | ---: |
| B. subtilis | 1 |
| Baumannii | 0 |
| Bradyrhizobium | 0 |
| Diphtheria | 1 |
| E. coli | 0 |
| Staphylococcus | 14 |

## 4. 同源簇统计与四组BH

分析配置在任何项目holdout扫描前冻结。每个纯标签同源簇以sequence_id字典序最小记录作为唯一统计代表，选代表不查看命中；2个含正负混合标签的簇排除于确认统计并登记，原始序列/位点不删除。这个保守规则避免一个簇贡献多个独立样本，因而site总数与统计代表命中数不同。

每promoter每motif主site按q最小、score最大、start最小，额外stop/strand只作完全同分破平局；中心=(start+stop)/2−61，gap=下游start−上游stop−1，重叠为负，完整原始位点仍保留。Fisher按独立簇代表出现/共现；OR、log-Wald95%CI，零格使用Haldane0.5仅计算CI，完全无命中OR记NA而非补造效应。位置/gap按5000次簇单位置换KS、至少每组20个独立命中簇，未达门槛记NA；合格时中位差CI用2000次簇bootstrap。

| BH family | planned | tested | BH q≤0.05 | descriptive only |
| --- | ---: | ---: | ---: | ---: |
| single motif presence | 17 | 17 | 0 | 0 |
| pair co-occurrence | 17 | 17 | 0 | 0 |
| TSS-relative position | 17 | 0 | 0 | 17 |
| pair gap | 17 | 0 | 0 | 17 |

Staphylococcus首个候选在统计代表中正类9/361、负类0/229，raw p=0.014404，BH q=0.244868；不报告确认富集。位置/间距样本不足，不能用探索raw-p位点替换严格位点来拼出显著性。零命中不等于motif不存在，也不能从当前阴性结果判断STREME候选没有功能。

## 5. 验收、复现与未完成项

全套147项测试通过（16项既有依赖警告），另新增5项测试覆盖BH顺序、gap重叠/嵌套、主site排序、全零效应和红字短片段解析。输入/候选/参考/工具/输出SHA256保存，17候选保留状态逐项核对。原始数据及v2未修改。

复现代码：`preprocessing.download_reference_boxes`（原始API快照在tmp，拒绝覆盖）；`preprocessing.build_reference_library_v1`；`experiments.run_frozen_m3_followup`（拒绝覆盖）；`experiments.analyze_frozen_m3_grammar`。完整命令见run_manifest和analysis_config。数据库来源记录/序列仅本地，公开汇总在`results/motif/m3_followup_20261008/`。

当前未执行：固定主PWM的dinucleotide背景FIMO扫描（已做的是STREME双背景和PWM相似比较）；RegulonDB sigma专门发现/250条外部验证；代表性MEME对照；SpaMo及机制专项；其他物种参考扩展。不能称整个课程项目已完成。本轮结果为真实短元件参考注释与冻结17候选的天然背景holdout检验。
