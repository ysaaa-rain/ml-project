# PR01-02 DNA Embedding 独立分支：首轮执行报告

日期：2026-10-09。总体状态 **PARTIAL**。已完成 53 次真实数据运行、2 次合成发现实验，首轮工程闭环可复现；尚未完成全部研究问题，也未确认新的功能 Motif。

## 本轮结论和停止扩展依据

Frozen ProkBERT-mini 可以形成可追溯的聚类/PWM候选，但当前 **Token→Base→Window** 方案没有显示超越 One-hot 或传统 STREME 的可靠额外 Motif 发现价值。它在天然对照上的富集较多、PWM跨种子相似度较高，却没有通过配对 Shuffle 的跨物种 BH 校正，而且精确实例定位稳定性很低。天然对照的组成差异是合理替代解释；这些事实不能证明该模型完全无用或Motif没有生物学功能。

合成恢复诊断显示，改用“窗口内完整token平均”能改善一个强植入模式的恢复。这支持继续研究表示聚合边界，尚不能支持生物学优越性或方法创新已经成立。

本轮按任务的失败回退规则，完成基线、聚合诊断与阴性结果归档后，**暂不启动SAE、继续预训练或微调**。基础候选尚未通过稳健性门槛，未开启Holdout，也未用辅助来源拼接外部显著性。后续研发应先改善边界/实例对齐，采用独立合成难度梯度，再决定扩展；这些属于尚未执行的工作。

## 原项目事实与新增工作

原仓库 HEAD：`55d0c96797c14ff9b53c8f47ee568d486192c67a`。原STREME 12组运行和17个固定主候选对应当前v2清单。原17候选的FIMO/参考后续已有独立记录；本轮仅复用PWM，未重跑或改写传统主线。扩展M3的本地运行manifest未找到，因此不宣称42组扩展已经完整执行。M1/M2三来源待补统计维持原状态。

已保留工作区原先三份参考说明文件的删除改动。全部新增内容在本目录；没有commit/push，没有重建、重拆分或改写主线文件。

## 数据、环境和模型验收

- 原始/派生文件 **1,305/1,305** SHA256一致；数据manifest `449dca1e3ba5efb3d49b4427515b9cb953fff4e1061ef311c081ddf681927608`。
- TJU正类9,708、天然对照8,641、Discovery正类6,745；关联组跨split为0。每个主/辅助层的正类及天然对照FASTA核对sequence master，Shuffle父ID和二核苷酸组成全量核验通过。
- 81bp、TSS第61位/Python索引60，采用PromLoop已发布方向，不再RC；未逐条重建TJU原始基因组坐标。两方向≤8编辑关联隔离是操作性规则，不能等同于生物学绝对独立。
- RegulonDB core 1,793、extended 2,020；DBTBS core 432。保留多σ关联及独立证据字段；外部候选按冻结关联图核查与TJU全部split无共享组，RegulonDB core 250、DBTBS 1。没有新做跨来源统计，也没有重新计算全对全比对。
- Python 3.13.5、Torch 2.7.1 CPU、Transformers 5.8.0；当前CUDA不可用。实际依赖固定在 `requirements.txt` 和 `audit/environment.json`。
- 模型 [`neuralbioinfo/prokbert-mini`](https://huggingface.co/neuralbioinfo/prokbert-mini/tree/feb2520a43cd9cdb5b3d8477e47209dbcb55d1dc)，版本 `feb2520a43cd9cdb5b3d8477e47209dbcb55d1dc`；实现代码版本 `366d9336b6b89d7c1feddf56d9b01be8065ccc43`。使用预训练Encoder，未采用promoter分类微调模型。
- 权重82,581,348字节、SHA256 `19ede01f2ef0420e111c0be8879357fba23332fad88a27571bd8d9b3b04bca3b`，与发布方LFS摘要相同。Encoder参数20,489,472；missing/unexpected/mismatched keys均为空，没有随机初始化回退。模型eval并冻结梯度。
- 两份参考PDF共31页已读取并检查首页；综述有34处未解析引用占位符，不能据此认定论文结论已经核实。[官方tokenizer说明](https://prokbert.readthedocs.io/en/latest/prokbert_tokenizer.html)用于核对6mer/shift=1映射。

## 实验设置与评价边界

六物种分别拟合，输入各自完整Discovery正类；PCA=24、MiniBatch K-Means K=24，初始L6/W10、种子20261005/06/07。发现无Motif标签，但启动子正类入口是标签筛选后的集合；不能称对全体无标签序列盲发现。传统STREME使用差异背景发现，任务信息与仅正类聚类并不完全相同。

token i对应0-based `[i,i+6)`；每次真实编码逐条检查官方词表ID与76个明确6mer一致，CLS/SEP排除。缓存覆盖6,745条Discovery正类，层3/6、窗口6/8/10/12/16。每簇每启动子只保留最近中心窗口；同长正向对齐，完整长度两方向cosine≥0.95去重。未做gapped alignment和精确可变边界修剪。

用固定PWM对Development正类、天然label0和对应Shuffle独立扫描。背景/阈值只用Discovery天然对照，扫描两链最大log-odds，经验95%分位点、严格大于阈值。**这是自定义序列级扫描，不是FIMO site p/q，也不复用原FIMO q作为这里的证据。**

天然对照按泄漏组presence聚合、共有正负组剔除后单侧Fisher；Shuffle按父组做配对单侧二项检验。初始每run有自己的BH；最终固定L6/W10在每种子、每类对照跨六物种与所有实际方法统一BH，STREME重复记录只计一次。参数/聚合诊断使用Development，所以全部结果为探索性；没有把这些数据再当最终盲验证。

## 六物种实测结果

下表为主种子20261005、统一家族校正后的**天然对照显著候选数**，不是功能Motif数量：

| 物种 | STREME | One-hot | Token→Base→Window | 窗口内Token |
|---|---|---|---|---|
| bacillus_subtilis | 3 | 3 | 10 | 4 |
| baumannii | 2 | 6 | 8 | 5 |
| bradyrhizobium | 2 | 0 | 4 | 7 |
| diphtheria | 3 | 2 | 10 | 10 |
| escherichia_coli | 2 | 2 | 6 | 8 |
| staphylococcus | 2 | 14 | 1 | 1 |

同一统一家族中，全部方法的Shuffle显著候选总数为 **0**。候选总数及p/q、OR/区间、共有组排除数、配对discordance完整保存于 `summary/all_development_enrichment.tsv`。一个主启动子的多个窗口/位点不作为独立生物学样本。

Development正类平均GC从0.2521–0.5972，天然对照从0.3414–0.6501；每个物种正类GC均低于天然对照。分物种分标签的实际计数见 `audit/counts_gc.tsv`，固定GC区间下命中率见 `summary/gc_composition_diagnostic.tsv`。这里没有新增GC-matched背景，也没有按GC聚类；诊断不能证明消除了组成偏差。

## 层级、多尺度与定位稳定性

枯草芽孢杆菌另跑L3/W10、L6五尺度（W10在基线中）、五尺度同中心拼接融合三个种子，以及同长W16 One-hot。多尺度表示是实际拼接并重新拟合PCA，不只是多存缓存。独立局部模型编码/局部与上下文双视图融合仍未实现；没有运行DNABERT-2。

| method | layer | width | seed | motifs | seconds |
|---|---|---|---|---|---|
| multiscale | 6 | 16 | 20261005 | 23 | 19.01 |
| multiscale | 6 | 16 | 20261006 | 22 | 18.0 |
| multiscale | 6 | 16 | 20261007 | 21 | 17.49 |
| onehot | 6 | 16 | 20261005 | 24 | 25.79 |
| prokbert | 3 | 10 | 20261005 | 22 | 16.92 |
| prokbert | 6 | 12 | 20261005 | 21 | 24.96 |
| prokbert | 6 | 16 | 20261005 | 23 | 18.34 |
| prokbert | 6 | 6 | 20261005 | 21 | 13.37 |
| prokbert | 6 | 8 | 20261005 | 16 | 11.09 |

跨种子一对一PWM匹配，full-width cosine≥0.8计为PWM稳定；位点Jaccard是精确(sequence_id,start0)集合重叠：

| method | matches | stable_fraction | median_cosine | median_instance_jaccard |
|---|---|---|---|---|
| multiscale | 64 | 0.9062 | 0.9118 | 0.0398 |
| onehot | 432 | 0.169 | 0.6683 | 0.0321 |
| prokbert | 323 | 0.8669 | 0.9159 | 0.0051 |
| tokenwindow | 64 | 0.7969 | 0.8935 | 0.13 |

ProkBERT较高PWM相似和很低位点重叠并存，提示部分候选只是常见组成模式的稳定聚类；不能把“聚类稳定”当作精确定位稳定。tokenwindow只有枯草芽孢杆菌完成三种子，其他五物种仅主种子，不伪称全量三种子消融。

## 合成发现恢复与失败诊断

每次Toy使用256条Discovery、128条Development，均匀随机81bp序列在0-based25插入强12bp `TTGACATATAAT`，另有独立随机校准/验证对照。发现阶段仍仅拟合无Motif标签的窗口表示/PCA/K-Means；真值只用于事后评价，不用于挑实例或建PWM。

恢复标准：至少10bp对齐cosine≥0.8；抽取实例与植入区间重叠≥8bp比例≥0.5；独立Development富集BH q≤0.05。原始聚合诊断记录在 `synthetic/`；新增tokenwindow诊断在 `synthetic_tokenwindow/`，保留全部阴性候选。

| method | candidates | recovered |
|---|---|---|
| onehot | 24 | 3 |
| prokbert | 24 | 0 |
| tokenwindow | 24 | 1 |

这些计数仅描述一个强模式Toy，不是全难度benchmark或功能验证。tokenwindow改善了这个例子，不能证明普遍优于One-hot。原ProkBERT聚合把覆盖窗口外侧的token也平均进来，可能模糊局部序列模式；该解释是结合消融的推断，仍需多植入位置/保守性/GC/种子实验。

## Grammar、参考与外部验证的实际状态

每次真实运行有Development主位点center/strand表与同序列两Motif的顺序、center distance、edge gap（负值为重叠）。目前是描述性产物，未完成聚类单位置换检验、可靠σ组比较或功能机制归属。

全部Embedding候选与同物种固定STREME PWM做两方向offset余弦比较，保存完整描述表，未将它当Tomtom显著性。原报告的4个已知短元件PWM仅存在历史本地tmp路径；当前工作树未找到矩阵，故本轮没有声称完成数据库PWM对照，公开来源元数据仍可追溯。其余四物种无合格参考的原限制保留。

RegulonDB外部统计、DBTBS一条外部留出、最终Holdout均未运行。传统主线的Holdout已被查看，即使未来锁定Embedding配置也要说明属于探索性再评价。预训练语料与目标基因组的重叠未审计。

## 工程验证、成本与失败记录

- 4项自动测试通过，覆盖token/base边缘、TSS映射、拒绝Holdout、正反链合成扫描、PWM/BH和真实Shuffle组成。
- 53 次运行输出哈希与 **789,842 条抽取实例** 逐条校验原始DNA子串、group/split、0/1-based/TSS坐标和PWM归一化，通过。实例数包含不同方法/种子，绝不能解释为独立样本规模。
- Logo单图及六物种示例联系图已目视检查； `summary/logo_contact_sheet.png` 仅示例，并非筛选后的有效Motif集合。
- 真实数据运行累计计时 1277.26 秒（约21.29分钟），并发执行，非端到端墙钟时间；下载、Toy、汇总/QA不在此累计内。缓存复用种子运行更便宜。分支落盘约11.86 GB，主要是两层五尺度缓存；GPU时间/显存不适用，CPU峰值内存未测。
- 首次默认权重下载停滞，改用固定revision的HTTP Range分块下载，校验Content-Range、完整长度和LFS SHA256。首次本地import缺少sys.modules登记导致官方代码的docstring检查报错，已修复；最终权重加载无缺失。第一次命令默认GBK读取UTF-8 JSON失败，现已显式UTF-8。失败尝试没有生成伪模型结果。

## 工作包与复现

`STATUS.json`逐项列明完成范围。EMB-01/02及基础EMB-03工程闭环完成；EMB-04/05/07为PARTIAL；EMB-06为PLANNED、因前序证据不足未启动。潜在多尺度双视图、Grammar感知排序、边界敏感SAE仍是待检验假设，没有写成创新成果。

配置与验收判据在 `README.md`/`config.json`，真实依赖在 `requirements.txt`。仓库根可执行：

```powershell
python -m experiments.pr01_02_embedding_v2.audit
python -m pytest experiments/pr01_02_embedding_v2/test_pipeline.py -q
python -m experiments.pr01_02_embedding_v2.run --sources tjupan_bacillus_subtilis tjupan_baumannii tjupan_bradyrhizobium tjupan_diphtheria tjupan_escherichia_coli tjupan_staphylococcus --representations prokbert onehot --seeds 20261005 20261006 20261007
python -m experiments.pr01_02_embedding_v2.ablate
python -m experiments.pr01_02_embedding_v2.run --sources tjupan_bacillus_subtilis tjupan_baumannii tjupan_bradyrhizobium tjupan_diphtheria tjupan_escherichia_coli tjupan_staphylococcus --representations tokenwindow --seeds 20261005
python -m experiments.pr01_02_embedding_v2.run --sources tjupan_bacillus_subtilis --representations tokenwindow --seeds 20261006 20261007
python -m experiments.pr01_02_embedding_v2.synthetic --tag synthetic_tokenwindow
python -m experiments.pr01_02_embedding_v2.summarize
python -m experiments.pr01_02_embedding_v2.verify
python -m experiments.pr01_02_embedding_v2.report
```

同名运行拒绝覆盖。复现需在新的项目内checkout或分支目录开展；不要删除现存结果后假装原运行历史不存在。每run存代码文件哈希、Git HEAD、输入/输出哈希、参数、seed、计时和实际状态。模型/特征缓存从Git忽略；可审阅的小型代码、报告及汇总留在独立目录。
