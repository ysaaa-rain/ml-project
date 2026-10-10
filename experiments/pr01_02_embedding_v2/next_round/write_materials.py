"""Build research materials from frozen result tables, without refitting."""
import json
import pandas as pd
from ..audit import BRANCH,ROOT,sha,save
from .common import NEXT,PRESENTATION,SOURCES

def table(d,digits=4):
    def cell(x):
        return f'{x:.{digits}f}' if isinstance(x,float) else str(x)
    return '\n'.join(['| '+' | '.join(map(str,d.columns))+' |',
        '| '+' | '.join(['---']*len(d.columns))+' |']+
        ['| '+' | '.join(cell(x) for x in row)+' |' for row in d.itertuples(index=False,name=None)])

def main():
    metrics=pd.read_csv(NEXT/'synthetic_benchmark/strict_metrics.tsv',sep='\t')
    mainmask=((metrics.width==10)&~metrics.scenario.str.startswith(('S5','S6')))|((metrics.width==6)&metrics.scenario.str.startswith(('S5','S6')))
    selected=metrics[(metrics.split=='test')&(metrics.tolerance_bp==2)&mainmask]
    summary=selected.groupby(['scenario','method']).f1.agg(['mean','min','max']).reset_index()
    summary.to_csv(NEXT/'synthetic_benchmark/strict_test_summary.tsv',sep='\t',index=False)
    means=summary.pivot(index='scenario',columns='method',values='mean')[['OH','A0','A1','A2']].reset_index()
    stability=pd.read_csv(NEXT/'stability_diagnostics/tolerance_stability.tsv',sep='\t')
    selectedst=stability[(stability.strategy=='full')&(stability.pwm_similarity>=.8)&(stability.method=='prokbert')]
    st=selectedst.groupby(['stage','tolerance_bp']).agg(paired_families=('pwm_similarity','size'),position_only_median=('physical_position_jaccard','median'),direction_aware_median=('direction_aware_aligned_jaccard','median'),uniform_null_median=('uniform_null_jaccard_mean','median')).reset_index()
    st.to_csv(NEXT/'stability_diagnostics/A0_stability_summary.tsv',sep='\t',index=False)
    gc=pd.read_csv(NEXT/'gc_diagnostics/gc_distribution_overlap.tsv',sep='\t')
    counts=pd.read_csv(NEXT/'real_data_comparison/primary_seed_counts.tsv',sep='\t')
    counttable=counts.pivot(index=['source','method'],columns='control',values='significant_new_family').reset_index()
    a2=counts[(counts.method=='A2')&(counts.control=='natural')][['source','candidates','significant_new_family']]
    a2['shuffle_significant']=0
    a2st=pd.read_csv(NEXT/'real_data_comparison/A2_seed_stability.tsv',sep='\t')
    shuffle=pd.read_csv(NEXT/'gc_diagnostics/shuffle_discordance_diagnostic.tsv',sep='\t')
    cmh=pd.read_csv(NEXT/'gc_diagnostics/motif_gc_adjustment.tsv',sep='\t')
    cmhcounts=cmh.groupby('method').cmh_q_exploratory.apply(lambda x:int((x<=.05).sum())).reset_index(name='exploratory_q_le_005')
    inventory=pd.read_csv(NEXT/'run_manifests/run_inventory_final.tsv',sep='\t')
    qa=json.loads((NEXT/'verification.json').read_text())
    report=f'''# PR01-02 Embedding 后续实验研究报告

日期：2026-10-09。版本：EMB-next-1；以本报告和 `next_round/` 运行证据补充首轮 `../REPORT.md`，首轮结果未被覆盖。

## 1. 本轮回答的问题与结论范围

传统 STREME 主线从序列保守性发现候选模式；本分支考察冻结 DNA 模型的上下文表示能否提供额外发现价值。研究单位是短 DNA 模式及其可追溯实例，不是启动子分类准确率。Embedding 聚类必须返回实际 DNA 片段、构建位置权重矩阵（PWM），再在独立序列中扫描与检验，才能讨论模式发现。

本轮完成定位稳定性诊断、GC 与 Shuffle 诊断、聚合权重代数核验、六物种 A2 三种子实验、A0 数值实现对照，以及九类合成条件的难度与窗口消融。当前没有证据支持 Embedding 稳定优于 One-hot 或传统 STREME 的总体额外发现价值；也没有新增功能确认的启动子调控 Motif。

三个较明确的认识是：（1）聚类抽取的建模片段与固定 PWM 的独立扫描位点是不同对象，稳定性必须分别报告；（2）覆盖比例加权 A2 在大多数内部窗口与 A0 数学等价，不能把它描述成全面去除窗口外信息；（3）改善依赖窗口宽度、模式退化和背景，强模式的小幅改善不能外推到真实数据整体性能。

## 2. 冻结输入与实现

复用 `data/processed/pr01_02_data_v2` 的六物种 TJU 主数据及已有 leakage_group 划分；不重下载、不重划分、不增加新天然负样本。共有 9708 条正序列、8641 条天然负序列，Discovery 正序列 6745 条。真实序列长度为 81 bp，TSS 是 1-based 第 61 位（0-based 索引 60）。上游 RegulonDB 的 101-bp Gold 是另一数据层，不能与本轮 TJU 输入混称。

Discovery 只用正序列构建表示、PCA、聚类和 PWM；天然 Discovery 负序列用于冻结零阶背景与阈值。Development 用于诊断和比较，因反复使用只能视为探索性验证。Embedding 分支未使用真实 Holdout 评价；图 F01 的 Holdout 仅为冻结数量。已有主线是否曾查看 Holdout 不改变本分支的声明范围。

模型是冻结 `neuralbioinfo/prokbert-mini`，revision `feb2520a43cd9cdb5b3d8477e47209dbcb55d1dc`，远程代码 revision `366d9336b6b89d7c1feddf56d9b01be8065ccc43`。首轮加载报告没有 missing/unexpected/mismatched keys，权重 SHA256 为 `19ede01f2ef0420e111c0be8879357fba23332fad88a27571bd8d9b3b04bca3b`。本轮复用按顺序 sequence_id 与哈希绑定的 L6 token 缓存。81-bp 序列有 76 个 stride=1 的 6-mer token，每个 hidden state 为 384 维。本轮没有微调、没有 SAE、没有 GPU 实验。

运行环境与版本见各 manifest：Python 3.13.5 / Torch 2.7.1 CPU / Transformers 5.8.0。真实实验窗口 10 bp、PCA 24 维、MiniBatchKMeans 24 簇、n_init=5、种子 20261005/06/07；每簇每条启动子只取离中心最近的一个窗口，至少 20 条启动子，构建 PWM 后按全宽 RC-aware cosine 0.95 去重。One-hot 与 Embedding 共用发现与验证流程；STREME 是已有固定主线候选，候选预算、发现目标不同，不按显著候选总数给方法排优劣。

## 3. 评价口径

真实数据采用自定义 PWM 最大 log-odds 双链扫描，每条序列每个 PWM 只保留最大分位点。固定对称零阶 Discovery 天然背景；阈值是 Discovery 天然负样本最大分分布的第 95 百分位，命中条件严格为 score > threshold。它不是 FIMO 的 p/q，也不等于统计显著的单个位点。

天然富集以 leakage_group 为统计单位，组内 presence 取最大；共享混合标签组排除；Fisher 单侧检验。Shuffle 对照通过冻结的原序列—二核苷酸重排关系配对，组级 discordant pairs 做单侧二项检验。两类背景单独 BH。`q_original` 保留旧固定检验族；`q_new_exploratory` 是追加 A2 后、每个种子和每类背景分别校正的新探索检验族，不替换正式主线或首轮 q。

PWM 跨种子采用 Hungarian 一对一家族匹配，同时保留全宽匹配与最小 6 列局部匹配。核心描述以全宽 cosine≥0.8 为筛选，各阶段必须有共同 sequence_id。位置比较分别给出忽略方向的物理起点、包含 PWM 对齐方向的起点、±0/1/2/3 bp 容差与区间 IoU。相同家族在不同种子配对中重复出现，这些配对不是独立生物学重复，不据此计算正式置信区间或显著性。

## 4. 定位稳定性：先区分测量对象

以下为 A0 六物种、三个随机种子、全宽 cosine≥0.8 的家族配对中位数。构建片段有 280 个合格家族配对，Development 扫描有 263 个；两阶段分母不同。

{table(st)}

建模片段的精确位置 Jaccard 约 0.0054，允许 ±2 bp 后约 0.0304；固定 PWM 的 Development 扫描分别约 0.1642 和 0.3913。加入方向约束会改变数值，因此不能只展示最有利的定义。小幅偏移解释了部分扫描不一致，却不能解释建模片段全部不稳定。PWM 相似度高并不保证来自同一批、同一位置的 DNA 片段。

诊断同时保存 ±全范围偏移、近 TSS（距离≤15 bp）与其他位置、扫描阈值、每侧位点数、原特征空间聚类中心 cosine/L2。聚类中心比较没有假定 hidden vector 存在已知 RC 变换，因此只作描述。200 次均匀随机位置和保持边缘位置分布的经验重排基线用于参照，不当作正式 p 值。One-hot 六物种三种子与 A1 枯草芽孢杆菌三种子也有逐家族数据；A1 其他物种首轮只有单种子，不能声称完成六物种 A1 三种子稳定性。

![定位稳定性](figures/F02_stability.png)

![偏移分布](figures/F03_offsets.png)

![实际 DNA 失败示例](figures/F11_instance.png)

图 F11 按固定字典序选择首个偏移>10 bp 的真实 Discovery 片段案例，仅用于解释失败，不代表总体。

## 5. 天然对照、Shuffle 与 GC

六物种正序列的 GC 均值均低于天然负序列。下表是 sequence-level 的组成描述，不是组级显著性检验：

{table(gc)}

天然对照回答候选是否区别于现有非启动子样本；二核苷酸 Shuffle 更接近询问超出原有二核苷酸组成后是否仍有模式效应。两者效应量和显著性不同不构成矛盾，也不能归结为“Shuffle 太严格”或“全部都是 GC 假象”。

按组级均值 GC 固定分层 [0,.2,.4,.6,.8,1]，每层至少 5 个正组和 5 个天然负组，计算探索性 CMH pooled OR/p，并单独 BH。保留组覆盖比例、分层效应与未校正结果。仍有少数候选保留关联：

{table(cmhcounts)}

CMH 只控制粗粒度 GC 分层，不能控制二核苷酸结构、局部背景、物种内部来源偏差或反复选择；没有“GC 已解决”的结论。

旧主种子 Shuffle 族有 {len(shuffle)} 个检验，discordant groups 中位数 {shuffle.discordant_groups.median():.0f}。其中 {int(shuffle.discordance_too_small_for_first_rank_bh.sum())} 个候选即使所有 discordant groups 都偏向正序列，条件下限 p=2^(-D) 仍大于首位 BH 阈值 .05/{len(shuffle)}。这只是可达到 p 值的条件边界；BH 后续排名的阈值不同，不能宣称这些候选绝无通过 BH 的可能，也不能把它称为事后检验功效或证明低功效解释了所有不显著。表中保存正独有 b、Shuffle 独有 c、(b-c)/N 效应量。

![GC 组成](figures/F04_gc.png)

![不同背景的效应量](figures/F05_controls.png)

## 6. Pooling 改进及其代数诊断

设 token t 覆盖碱基区间 I_t，目标窗口 W。A0 先把 token 平均投射到每个碱基，再对 W 平均；其权重为 w0(t)=Σ[b∈W∩I_t] 1/(|W|·c_b)，c_b 为覆盖碱基 b 的 token 数。A1 只平均完整落在窗口内的 token。A2 按覆盖比例 |W∩I_t|/6 加权，再归一化。

当 W 中每个碱基均由六个 token 覆盖，w0(t)=|W∩I_t|/(6|W|)，与归一化的 A2 完全一致。81-bp 输入、10-bp 窗口共有 72 个位置，其中 62 个内部位置等价，只有两端各 5 个位置有真正权重差异。宽度 6 和 16 也通过代数与单元检查。

A1 限制的是参与聚合的 token 区间；上下文化 hidden state 仍可含窗口外信息。A2 同样不能普遍隔绝上下文。因此“定位不稳由窗口外信息导致”仍是假设，现有实验没有建立单一因果解释。边缘权重变化可以通过整体 PCA/聚类改变内部候选，不应把候选差异直接解释成局部窗口语义改善。

新增六物种 A0_recalc 主种子对照使用与 A2 相同的 einsum 路径重新实现 A0。六物种候选数均与原 A0 一致，匹配 PWM 相似度中位数均为 1.0；这降低了本组真实结果仅由实现路径数值差异造成的疑虑，不能作为生物学提升的证据。

![聚合权重原理](figures/F06_pooling.png)

## 7. 合成难度梯度与独立 Test

九类条件：强模式固定位置；15%/30% 逐位替换；强模式随机位置；GC=.25/.75；两个 6-bp 模式固定 gap=17；gap=10–24；一阶持续性=.7 的 Markov 背景。每类三个种子 20261010/11/12，共 27 个数据条件。每个条件独立生成 Discovery 256、Development 128、Test 128 条正序列及对应空白对照，另有 256 条阈值校准对照。方法共享同一份序列；无反向链植入难度梯度。背景 GC 是生成参数，植入后的实际全序列 GC 可能变化。

单模式主窗口 10 bp、双模式 6 bp；随机强模式额外测 6/16 bp，四方法共 132 次拟合。相同簇预算 24、最低 20 条序列、每簇每序列一个建模窗口；合成实验不做真实数据的候选去重，故候选数量不能直接与真实数据比。

真值位置不进入 PCA/聚类/片段选择/PWM 构建。评价时用合成参考 PWM 与候选一对一匹配，cosine≥0.8 才计入恢复；这属于参考辅助的模式恢复评价，不能解释成无人监督发现的绝对 precision。参考 alignment 用于把短窗口扫描位点投射为植入模式起点；6-bp 窗口匹配 12-bp 模式只证明部分锚点定位，不证明完整 12-bp 边界恢复。

初始评价允许最小 6 列局部相似；在汇总指标前修正为至少覆盖 min(候选宽度,参考宽度)，避免把局部相似当作长模式恢复。所有方法统一重新评价已冻结 PWM，不重新拟合；旧宽松输出与初始源文件保留在 configs，严格规则见 `configs/strict_evaluation.json` 和 manifest。后续若利用这些 Test 结果改方法，应另用新种子作为独立 Test，不能继续把本 Test 当未查看的最终集。

TP 要求预测方向正确且起点误差不超过指定容差；precision 的分母包含植入正序列上的全部阳性预测及空白对照上的预测，位置错误记 FP 且对应真值未找回，recall 分母为所有植入位点。每个家族每条序列最多一个最大分预测。PWM 未过门槛记为 0；这不是 hidden state 完全没有相关信息的证明。

独立 Test、±2 bp、三个种子（双模式再对两家族平均）的 F1：

{table(means,3)}

随机强模式、10-bp 窗口中 A2 平均 F1 约 .992，A0 约 .970，One-hot 约 .981；但退化模式 A0/A1/A2 均未过严格 PWM 门槛，One-hot 在 15% 与 30% 替换条件约 .840/.448。16-bp 随机强模式 A0/A2 为 0，而 A1 约 .624、One-hot 约 .955。GC=.75 和 Markov 条件部分 Embedding 种子失败，均值必须结合 min/max 查看。只用一个强模式场景讲“改进成功”不成立。

![严格 Test 难度矩阵](figures/F07_synthetic.png)

![窗口消融](figures/F08_width.png)

双模式 grammar 的总体分母为每方法每条件 384 条独立 Test 启动子。图 F09 只画双模式同时检出的点；漏检者仍在 `grammar_summary.tsv` 的总体正确顺序率和 gap±3 率分母内。固定 gap 条件 A0 未同时检出两模式；变化 gap 的较高描述性恢复不能证明真实启动子中存在相同调控机制，也不能假定固定 gap 一定更容易。

{table(pd.read_csv(NEXT/'synthetic_benchmark/grammar_summary.tsv',sep='\t'))}

![双模式间距](figures/F09_pair.png)

本轮合成矩阵没有重新执行 STREME；因此只能比较四种表示构造流程。真实数据 STREME 对照来自冻结主线。这是总体方法排名仍未闭合的一个限制。首轮 toy 的“合格候选数/纯度”和本轮的参考辅助定位 F1、数据种子不同，不能直接合并成单一提升曲线。

## 8. 六物种真实数据新增 A2 结果

新增 A2 六物种×三种子 18 次真实实验，另六次 A0_recalc 数值控制。主种子在新探索检验族下：

{table(a2)}

自然富集候选数不能直接评价额外发现价值，尤其候选间相互依赖且背景组成存在差异。A2 三种子全宽一对一 PWM 相似度中位数 {a2st.full_pwm_similarity.median():.4f}，建模片段精确位置 Jaccard 中位数 {a2st.instance_exact_jaccard.median():.4f}；该汇总包含所有匹配家族，未套用前文 A0 cosine≥.8 筛选，因此不能据数值大小声称 A2 定位改善。

全部方法主种子新探索 q≤.05 的候选数（天然与 Shuffle 独立校正；只作描述）：

{table(counttable)}

对六物种 A2，各种子 Shuffle 均未得到 q≤.05 的候选。结合 One-hot、传统 STREME、定位稳定性和合成退化实验，当前只能说 A2 生成了可追溯 PWM，未证明可靠的额外发现价值。没有从真实 Development 得出召回率或准确率，因为不存在完整的真实位点真值。

## 9. 传统基线与生物学解释边界

传统主线已有双背景 STREME 12 个运行、17 个天然背景冻结主候选；后续严格 FIMO/参考库/grammar 报告优先于更早“下一步待做”的文案。本轮另用同一自定义最大分口径绘制 Development 的传统 PWM TSS 位置与配对 gap，仅为描述，不能称为严格 FIMO q 通过、−35/−10 调控确认或显著共现。

旧参考库的 4 个短参考矩阵文件当前缺失，只有元数据与历史报告可用。没有重新制造矩阵，也没有用大肠杆菌共识冒充六物种已知真值；新的数据库匹配结论暂缺。原主线严格检验没有支持显著 presence/cooccurrence 或可确认的位置间距机制，本轮图示不能推翻这一结论。

![传统候选位置与间距](figures/F10_grammar.png)

![原 STREME 候选](figures/F12_streme_logo.png)

## 10. 可复现性、运行资源与核查

`protocol.json`、各 manifest、旧源文件快照、方法代码、逐实例与扫描表、参考匹配表均保留。manifest 记录真实输入/输出 SHA256、代码 SHA256、Git HEAD、依赖、CPU、开始时间、耗时与进程内存高水位。内存是顺序进程累计峰值，不能当作单个方法隔离峰值；真实 A2 不含首轮模型编码成本，不能拿这些秒数与端到端 STREME 时间直接比较。各记录耗时总和约 {inventory.seconds.sum():.1f} 秒，亦不是本任务从开始到交付的墙钟时间。

复核结果：{qa['hashes_checked']} 项输入/输出哈希通过，53 个首轮运行输出未改；24 个新真实运行的 {qa['new_real_instances']} 条片段与坐标通过；132 个合成拟合的 {qa['synthetic_instances']} 条片段及 {qa['synthetic_truth_sites']} 个植入真值通过；27 个合成条件跨 Discovery/Development/Test 与控制组无完全相同序列。7 项方法测试通过。哈希核查证实保存文件一致，不能独立重建每一次历史执行；没有新跑全量同源比对。

图表全部来自保存表格或明确标记的权重原理图；12 张已进行实际图片检查，F07 标签重叠已修正。图索引带源路径、统计单位、局限和图像哈希。

可复核命令（项目根目录）：

```powershell
python -m experiments.pr01_02_embedding_v2.next_round.verify_results
python -m pytest experiments/pr01_02_embedding_v2/test_pipeline.py experiments/pr01_02_embedding_v2/next_round/test_methods.py -q
python -m experiments.pr01_02_embedding_v2.next_round.figures
python -m experiments.pr01_02_embedding_v2.next_round.write_materials
```

原始拟合入口是 `next_round.real`、`next_round.benchmark`；它们拒绝覆盖已存在的 manifest/目录。重做实验应先复制代码到新的版本目录、指定新输出位置；不要删除当前证据直接重跑。仅绘图/汇总脚本允许重新生成派生材料。最终 inventory 为 `run_inventory_final.tsv`；早期 `run_inventory.tsv` 是汇总程序完成前快照，不能据其 RUNNING 判断实际状态。

## 11. 失败记录、未完成项与后续判据

已观察的限制包括：建模片段不稳定；严格退化模式未恢复；A2 主要内部等价；高 GC/相关背景下种子失败；双模式漏检；Shuffle 无可靠通过；真实位点真值不全；缺失参考矩阵；真实 Development 已反复查看。

A3 中心加权仅有实现定义，未运行；A4 局部片段重新编码、SAE、微调、反向链植入梯度、更多独立种子与新的 STREME 合成矩阵均未执行，不能写为成果。首轮部分状态仍为 PARTIAL；本轮完成的是这组诊断与比较，项目整体发现结论仍是部分完成。

后续优先做局部重新编码与全序列上下文的成对控制，保持窗口、片段选择和评价统一，以测试上下文来源假设；其次把初始表示选择与聚类选择分开诊断，并用新合成 Test 种子复核退化、宽度和种子失败。若要论证真实额外发现，应先恢复物种可靠参考矩阵/位点来源、预先锁定候选与独立验证族，再使用未查看的合适验证资源。增加模型数量不是充分的推进理由。

## 12. 两分钟汇报的准备边界

依据最新要求，主体只保留三页、总时长≤120秒。当前先交付完整研究基础、图表、问答与结论证据表，不把技术细节塞进课堂主体。之后三页压缩为：为什么尝试上下文表示；实际候选与定位/背景问题；已完成聚合诊断与条件性结果。旧“计划开展覆盖加权”应替换为本轮实测与等价性认识，结尾回到可靠启动子 Motif 的发现与验证。
'''
    (NEXT/'REPORT.md').write_text(report,encoding='utf-8')
    backup='''# 备用材料与教师问答（不进入三页主体）

1. **为什么尝试 Embedding？** 传统方法从序列保守性提取模式。预训练模型提供上下文表示，值得检验是否能补充；这是待验证假设，并非预设深度学习更优。
2. **Embedding 和 One-hot 差在哪？** One-hot 直接编码碱基身份；Embedding 来自冻结模型并含上下文。共用聚类和 PWM 流程，方便检查表示本身的适用性。
3. **聚类就是 Motif 吗？** 不是。需要找到实际 DNA 片段、构建 PWM，并在独立序列及背景上验证；功能解释还需要外部位点与生物学依据。
4. **PWM 相似为何位置不稳？** 相似矩阵可由不同位置、不同实例构成。随机种子与簇选择可能改变建模片段；固定 PWM 独立扫描的稳定性是另一指标。
5. **容差是否掩盖问题？** 同时保留精确与 ±1/2/3 bp、方向与位置两种定义。扫描结果对小偏移敏感，建模片段即使容差后仍不稳定。
6. **天然对照显著，Shuffle 不显著，为何？** 两类对照检验不同背景假设；组成、局部序列结构和配对样本量都可能影响。不能选择对自己有利的一类报告。
7. **GC 混杂解决了吗？** 没有。已画六物种分布、保存分层效应并做探索性 CMH。粗 GC 分层不能控制所有序列背景和来源差异。
8. **Shuffle 检验是不是功效不足？** 已保存 discordant counts 和效应量，可达到 p 下限只是条件描述，不能证明全部不显著都源于功效，也不能把首位 BH 阈值当全部候选的阈值。
9. **覆盖加权确实改善了局部表示吗？** 62/72 个内部位置与原聚合数学等价；差异主要在输入边缘。候选变化可能经整体 PCA/聚类传播，不能宣称普遍消除了窗口外信息。
10. **窗口内 token 就没有上下文了吗？** 不是。参与聚合的 token 可在窗口内，但 hidden state 仍由全序列上下文化产生。局部重新编码尚未运行。
11. **合成 F1=0 表示模型完全无信息吗？** 不是。这是 PWM 参考匹配与定位联合评价，未过 0.8 匹配门槛记零。它显示当前发现流程未恢复合格模式。
12. **合成结果可证明总体提升吗？** 随机强模式某些宽度改善，退化、宽度16及某些背景失败。三个种子范围已保留，无总体提升结论。
13. **用了真值是不是泄漏？** 真值不进入拟合或片段选择，仅用于合成参考家族映射与定位评价。该指标是参考辅助恢复能力；本轮 Test 已查看，后续改方法应使用新的独立种子。
14. **为什么没有新的真实召回率？** 真实六物种没有完整可靠位点真值；只报告独立扫描、背景富集和稳定性，不能虚构召回率。
15. **是否已经比 STREME 更好？** 没有证据。真实数据有冻结 STREME 对照，但候选预算不同；本轮合成未运行 STREME，不能给总体排名。
16. **这算算法创新吗？** 目前是表示构造的可检验尝试和适用性诊断。没有宣称新的网络架构、首创算法或生物学发现。
17. **为什么没有做 SAE？** 本轮优先厘清表示、定位与背景问题。SAE 未运行，没有把它当已完成成果。
18. **双模式间距图为什么漂亮但恢复率有限？** 散点只显示双检出者；384 条 Test 的总体分母含漏检，必须同时看 grammar_summary，不能只看对角线。
19. **真实 −35/−10 机制得到确认了吗？** 没有。当前传统图是自定义最大分扫描的描述性位置/间距，不等于严格 FIMO q 或机制证据。
20. **下一步最值得做什么？** 成对比较局部重新编码与全序列表示，用新种子验证退化和宽度问题；先证明可靠的位点恢复，再追求额外发现。

## 查图入口

- F01：固定六物种正序列划分，Holdout 仅数量。
- F02/F03/F11：定位稳定性、偏移分布、真实失败实例。
- F04/F05：GC 分布与两类背景效应。
- F06：聚合权重原理；内部 A0=A2。
- F07/F08：严格合成难度矩阵、宽度消融；零值和三种子范围需解释。
- F09：双模式间距，另表保留总体漏检分母。
- F10/F12：传统冻结候选及描述性位置，不能解释成功能确认。

所有图在 `next_round/figures/`；实际源文件与图像哈希见 `evidence/figure_index.csv`。完整数字和限制见 `../next_round/REPORT.md`，逐指标数据不进入两分钟主体。
'''
    (PRESENTATION/'06_backup_materials.md').write_text(backup,encoding='utf-8')
    claims=[
      ('C01','Embedding能生成可追溯PWM，未证明总体额外发现价值','real_data_comparison/all_methods_enrichment.tsv','real_A2_18_runs','OH/A0/A1/A2/STREME','Development','background enrichment','观察与有限结论','无完整真实位点真值；不同候选预算',2),
      ('C02','A0建模片段精确/±2bp位置Jaccard中位数约.0054/.0304','stability_diagnostics/A0_stability_summary.tsv','historical_stability','A0','Discovery','position-only Jaccard median','已运行','全宽cos≥.8；280家族配对；忽略方向',2),
      ('C03','A0独立扫描精确/±2bp位置Jaccard中位数约.1642/.3913','stability_diagnostics/A0_stability_summary.tsv','historical_stability','A0','Development','position-only Jaccard median','已运行','263家族配对，与建模片段分母不同',2),
      ('C04','天然与Shuffle效应不同，GC组成需考虑','gc_diagnostics/gc_distribution_overlap.tsv','historical_gc_and_shuffle','all','Discovery/Development','GC mean/distribution','已运行；解释为假设','不能证明GC导致所有差异或已消除混杂',2),
      ('C05','10bp窗口72位置中62个内部位置A0=A2','pooling_ablation/weight_equivalence.tsv','pooling_algebra','A0/A2','analytic','normalized token weights','代数核验完成','6mer stride1、81bp；hidden state仍含上下文',3),
      ('C06','强随机10bp合成Test定位F1 A0约.970，A2约.992','synthetic_benchmark/strict_test_summary.tsv','synthetic_strict_evaluation','A0/A2','synthetic Test','reference-gated ±2bp F1 mean3seeds','已运行','cos≥.8；非整体提升；未比较合成STREME',3),
      ('C07','15%退化条件One-hot约.840，Embedding未过严格PWM门槛','synthetic_benchmark/strict_test_summary.tsv','synthetic_strict_evaluation','OH/A0/A1/A2','synthetic Test','reference-gated ±2bp F1','已运行','零值不等于hidden state无信息',3),
      ('C08','真实A2各六物种三种子Shuffle无q≤.05候选','real_data_comparison/all_methods_enrichment.tsv','second_round_collection','A2','Development','q_new_exploratory','已运行','新的探索检验族；非FIMO单个位点q',3),
      ('C09','局部聚合是否混入上下文需进一步成对控制','pooling_ablation/weight_equivalence.tsv','pooling_algebra','A0/A1/A2','analytic','hypothesis','未建立因果解释','局部重新编码未运行',3),
      ('C10','A0数值重实现对照六物种PWM中位相似度均1','real_data_comparison/A0_numerical_control.tsv','A0_recalc_6_runs','A0','Discovery','matched PWM cosine','已运行','只控制实现路径，不是生物学证据','backup'),
      ('C11','双模式间距恢复需要包含未检出者','synthetic_benchmark/grammar_summary.tsv','synthetic_strict_evaluation','OH/A0/A1/A2','synthetic Test','gap/order rate all384','已运行','双检出散点不能代替总体恢复率','backup'),
      ('C12','CMH只是粗GC分层探索，未解决全部混杂','gc_diagnostics/motif_gc_adjustment.tsv','historical_gc_and_shuffle','all','Development','CMH OR/q separate family','已运行','分层覆盖与来源差异仍有局限','backup'),
      ('C13','新增24真实运行、132合成拟合；真实Holdout未评价','verification.json','verification','all','Discovery/Development/synthetic Test','artifact audit','已核查','不代表所有后续方法完成','backup')
    ]
    fields=['claim_id','text','source_file','run_id','method','split','metric','evidence_status','limitation','planned_slide']
    frame=pd.DataFrame(claims,columns=fields)
    frame['source_file']=frame.source_file.map(lambda p:str((NEXT/p).relative_to(ROOT)))
    frame['source_sha256']=frame.source_file.map(lambda p:sha(ROOT/p))
    frame.to_csv(PRESENTATION/'07_claim_evidence_matrix.csv',index=False,encoding='utf-8-sig')
    save(NEXT/'STATUS.json',{'version':'EMB-next-1','research_round':'COMPLETED_WITH_LIMITS','overall_motif_discovery':'PARTIAL',
      'completed':['localization diagnosis','GC/Shuffle diagnosis','pooling algebra','18 A2 real runs','6 A0 numerical controls','27 synthetic conditions / 132 fits','strict reevaluation','research report','12 figures','backup QA','claim evidence matrix'],
      'not_executed':['A3 fitting','A4 local reencoding','SAE','fine tuning','new synthetic STREME','real Holdout evaluation'],
      'presentation':{'status':'DEFERRED_UNTIL_RESEARCH_FOUNDATION_REVIEW','slides':3,'max_seconds':120,'script_hanzi':[250,320]},
      'holdout_used':False,'verification':qa})
    print('Written research report, backup materials, claim matrix, and status.')

if __name__=='__main__':main()
