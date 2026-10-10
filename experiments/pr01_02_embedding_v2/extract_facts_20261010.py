"""Read existing artifacts and extract report facts. No fitting or scanning."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
from .discovery import read_meme

BRANCH=Path(__file__).resolve().parent
ROOT=BRANCH.parents[1]
NEXT=BRANCH/'next_round'
OUT=BRANCH/'facts_20261010'

def clean(value):
    if isinstance(value,dict):return {k:clean(v) for k,v in value.items()}
    if isinstance(value,list):return [clean(x) for x in value]
    if isinstance(value,np.generic):return clean(value.item())
    if isinstance(value,float) and not np.isfinite(value):return None
    return value

def distribution(values):
    a=np.asarray(values,float)
    if not len(a):return {'n':0,'summary':'no observed hits'}
    return {'n':len(a),'min':a.min(),'q25':np.quantile(a,.25),'median':np.median(a),
            'q75':np.quantile(a,.75),'max':a.max(),'mean':a.mean()}

def link(path,label=None):
    p=Path(path)
    if not p.is_absolute():p=ROOT/p
    return f'[{label or p.name}]({p.as_posix()})'

def mdtable(headers,rows):
    return '\n'.join(['| '+' | '.join(headers)+' |','| '+' | '.join(['---']*len(headers))+' |']+
      ['| '+' | '.join(str(x) for x in row)+' |' for row in rows])

def number(x):
    return f'{x:.6g}'

def report(r):
    b='experiments/pr01_02_embedding_v2/'
    n=b+'next_round/'
    counts=pd.DataFrame(r['new_counts'])
    methods=['OH','A0','A1','A2','STREME']
    names=['枯草芽孢杆菌','鲍曼不动杆菌','慢生根瘤菌','白喉杆菌','大肠杆菌','葡萄球菌']
    sources=['tjupan_bacillus_subtilis','tjupan_baumannii','tjupan_bradyrhizobium','tjupan_diphtheria','tjupan_escherichia_coli','tjupan_staphylococcus']
    rows=[]
    for source,name in zip(sources,names):
        cells=[]
        for method in methods:
            a=counts[(counts.source==source)&(counts.method==method)&(counts.control=='natural')].iloc[0]
            s=counts[(counts.source==source)&(counts.method==method)&(counts.control=='shuffle')].iloc[0]
            cells.append(f'{a.candidates}/{a.significant_new_family}/{s.significant_new_family}')
        rows.append([name]+cells)
    totals=[]
    for method in methods:
        a=counts[(counts.method==method)&(counts.control=='natural')]
        s=counts[(counts.method==method)&(counts.control=='shuffle')]
        totals.append(f'{a.candidates.sum()}/{a.significant_new_family.sum()}/{s.significant_new_family.sum()}')
    rows.append(['合计']+totals)
    text=f'''# PR01-02 Embedding 分支：实验事实汇总

日期：2026-10-10。仅从已有代码、报告和运行表读取、汇总，未新增训练、扫描或统计检验，未修改旧产物。下文 A0=`prokbert`，A1=`tokenwindow`，OH=`onehot`。

**核实结论：A2 主种子共116个保留PWM，45个天然对照显著、0个Shuffle显著，准确。** 范围为 EMB-next-1、六物种、L6/W10、seed=20261005、Development组级检验；显著性使用加入A2后的跨物种探索检验族 `q_new_exploratory≤.05`，每类背景522个检验。不是18次A2运行的总数，不是功能Motif数量。

## 1. 技术路线与实际实现

参考材料 {link('docs/参考资料/motif-embedding-methods.pdf')} 第4–5页列出预训练表示的注意力提取、输入重建、Embedding聚类、SAE解释等路线；说明表示模型本身不完成Motif提取。后续章节讨论注意力不等同于因果贡献，以及SAE需要额外训练和特征验证。这些是**路线选项与研究动机**，并没有给本项目指定ProkBERT-mini、PCA=24或K=24。文档中的外部论文案例和性能数字未作为本项目实验事实。

**我们自行确定的方案**是冻结ProkBERT-mini，先实现局部表示→聚类→真实片段→PWM→独立验证的基础闭环，与One-hot共用后端，便于隔离表示差异、控制CPU成本、保留坐标追溯。SAE作为基础候选通过稳定性和富集门槛后的扩展；基础证据不足，因此未启动。注意力解释也未执行。没有注意力/SAE对照，不能宣称这两类方法更差。来源：{link(b+'README.md')}、{link(b+'config.json')}、{link(b+'encoder.py')}、{link(b+'run.py')}。

| 项目 | 实际设置 |
| --- | --- |
| 编码器 | Frozen `neuralbioinfo/prokbert-mini`；eval、关闭梯度，无分类微调 |
| 模型版本 | revision `feb2520a43cd9cdb5b3d8477e47209dbcb55d1dc`，加载证据在 `audit/model.json` |
| 真实输入 | TJU六物种81bp，TSS 1-based=61 / 0-based=60；6745条Discovery正序列、967条Development正序列；上下游方向沿用冻结输入 |
| Token | 6-mer、stride=1；81bp有76个真实Token，排除首尾特殊Token；hidden维度384 |
| 主设置 | hidden layer6、窗口10bp、步长1；每序列72个窗口 |
| 降维与聚类 | PCA24（randomized、给定seed）；实际是MiniBatchKMeans，K=24、n_init=5、batch_size=2048；CPU4线程 |
| PWM | 同长片段，伪计数每位置每碱基0.5；最低20个启动子ID/簇 |
| 种子 | 真实OH/A0六物种各3种子20261005/06/07；A1六物种主种子，额外2种子只在枯草芽孢杆菌；A2六物种各3种子 |

**实际窗口范围：** 首轮六物种主比较、A1诊断、A2真实对照、首轮toy均默认10bp。首轮枯草芽孢杆菌A0主种子额外运行6/8/12/16bp，L3/W10层消融，以及OH/W16。`multiscale`实际在该物种3种子运行：6/8/10/12/16按共同中心对齐后拼接（1920维），再PCA；用16bp外包络取DNA。缓存所有尺度不代表六物种都执行了全部宽度实验。后续合成单模式主宽10bp、双模式6bp；S3随机强模式另测6/16bp。来源：{link(b+'summary/run_summary.tsv')}、{link(n+'benchmark.py')}。

设目标窗口W=[s,s+w)，Token区间I_t=[t,t+6)，向量h_t，覆盖数c_b，重合长度o_t=|W∩I_t|：

- **A0**：先碱基平均 g_b=Σ[t:b∈I_t]h_t/c_b，再窗口平均 z_W=Σ[b∈W]g_b/w。Token权重 α_t=Σ[b∈W∩I_t]1/(w·c_b)。
- **A1**：仅取I_t完整包含于W内的Token，z_W=Σ[t:s≤t,t+6≤s+w]h_t/(w−5)。10bp窗口平均5个Token，6bp窗口仅1个。
- **A2**：α_t=(o_t/6)/Σ_u(o_u/6)，z_W=Σ_tα_th_t。实现为归一化权重矩阵与Token数组einsum；不是重新编码局部片段。

所有h_t来自整条81bp输入的上下文化表示；A1/A2不能宣称消除了窗口外上下文。公式与实现见 {link(n+'common.py')}、{link(b+'encoder.py')}。

**聚类到候选的真实步骤：** 对全部Discovery正类窗口拟合PCA及MiniBatchKMeans → 计算每窗口在PCA空间到所属簇中心的欧氏距离 → 每簇按距离升序稳定排序 → 同一启动子ID只保留最近的一个窗口 → 少于20个启动子则弃簇 → 按原始DNA坐标截取同长正向片段，无重新对齐/变长边界 → 计算P(i,b)=(N(i,b)+0.5)/(N+2) → 按簇编号顺序与先保留候选比较，全宽、允许反向互补的中心化cosine≥.95则映射到代表，否则保留 → 保存 `motifs.meme`、Logo、实例、PFM及merge_map。

这里“合并”只是**丢弃相似矩阵、保留较早代表及映射**，没有合并两簇片段重新估计PWM。`instances.tsv`仍含被丢弃候选的实例，首轮 `candidates.tsv`也为去重前候选；最终数必须数 `motifs.meme` 或manifest的motifs。一个簇不保证对应一个最终PWM。同一启动子可进入不同簇；拟合片段按启动子ID去重，未按leakage_group降权，也未按相同DNA字符串去重。

## 2. 六物种主要结果

每格为**保留候选数 / 天然显著数 / Shuffle显著数**。统一主种子20261005；OH/A0/A1取首轮冻结结果，A2取后续轮次，STREME取20261008冻结天然主候选，并在本分支以相同自定义扫描口径重评价。所有显著数使用同一个新增探索族，天然与Shuffle分别BH，每族522个。

{mdtable(['物种','One-hot','ProkBERT A0','A1','A2','STREME'],rows)}

来源：{link(n+'real_data_comparison/primary_seed_counts.tsv')}、{link(n+'real_data_comparison/all_methods_enrichment.tsv')}。

旧主种子族只有406个检验（未加入A2），天然显著数为OH27、A0 39、A1 35、STREME14，Shuffle均0。新增族变成OH28、A0 40、A1 36、A2 45、STREME14，**差异来自校正族变化，不是重跑导致性能提升**。旧证据：{link(b+'summary/baseline_significant_counts.tsv')}。

A2三个种子分别保留116/109/113个、天然显著45/32/37个、Shuffle均0；后两个种子的新增族分别402/409个，因为A1其他五物种没有重复种子。因此不能把116说成全部18次运行的候选总数，也不能把跨种子相加当不同新Motif。

STREME原发现为6物种×天然/二核苷酸背景12个运行，搜索宽5–15、order=2、nmotifs=10；只从天然发现结果按内部holdout p≤.05排序，每物种最多5个，精确PWM或精确RC重复去重，合计17个。**STREME内部holdout是Discovery内部10%划分，不是项目Holdout**。其候选规则与Embedding的24簇/≥.95相似去重不同；17个不是双背景所有发现的总数。依据：{link('results/motif/tju_streme_v2_20261008/public_evidence/run_manifest.json')}。

数量不能直接排名：预算、长度、候选去重与发现目标不同；天然背景有组成差异；一个候选可高度依赖其他候选；q是背景关联而非功能证据。本分支对STREME的Development计数也不是原STREME内部p、严格FIMO位点q或其他主线阶段的计数。

## 3. 三个真实候选的完整例子

选例全部来自**枯草芽孢杆菌、主种子20261005、L6/W10**，用于对比证据类型，不代表六物种总体。共识为每列最大概率碱基（ACGT顺序处理并列），不是IUPAC，也不要求某条实例完全等于共识。位置为Motif中心相对TSS，负值表示上游；以下分布是描述性分位数，无位置富集p/q。
'''
    labels=['PWM相似而建模位置不一致，GC分层后不支持关联','天然富集明显，但Shuffle不支持','Shuffle原始p较小，BH后仍未获支持']
    for i,(x,label) in enumerate(zip(r['examples'],labels),1):
        ev={a['control']:a for a in x['enrichment']}
        dev={a['control']:a for a in x['development']}
        p=ev['natural'];s=ev['shuffle'];sample=x['first_five_fragments'][0]
        text+=f'''\n### 例{i}：{x['motif']}（{x['method']}）

**证据类型：{label}。** 长度{x['width']}bp，共识 `{x['consensus_argmax']}`；构建片段{x['instance_fragments']}条、启动子ID{x['unique_promoter_ids']}个、关联组{x['unique_leakage_groups']}个；本代表无其他候选合并到它。Logo：{link(x['logo'])}。

真实建模片段示例：`{sample['sequence']}`，父ID `{sample['sequence_id']}`，0-based区间[{sample['start0']},{sample['end0_exclusive']})，中心相对TSS={sample['center_tss']}bp。完整片段与PWM：{link(x['folder']+'/instances.tsv')}、{link(x['folder']+'/motifs.meme')}。

| Development对照 | 命中序列/扫描序列 | 命中组/评价组 |
| --- | --- | --- |
| 正启动子 | {dev['positive']['hit_sequences']}/{dev['positive']['scanned_sequences']} | {int(p['positive_hits'])}/{int(p['positive_groups'])} |
| 天然负例 | {dev['natural']['hit_sequences']}/{dev['natural']['scanned_sequences']} | {int(p['control_hits'])}/{int(p['control_groups'])} |
| 配对Shuffle | {dev['shuffle']['hit_sequences']}/{dev['shuffle']['scanned_sequences']} | {int(s['control_hits'])}/{int(s['control_groups'])} |

命中阈值={number(p['threshold'])} bits，严格大于阈值；Development正类命中方向：{dev['positive']['hit_strands']}。

| 检验 | 效应量 | 原始p | 新探索族q（522） |
| --- | --- | --- | --- |
| 天然单侧Fisher | OR={number(p['odds_ratio'])}；命中率 {int(p['positive_hits'])}/{int(p['positive_groups'])} 对 {int(p['control_hits'])}/{int(p['control_groups'])} | {number(p['p'])} | {number(p['q_new_exploratory'])} |
| Shuffle配对单侧二项 | 正独有b={int(s['discordant_pos_only'])}，null独有c={int(s['discordant_null_only'])}；率差={(s['discordant_pos_only']-s['discordant_null_only'])/s['positive_groups']:.6f} | {number(s['p'])} | {number(s['q_new_exploratory'])} |
'''
        if x['method']=='A0':
            text+=f"\n本候选旧406族天然q={number(p['q_original'])}、Shuffle q={number(s['q_original'])}；单运行q分别{number(p['q'])}/{number(s['q'])}。\n"
        else:
            text+=f"\n本候选没有旧406族q；A2单运行19候选族q分别为{number(p['q_run'])}/{number(s['q_run'])}，不能冒充跨六物种q。\n"
        positions=[]
        for name,dist in [('Discovery建模片段',x['construction_center_tss']),('Development正类命中',dev['positive']['hits_center_tss'])]:
            positions.append([name,dist['n'],f"{dist['median']:.2f}",f"[{dist['q25']:.2f}, {dist['q75']:.2f}]",f"[{dist['min']:.2f}, {dist['max']:.2f}]"])
        text+='\n'+mdtable(['位置分布','n','中位数bp','Q25–Q75 bp','范围bp'],positions)+'\n'
        if x['method']=='A0':
            stab=pd.DataFrame(x['seed_stability'])
            stabilityrows=[]
            for seed in [20261006,20261007]:
                d=stab[stab.seed_b==seed]
                a=d.iloc[0]
                vals=[d[(d.stage==stage)&(d.tolerance_bp==t)].physical_position_jaccard.iloc[0] for stage in ['discovery_instances','development_scans'] for t in [0,2]]
                stabilityrows.append([seed,a.motif_b,f'{a.pwm_similarity:.4f}']+[f'{v:.4f}' for v in vals])
            text+='\n'+mdtable(['匹配种子','匹配PWM','全宽cosine','建模J精确','建模J±2','扫描J精确','扫描J±2'],stabilityrows)+'\n'
            cmh=x['gc_adjustment'][0];best=x['streme_alignment'][0]
            text+=f"\nGC分层CMH探索：OR={number(cmh['cmh_or_exploratory'])}，p={number(cmh['cmh_p_exploratory'])}，独立CMH族q={number(cmh['cmh_q_exploratory'])}，正组覆盖{cmh['cmh_positive_coverage']:.4f}、天然组覆盖{cmh['cmh_natural_coverage']:.4f}。这提示背景组成需考虑，不证明全部关联由GC造成。\n"
            text+=f"\n与STREME已实际做描述性矩阵比对：最高对应 `{best['streme']}`，cosine={best['cosine_similarity']:.4f}、offset={best['offset']}、方向{best['strand']}；该对齐仅5列重合（默认最小重合5），无TOMTOM显著性，不能叫同一个功能Motif。{link(x['folder']+'/streme_pwm_similarity.tsv')}。该Embedding候选的已知位点、结合蛋白或功能检验：当前没有证据。\n\n**最终定位：低复杂度A-rich PWM候选；有天然背景关联及矩阵重复性，未获得GC分层和Shuffle的独立支持，不认定新功能模式。**\n"
        else:
            sr=[[a['run_b'].split('__s')[-1],a['motif_b'],f"{a['full_pwm_similarity']:.4f}",f"{a['instance_exact_jaccard']:.4f}"] for a in x['seed_stability']]
            text+='\n'+mdtable(['匹配种子','匹配PWM','全宽cosine','建模精确Jaccard'],sr)+'\n'
            text+='\n本A2候选的扫描位点跨种子容差稳定性、候选级CMH、STREME矩阵比对：**未计算**；已知位点/功能对应：**当前没有证据**。已有A2总体矩阵/建模实例匹配不能替代这些检验。\n'
            if x['motif']=='A2_C10':text+='\n**最终定位：天然背景富集的AT-rich PWM候选；没有超过二核苷酸背景的显著支持，建模精确定位在两次匹配中均为0，不认定新功能模式。**\n'
            else:text+='\n**最终定位：未校正Shuffle富集趋势的探索候选；两类跨物种q均未通过，不能以p≈.011称为显著发现。**\n'
        text+=f"\n扫描/检验原表：{link(x['folder']+'/development_scans.tsv')}、{link(x['folder']+'/development_enrichment.tsv')}；跨种子证据：{link(n+('stability_diagnostics/tolerance_stability.tsv' if x['method']=='A0' else 'real_data_comparison/A2_seed_stability.tsv'))}。\n"
    text+=f'''
## 4. PWM验证的实际执行顺序与结论边界

1. **发现/校准分离。** Discovery正启动子用于PCA、聚类、实例和PWM；标签已用于筛选启动子，但无位点标签指导聚类。Discovery天然负例用于对称零阶背景b（各碱基总数加0.5，再令A=T、C=G）及阈值；不是最终独立检验。真实Development正/天然/Shuffle用于固定PWM验证，反复查看后属于探索性。真实Holdout在本分支未评价；不能据此宣称最终泛化或对照全部方法的最终性能。
2. **双链最大分扫描。** 对每个PWM、每条序列所有合法起点算 S=Σ_i log₂[P(i,x_i)/b(x_i)]，同时扫描反向互补PWM；每序列仅保留全局最高分位点。并列默认先正链、再较早位置。阈值为该PWM在Discovery天然负序列最大分的95%分位（`method='higher'`），命中严格score>threshold。它是经验序列级校准，不是FIMO单个位点p/q，也不验证同一序列内所有实例。
3. **天然富集。** 每leakage_group的presence取组内最大值；两类共享组从双方排除。2×2命中/未命中表做Fisher单侧greater，保存OR及0.5校正近似区间。支持给定背景下的出现率关联，不支持因果功能；序列级命中数不等于组数。
4. **二核苷酸Shuffle。** 复用冻结FASTA，每正序列一个null，ID追加`__dinucleotide_null`并继承父组；原重建seed20261005。随机Euler trail重排相邻碱基边，保持长度、单碱基和二核苷酸计数；不是所有重排均匀抽样。原构建拒绝原序列不变，最多128尝试，Discovery另排除与测试端关联的重排。组级配对后，正独有b、null独有c，p=Pr[Binomial(b+c,.5)≥b]，b+c=0时p=1。支持超出该二核苷酸背景的关联，仍不是功能证明。来源：{link('preprocessing/shuffle.py')}、{link('preprocessing/rebuild_audited_data.py')}。
5. **GC诊断。** Discovery/Development按物种/标签看GC均值和20bin分布重合；Development组GC取均值、presence取最大，排除混合组后固定[0,.2,.4,.6,.8,1]分层，层内两类各≥5组才纳入CMH。CMH检验的原假设为共同OR=1，使用双侧检验，独立对406个有限p做BH，不替换原q。保存层效应和覆盖比例。它检查粗GC组成的替代解释，未解决二核苷酸、局部背景和样本来源混杂；只覆盖首轮主种子OH/A0/A1/STREME，没有A2候选级CMH。
6. **BH范围。** 原单运行`q`：同物种/方法/seed内，本方法保留候选加固定STREME，按天然/Shuffle分开。A2单运行`q_run`仅该运行A2候选。首轮`q_six_species`/后续`q_original`：每seed六物种主配置OH/A0/A1+一份STREME，背景分开；主种子m=406。新增`q_new_exploratory`追加A2、排除A0_recalc，主种子m=522，其他seed m=402/409；其他层/宽消融不入这个主族。各族均不等于把所有尝试过的实验合成一个全局校正，因此整个探索仍需保留选择偏差限制。
7. **跨种子匹配。** PWM逐列减.25后展平成向量，搜索相对偏移及正/反向互补，cosine取最大；核心同宽10bp比较要求全宽重合，局部诊断另用最低6列。Hungarian一对一匹配最大总相似度；它是描述性匹配，没有矩阵相似p/q。
8. **位点Jaccard。** 首轮建模集合为(sequence_id,start0)，J=交集/并集。后续在共有ID上一对一比较起点，容差t下匹配m，J=m/(n_a+n_b−m)，未共有ID仍留在分母。物理位置版忽略链；方向版按PWM匹配方向翻转B链并按offset平移。独立扫描只取Development正类present位点。不是逐碱基交并，也不是关联组集合；家族种子配对不是独立生物学重复。200次均匀/经验位置参照没有转换为正式p。

实现依据：{link(b+'discovery.py')}、{link(b+'run.py')}、{link(n+'real.py')}、{link(b+'summarize.py')}、{link(n+'collect.py')}、{link(n+'diagnostics.py')}。扫描位置支持描述性位置/间距观察，未计算候选例子的TSS显著性或组合功能。

## 5. 最重要的结果与适用范围

**矩阵相似而实例不稳。** A0主配置六物种3种子共{r['first_round_A0_stability']['pwm_family_pairs']}个全宽家族配对，全部匹配的PWM cosine中位数{r['first_round_A0_stability']['all_family_median_pwm_cosine']:.4f}，建模精确Jaccard中位数{r['first_round_A0_stability']['all_family_median_exact_jaccard']:.6f}。A2对应{r['A2_all_family_stability']['pairs']}对，分别{r['A2_all_family_stability']['median_pwm_cosine']:.4f}/{r['A2_all_family_stability']['median_exact_jaccard']:.6f}，不能凭小数差异宣称改善。

进一步限定A0全宽cosine≥.8且有共同ID后：
'''
    st=pd.DataFrame(r['A0_tolerance_summary'])
    rows=[]
    for stage,name in [('discovery_instances','Discovery建模片段'),('development_scans','Development独立扫描')]:
        d=st[st.stage==stage]
        rows.append([name,int(d.iloc[0].paired_families)]+[f'{d[d.tolerance_bp==t].position_only_median.iloc[0]:.4f}' for t in [0,1,2,3]])
    text+=mdtable(['对象','家族配对数','J精确','J±1bp','J±2bp','J±3bp'],rows)
    text+=f'''

上述忽略方向；方向版精确/±2bp分别为建模.0011/.0033、扫描.1579/.3600。两阶段分母和抽取逻辑不同，不能解释成同一组位点得到修复。图：{link(n+'figures/F02_stability.png')}；表：{link(n+'stability_diagnostics/A0_stability_summary.tsv')}。

**A2等价性。** 窗口6/10/16bp的等价位置分别66/76、62/72、56/66。10bp时s=5…66等价，只有两端各5个位置不同；因为内部每碱基被6个Token覆盖，A0权重变成o_t/(6w)，与归一化A2相同。这是聚合权重事实，不是“模型完全等价”或“定位问题已解决”。边缘变化可通过全局PCA/聚类影响内部候选。六物种主种子的A0同einsum重实现对照候选数均未变、匹配PWM中位cosine均1.0，属于实现控制。依据：{link(n+'pooling_ablation/weight_equivalence.tsv')}、{link(n+'real_data_comparison/A0_numerical_control.tsv')}、{link(n+'figures/F06_pooling.png')}。

**合成结果。** 九类条件×3种子；每条件Discovery正256、Development正128、独立合成Test正128及相应空白对照，另256条校准对照；共132个表示/窗口拟合。下表是严格参考匹配后的Test±2bp定位F1，3种子均值，双模式再平均2个家族。
'''
    synthetic=pd.DataFrame(r['synthetic_test_summary'])
    conditions=[('S1_strong_fixed','强模式固定位置'),('S2_mut15_fixed','固定位置15%替换'),('S2_mut30_fixed','固定位置30%替换'),('S3_strong_random','强模式随机位置'),('S4_gc25_random','随机位置，背景GC=.25'),('S4_gc75_random','随机位置，背景GC=.75'),('S5_pair_fixed','双模式gap=17'),('S6_pair_variable','双模式gap=10–24'),('S7_markov_random','随机位置，Markov持续=.7')]
    rows=[]
    for sid,label in conditions:
        rows.append([label,6 if sid.startswith(('S5','S6')) else 10]+[f'{synthetic[(synthetic.scenario==sid)&(synthetic.method==m)].iloc[0]["mean"]:.3f}' for m in ['OH','A0','A1','A2']])
    text+=mdtable(['条件','窗口bp','OH','A0','A1','A2'],rows)
    text+=f'''

指标限定：候选与预设真值PWM一对一匹配，cosine≥.8且至少覆盖较短矩阵全宽；扫描方向正确、起点投射误差≤2bp记TP。precision分母包含正序列全部阳性预测和空白对照误报，位置错误计FP且漏真值，recall分母为全部植入真值；每家族每序列仅一个最大分预测。**未过PWM匹配门槛也记0，不等于hidden state无信息。** 真值只用于家族映射和评价、不进拟合；这是参考辅助恢复，非真实功能验证。双模式零值范围可由某家族失败产生，需与逐家族表一起读。

宽度消融（S3随机强模式）：6bp OH/A0/A1/A2=.981/.649/.963/.977；10bp=.981/.970/.982/.992；16bp=.955/0/.624/0。6bp匹配12bp真值只证明部分锚点定位，不能称完整边界恢复。高GC和Markov部分Embedding种子失败，完整min/max在原表，均值非置信区间。

双模式整体gap±3bp恢复率（分母每方法384条Test，包含漏检）：固定gap OH/A0/A1/A2=.6458/0/.6458/.2969；可变gap=.6536/.6536/.6536/.6276。该指标与上表的家族定位F1不同；只画双检出者会漏掉失败，不能据间距图证明真实启动子grammar。

首轮强模式toy的合格候选数：OH3、A0 0；窗口内Token版本OH3、A0 0、A1 1。旧准则是建模片段重叠≥8bp的purity≥.5、PWM全宽相似≥.8且Development q≤.05，**与后续Test定位F1不同，不能拼成提升曲线**。后续严格评价版本保留了原宽松6列局部匹配结果，所有方法统一改为较短矩阵全宽再评价冻结PWM；本次直接使用已保存严格结果，没有改规则。

来源：{link(n+'synthetic_benchmark/strict_test_summary.tsv')}、{link(n+'synthetic_benchmark/strict_metrics.tsv')}、{link(n+'synthetic_benchmark/grammar_summary.tsv')}、{link(n+'configs/strict_evaluation.json')}、{link(n+'benchmark.py')}、{link(b+'synthetic/manifest.json')}、{link(b+'synthetic_tokenwindow/manifest.json')}。已有图：{link(n+'figures/F07_synthetic.png')}、{link(n+'figures/F08_width.png')}、{link(n+'figures/F09_pair.png')}。

## 6. 能支持与尚未完成的结论

可以说：在本项目冻结ProkBERT-mini+窗口/PCA/聚类/PWM流程和既有合成条件下，One-hot在退化模式及部分宽度/背景条件更可靠；Embedding能产生真实可追溯候选，强模式部分条件表现好，但没有可靠整体额外发现价值。不能说：One-hot在所有启动子/所有预训练表示上都更优，或预训练模型完全无用。本轮合成没有STREME，真实数据没有完整位点真值，不具备总体真实召回率排名。

**本分支当前没有获得独立支持的新功能Motif。** 已完成的是候选发现、PWM恢复、背景关联评价、定位/聚合诊断和受控合成恢复能力分析。天然q通过不等于功能确认，STREME局部相似不等于已知功能对应，候选数量不等于新调控元件数量。

尚未完成：真实Holdout/独立来源的最终候选验证；Embedding候选可靠已知位点/数据库显著匹配；完整真实位点召回率；功能扰动或结合蛋白确认；显著位置/间距/组合机制；A2候选级CMH与独立扫描跨种子容差诊断。A1六物种3种子未齐，A3拟合、A4局部重新编码、SAE、模型微调均未执行。本次不安排补做。

完整研究依据保留在 {link(b+'REPORT.md')} 和 {link(n+'REPORT.md')}。本次机器可读提取为 {link(OUT/'事实数据.json')}，只汇总原表；候选位置分位数及共识属于允许的简单描述统计，并非新显著性检验。
'''
    (OUT/'实验事实汇总.md').write_text(text,encoding='utf-8')

def main():
    OUT.mkdir(exist_ok=True)
    allresults=pd.read_csv(NEXT/'real_data_comparison/all_methods_enrichment.tsv',sep='\t')
    newcounts=pd.read_csv(NEXT/'real_data_comparison/primary_seed_counts.tsv',sep='\t')
    oldcounts=pd.read_csv(BRANCH/'summary/baseline_significant_counts.tsv',sep='\t')
    seedstable=pd.read_csv(NEXT/'stability_diagnostics/tolerance_stability.tsv',sep='\t')
    a2stable=pd.read_csv(NEXT/'real_data_comparison/A2_seed_stability.tsv',sep='\t')
    examples=[]
    chosen=[('A0','EMB_prokbert_C22'),('A2','A2_C10'),('A2','A2_C18')]
    source='tjupan_bacillus_subtilis'
    for method,motif in chosen:
        folder=(BRANCH/'runs'/f'{source}__prokbert__L6W10__s20261005') if method=='A0' else NEXT/'real_data_comparison'/f'{source}__A2__s20261005'
        pwm=read_meme(folder/'motifs.meme')[motif]
        inst=pd.read_csv(folder/'instances.tsv',sep='\t');inst=inst[inst.motif==motif]
        scans=pd.read_csv(folder/'development_scans.tsv',sep='\t');scans=scans[scans.motif==motif]
        ev=allresults[(allresults.source==source)&(allresults.method==method)&(allresults.motif==motif)&(allresults.seed==20261005)]
        record={'method':method,'motif':motif,'source':source,'seed':20261005,'folder':str(folder.relative_to(ROOT)),
          'width':len(pwm),'consensus_argmax':''.join('ACGT'[i] for i in pwm.argmax(axis=1)),
          'logo':str((folder/f'{motif}.png').relative_to(ROOT)),
          'instance_fragments':len(inst),'unique_promoter_ids':inst.sequence_id.nunique(),
          'unique_leakage_groups':inst.leakage_group.nunique(),'construction_center_tss':distribution(inst.center_tss),
          'first_five_fragments':inst[['sequence_id','start0','end0_exclusive','sequence','center_tss']].head(5).to_dict('records'),
          'development':[],'enrichment':ev.to_dict('records'),
          'merge_map_for_candidate':[r for r in json.loads((folder/'merge_map.json').read_text()) if r['motif']==motif or r['representative']==motif]}
        assert (folder/f'{motif}.png').is_file()
        for control,g in scans.groupby('control'):
            h=g[g.present]
            record['development'].append({'control':control,'scanned_sequences':len(g),'hit_sequences':len(h),
              'hits_center_tss':distribution(h.center_tss),'hit_strands':h.strand.value_counts().to_dict()})
        if method=='A0':
            s=seedstable[(seedstable.source==source)&(seedstable.method=='prokbert')&(seedstable.strategy=='full')&
                         (seedstable.seed_a==20261005)&(seedstable.motif_a==motif)&seedstable.tolerance_bp.isin([0,2])]
            record['seed_stability']=s[['stage','seed_a','seed_b','motif_a','motif_b','pwm_similarity','alignment_orientation',
                 'tolerance_bp','n_a','n_b','common_sequence_ids','physical_position_jaccard','direction_aware_aligned_jaccard']].to_dict('records')
            cmp=pd.read_csv(folder/'streme_pwm_similarity.tsv',sep='\t')
            record['streme_alignment']=cmp[cmp.embedding==motif].sort_values('cosine_similarity',ascending=False).to_dict('records')
            gc=pd.read_csv(NEXT/'gc_diagnostics/motif_gc_adjustment.tsv',sep='\t')
            record['gc_adjustment']=gc[(gc.source==source)&(gc.method=='prokbert')&(gc.motif==motif)].to_dict('records')
        else:
            s=a2stable[(a2stable.source==source)&a2stable.run_a.str.endswith('s20261005')&(a2stable.motif_a==motif)]
            record['seed_stability']=s.to_dict('records')
            record['streme_alignment']='not calculated for A2 in existing artifacts'
            record['gc_adjustment']='not calculated for A2; existing CMH used first-round main seed candidates'
        record['known_function_support']='No independent known-site/function validation for this candidate in branch artifacts'
        examples.append(record)
    firstst=pd.read_csv(BRANCH/'summary/seed_stability.tsv',sep='\t')
    a0=firstst[(firstst.method=='prokbert')&(firstst.layer==6)&(firstst.width==10)]
    algebra=pd.read_csv(NEXT/'pooling_ablation/weight_equivalence.tsv',sep='\t')
    result={'primary_seed':20261005,'main_width':10,'main_embedding_layer':6,
      'new_family_sizes':allresults[allresults.method!='A0_recalc'].groupby(['seed','control']).size().reset_index(name='m').to_dict('records'),
      'new_counts':newcounts.to_dict('records'),'original_counts':oldcounts.to_dict('records'),
      'examples':examples,'first_round_A0_stability':{'pwm_family_pairs':len(a0),
        'all_family_median_pwm_cosine':a0.full_width_cosine.median(),
        'all_family_median_exact_jaccard':a0.exact_instance_jaccard.median()},
      'A0_tolerance_summary':pd.read_csv(NEXT/'stability_diagnostics/A0_stability_summary.tsv',sep='\t').to_dict('records'),
      'A2_all_family_stability':{'pairs':len(a2stable),'median_pwm_cosine':a2stable.full_pwm_similarity.median(),
        'median_exact_jaccard':a2stable.instance_exact_jaccard.median()},
      'pooling_equivalence':algebra.groupby('width').agg(windows=('equivalent','size'),equivalent=('equivalent','sum')).reset_index().to_dict('records'),
      'synthetic_test_summary':pd.read_csv(NEXT/'synthetic_benchmark/strict_test_summary.tsv',sep='\t').to_dict('records'),
      'synthetic_grammar_summary':pd.read_csv(NEXT/'synthetic_benchmark/grammar_summary.tsv',sep='\t').to_dict('records')}
    (OUT/'事实数据.json').write_text(json.dumps(clean(result),ensure_ascii=False,indent=2),encoding='utf-8')
    report(clean(result))
    for r in examples:
        print(json.dumps(clean({k:v for k,v in r.items() if k not in ['first_five_fragments','enrichment']}),ensure_ascii=False))
    print('First-round A0 stability:',result['first_round_A0_stability'])
    print('A2 all-family stability:',result['A2_all_family_stability'])

if __name__=='__main__':main()
