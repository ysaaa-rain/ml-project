"""Data-seed and initialization-stratified report of the frozen-method replication."""
import numpy as np
import pandas as pd

from .pipeline import ART,HERE,PROTOCOL,SUMMARY,read


def table(frame):
    frame=frame.copy()
    for name in frame:
        if frame[name].dtype.kind=='f':frame[name]=frame[name].map(lambda x:f'{x:.3f}')
    return '\n'.join(['| '+' | '.join(map(str,frame.columns))+' |','| '+' | '.join(['---']*len(frame.columns))+' |',
        *['| '+' | '.join(map(str,row))+' |' for row in frame.itertuples(index=False,name=None)]])


def figure(paired,initialization,cfg):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=paired[(paired.scenario=='mut30_random')&(paired.view=='MATCH_SHIFT')].sort_values('data_seed')
    effects=initialization[initialization.scenario=='mut30_random'].sort_values('algorithm_seed')
    fig,axes=plt.subplots(1,2,figsize=(13,5.3))
    colors=['#315f9b','#d06329','#49976c','#8d62a1','#b88822']
    for row,color in zip(rows.itertuples(),colors):
        axes[0].plot([0,1],[row.f1_baseline,row.f1_shift],color=color,alpha=.8,marker='o',markersize=5,label=str(row.data_seed))
    axes[0].plot([-.08,.08],[rows.f1_baseline.mean()]*2,color='black',linewidth=2.5)
    axes[0].plot([.92,1.08],[rows.f1_shift.mean()]*2,color='black',linewidth=2.5)
    axes[0].set_xticks([0,1],['RAW','SHIFT']);axes[0].set_ylabel('Reference-assisted localization F1 (+/-2bp)')
    axes[0].set_ylim(-.03,1.04);axes[0].grid(axis='y',alpha=.2);axes[0].legend(title='Data seed',loc='upper left',fontsize=8,title_fontsize=8)
    x=np.arange(len(effects));bars=['#8d949d' if seed in cfg['reused_initialization_seeds'] else '#247fa1' for seed in effects.algorithm_seed]
    axes[1].bar(x,effects.mean_delta,color=bars,width=.72);axes[1].axhline(0,color='black',linewidth=.8)
    for position,value in zip(x,effects.mean_delta):
        axes[1].annotate(f'{value:+.3f}',(position,value),xytext=(0,3 if value>=0 else -11),
                         textcoords='offset points',ha='center',fontsize=8)
    axes[1].set_xticks(x,[str(n) for n in effects.algorithm_seed],rotation=35,ha='right')
    axes[1].set_ylabel('SHIFT - RAW F1, averaged across 5 data seeds');axes[1].grid(axis='y',alpha=.2)
    from matplotlib.patches import Patch
    axes[1].legend(handles=[Patch(color='#8d949d',label='Previously used'),Patch(color='#247fa1',label='New')],loc='best',fontsize=8)
    fig.suptitle('30% mutation, matched budgets: data-seed pairs and algorithm-seed effects',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.93));fig.savefig(SUMMARY/'initialization_replication.png',dpi=170);plt.close(fig)


def main():
    cfg=read(PROTOCOL);manifest=read(SUMMARY/'manifest.json');verified=read(SUMMARY/'verification.json')
    qa=read(SUMMARY/'independent_qa.json');cross=read(SUMMARY/'cross_version_qa.json');decision=read(SUMMARY/'decision.json')
    assert manifest['status']=='COMPLETED' and verified['status']==qa['status']==cross['status']=='PASSED'
    summary=pd.read_csv(SUMMARY/'test_summary.tsv',sep='\t');paired=pd.read_csv(SUMMARY/'paired_deltas.tsv',sep='\t')
    init=pd.read_csv(SUMMARY/'initialization_effects.tsv',sep='\t');filters=pd.read_csv(SUMMARY/'filter_decisions.tsv',sep='\t')
    metrics=pd.read_csv(SUMMARY/'metrics.tsv',sep='\t');classes=pd.read_csv(SUMMARY/'classification.tsv',sep='\t')
    figure(paired,init,cfg)
    chosen=paired[(paired.scenario=='mut30_random')&(paired.view=='MATCH_SHIFT')]
    weak_oh=summary[(summary.scenario=='mut30_random')&(summary.view=='MATCH_OH')].set_index('representation')
    weak_oh_pairs=paired[(paired.scenario=='mut30_random')&(paired.view=='MATCH_OH')]
    weak_largest=weak_oh_pairs.loc[weak_oh_pairs.delta.abs().idxmax()]
    weak_remaining=weak_oh_pairs[weak_oh_pairs.data_seed!=weak_largest.data_seed]
    mixed_oh=paired[(paired.scenario=='mixed_strand_mut15')&(paired.view=='MATCH_OH')]
    largest=mixed_oh.loc[mixed_oh.delta.abs().idxmax()];remaining=mixed_oh[mixed_oh.data_seed!=largest.data_seed]
    caps=summary[summary.view=='CAP5'].pivot(index='scenario',columns='representation',values='f1_mean').reset_index()
    main_summary=summary[(summary.scenario=='mut30_random')&(summary.view.isin(['FULL','CAP5','MATCH_SHIFT','MATCH_OH']))]
    checks=pd.DataFrame([{'预设门槛':k,'结果':'通过' if v else '未通过'} for k,v in decision['gate_checks'].items()])
    classes=classes[classes.split=='test'].groupby(['scenario','data_seed']).agg(auroc=('auroc','mean'),ap=('average_precision','mean')).reset_index()
    classes=classes.groupby('scenario').agg(auroc=('auroc','mean'),average_precision=('ap','mean')).reset_index()
    weak_init=init[init.scenario=='mut30_random'].copy()
    weak_init['cohort']=weak_init.algorithm_seed.map(lambda s:'reused' if s in cfg['reused_initialization_seeds'] else 'new')
    individual=metrics[(metrics.scenario=='mut30_random')&(metrics.view=='MATCH_SHIFT')&(metrics.split=='test')&(metrics.tolerance_bp==2)]
    eligible=individual.pivot(index=['data_seed','algorithm_seed'],columns='representation',values='reference_eligible')
    values=individual.pivot(index=['data_seed','algorithm_seed'],columns='representation',values='f1')
    gain=((~eligible.TIED_RAW)&eligible.TIED_SHIFT);loss=(eligible.TIED_RAW&(~eligible.TIED_SHIFT))
    unchanged=np.isclose(values.TIED_RAW,values.TIED_SHIFT,atol=1e-12,rtol=0)
    edges=filters.groupby('scenario').agg(activation_selected=('selected','sum'),edge_excluded=('excluded_edges','sum'),retained_filters=('retained','sum')).reset_index()
    edges['excluded_fraction']=edges.edge_excluded/edges.activation_selected.clip(lower=1)
    shift_counts=filters[filters.retained].groupby(['scenario','shift']).size().reset_index(name='retained_filters')
    membership=[]
    for folder in sorted(ART.iterdir()):
        data=read(folder/'data.json')
        for parent in sorted(folder.glob('a*')):
            members={r['representation']:set(r['candidates']) for r in read(parent/'views.json') if r['view']=='MATCH_SHIFT'}
            membership.append({'scenario':data['scenario']['id'],'same_set':members['TIED_RAW']==members['TIED_SHIFT'],
                'overlap':len(members['TIED_RAW']&members['TIED_SHIFT'])})
    weak_members=pd.DataFrame(membership);weak_members=weak_members[weak_members.scenario=='mut30_random']
    status='通过预设进一步合成研究门槛' if decision['advance_gate_passed'] else '未通过预设进一步合成研究门槛'
    guidance=('优先补充不同人工模式及背景条件，继续用独立数据评估定位、误报和初始化稳定性；同时保留OH48对照。'
              if decision['advance_gate_passed'] else
              '将边界处理保留为候选诊断，暂不把它作为稳定提升；若要厘清随机性，先用新协议和新Test独立改变模型权重初始化与训练批次顺序，再扩展人工模式及背景。')
    text=f'''# CNN边界处理：多算法种子独立复核

日期：2026-10-27。协议：{cfg['version']}。本轮为新合成数据和扩展初始化的复核，完整阴性结果保留。

## 主要结果与下一步

30%退化、匹配候选预算时，RAW定位F1 {chosen.f1_baseline.mean():.3f}，SHIFT {chosen.f1_shift.mean():.3f}，配对差{decision['mean_delta']:+.3f}；5组新数据中SHIFT胜出{decision['wins']}组、落后{decision['losses']}组，4个全新初始化中有{decision['new_initializations_positive']}个跨5组数据的平均差为正。候选集空白误报均值为{chosen.family_fpr_baseline.mean():.3f}/{chosen.family_fpr_shift.mean():.3f}。**{status}**。

{guidance} 门槛只决定是否值得进一步合成研究；CNN使用序列级监督，OH48没有同样的分类训练，比较不能直接归因为embedding本身。真实Holdout未用，不能据此宣称真实功能motif或总体模型优势。

匹配5候选条件下，30%退化的SHIFT为{weak_oh.loc['TIED_SHIFT','f1_mean']:.3f}，OH48为{weak_oh.loc['OH48_BG','f1_mean']:.3f}；SHIFT在5组数据中胜出{int((weak_oh_pairs.delta>0).sum())}组、落后{int((weak_oh_pairs.delta<0).sum())}组。混合反向链对OH48的平均差为{mixed_oh.delta.mean():+.3f}，其中最大绝对差来自数据种子{int(largest.data_seed)}（{largest.delta:+.3f}），其余4组差范围为{remaining.delta.min():+.3f}到{remaining.delta.max():+.3f}；均值是否受单组影响须结合逐数据结果。

30%退化对OH48的平均差为{weak_oh_pairs.delta.mean():+.3f}，其中数据种子{int(weak_largest.data_seed)}的差为{weak_largest.delta:+.3f}，其余4组范围{weak_remaining.delta.min():+.3f}至{weak_remaining.delta.max():+.3f}；这是带序列级监督的CNN与无同等分类训练的OH48在5组人工数据上的流程对照，不构成一般性或统计显著优势。

## 固定方法与数据

复用[上一轮已冻结边界方法](../cnn_boundary_20261027/REPORT.md)：双向共享CNN按Discovery序列级正负标签训练，两个PWM方案共用同一模型和同一批边界有效的高激活实例。RAW保持原窗口，SHIFT对每个filter只凭Discovery相对背景信息量在−2到+2bp中选一个统一平移。序列级标签可见；植入位置、参考PWM、Development/Test均不进入拟合、平移或排名。边缘±2bp不完整的实例在两方案均排除，至少20条；参照OH48不使用同样的监督与方向共享。

3条件（强模式、30%退化、混合反向链15%退化）×5个新数据种子202610278–202610282×6初始化20261112–20261117；其中20261112/20261113是上一轮使用过的初始化，20261114–20261117为新增。每条件Discovery植入/空白各256，Development/Test各128/128，另校准空白256。随机种子之间共享同一数据；同一数据种子的不同难度共享空白随机流，不把初始化或场景视作额外独立数据。人工连续12bp模式与真实启动子元件语法不同，窗口宽度预先已知。

沿用上一轮的训练实现：CNN算法种子同时设定初始权重和训练批次顺序；本轮比较的是两项共同改变的稳定性，不能将差异只归因于初始权重。OH48使用同编号作为其聚类随机种子；这两个方法并未共享随机算法或监督方式。

实际{manifest['cnn_trains']}次CNN训练、{manifest['cluster_fits']}次OH48聚类、{manifest['pwm_pools']}组PWM池、{manifest['evaluation_sets']}组预算评价。RAW/SHIFT复用模型和同一实例，不算额外训练。FULL、CAP5和两组匹配预算在Discovery排序后固定，零候选保留0分。全部模型、PWM、平移和候选集合先锁定，再一次性评分Development/Test；不按新Test选择种子、方法、阈值或子集。

## 主要终点逐数据种子

{table(chosen[['data_seed','f1_baseline','f1_shift','delta','family_fpr_baseline','family_fpr_shift','budget_baseline','budget_shift','zero_budget_count']])}

![Initialization replication](summary/initialization_replication.png)

左图每条线为同一数据种子在6个算法种子上的平均；右图每根柱为同一算法种子在5组数据的平均差。每个算法种子同时控制权重初始化和训练批次顺序，所有算法种子共享数据，柱数不是独立数据样本数。

{table(checks)}

符号翻转参照p={decision['sign_flip_p']:.4f}，按5个数据均值枚举32种符号分配；依赖零假设下符号可交换。5组数据双侧p最小0.0625，不宣称p<.05。

## 逐算法种子拆解

{table(weak_init[['algorithm_seed','cohort','mean_delta','data_seeds','data_wins','data_losses']])}

当前30%退化的{len(values)}次算法运行中，{int(gain.sum())}次由RAW未通过全宽参考门槛转为SHIFT通过，{int(loss.sum())}次由通过转为未通过；{int(unchanged.sum())}次RAW/SHIFT F1相同。这是对已固定评分的事后描述；5组数据×6初始化并非30组独立数据。

新增算法种子中，平均差为正的种子其正负数据组数见上表；较少数数据上的变化可以驱动单种子平均差，不应只凭其符号认定跨数据稳定。

RAW与SHIFT分别按Discovery重排序。30%退化的{len(weak_members)}次算法运行中有{int(weak_members.same_set.sum())}次匹配预算的候选集合相同，平均交集{weak_members.overlap.mean():.2f}个。身份相同不保证评分时选中同一个filter；完整结果也须看FULL和CAP5。

## 条件、预算、误报和分类

{table(caps)}

{table(main_summary.drop(columns='scenario'))}

严格全12列参考cosine≥0.8才计算参考辅助F1，命中方向正确、起点误差≤2bp；参考只从冻结候选中选择家族，不用于生成/调整候选。各PWM阈值来自独立校准空白最大分95百分位。候选集误报是Test空白中任一固定候选超阈值的比例，不是FIMO q或正式FDR。

{table(classes)}

分类AUROC/AP只描述序列级分类，不等于定位成功或motif功能证据。

## 平移选择与边界排除

{table(edges)}

{table(shift_counts)}

剔除数按filter实例计，同一序列可以被多个filter计入；RAW/SHIFT共享每个filter的剩余实例。偏移由Discovery信息量决定，频次本身不证明找对真实边界。

## 核验与限制

18项预拟合相关测试通过；{verified['instances_checked']}条实例和{verified['pwms_rebuilt']}个PWM由输入DNA/坐标重建，复算{verified['cnn_activation_maps_recomputed']}组CNN激活、{verified['filter_shift_decisions_recomputed']}项平移、{verified['rankings_recomputed']}组Discovery排序及{verified['evaluation_sets']}组预算集合；{verified['historical_files_unchanged']}个历史保护文件哈希不变。另核验{qa['initialized_models_recomputed']}组模型初始化、训练标签与权重变化，并重建{qa['OH48_pools_reconstructed']}组OH48背景筛选。

每个条件内各分区的exact/RC重复为0，并与前3轮合成实验比较，跨版本exact/RC重合均为0（见 `cross_version_qa.json`）；没有全面近同源审计，不能称绝对生物学独立。{verified['truth_sites_checked']}个植入真值仅用于评分和核验。没有独立重训、真实数据拟合、SAE、FIMO或功能验证。

Python {manifest['python']}，Torch {manifest['dependencies']['torch']}，CPU4线程；整批拟合含事前历史哈希核验约{manifest['fit_seconds']:.1f}秒，不当作单模型独立计时。真实Holdout未使用。

## 复现

```powershell
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.pipeline fit
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.pipeline score
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.pipeline qa
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.qa_cross_version
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.report --figure-only
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.pipeline verify
python -m experiments.pr01_02_embedding_v2.cnn_seed_validation_20261027.report
```

拟合拒绝覆盖已有artifacts。人工序列、模型、实例、聚类对象保留本地，公开协议、代码、聚合/分数据与分初始化统计、边界记录、图表和哈希。原始字节哈希只对应原环境，重训请另开输出目录并说明环境差异。
'''
    (HERE/'REPORT.md').write_text(text,encoding='utf8');print('CNN initialization replication report and figure generated')


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--figure-only',action='store_true')
    if parser.parse_args().figure_only:
        figure(pd.read_csv(SUMMARY/'paired_deltas.tsv',sep='\t'),pd.read_csv(SUMMARY/'initialization_effects.tsv',sep='\t'),read(PROTOCOL))
        print('Replication figure generated')
    else:main()
