"""Data-seed-level report and figure for the preregistered CNN boundary trial."""
import numpy as np
import pandas as pd

from .pipeline import ART,HERE,PROTOCOL,SUMMARY,read


def table(frame):
    frame=frame.copy()
    for name in frame:
        if frame[name].dtype.kind=='f':frame[name]=frame[name].map(lambda x:f'{x:.3f}')
    return '\n'.join(['| '+' | '.join(map(str,frame.columns))+' |','| '+' | '.join(['---']*len(frame.columns))+' |',
        *['| '+' | '.join(map(str,row))+' |' for row in frame.itertuples(index=False,name=None)]])


def figure(paired):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows=paired[(paired.scenario=='mut30_random')&(paired.view=='MATCH_SHIFT')].sort_values('data_seed')
    fig,axes=plt.subplots(1,2,figsize=(11.5,5))
    colors=['#315f9b','#d06329','#49976c','#8d62a1','#b88822']
    for axis,columns,ylabel in [(axes[0],['f1_baseline','f1_shift'],'Reference-assisted localization F1 (+/-2bp)'),
                                (axes[1],['family_fpr_baseline','family_fpr_shift'],'Any-candidate blank FPR')]:
        for row,color in zip(rows.itertuples(),colors):
            axis.plot([0,1],[getattr(row,columns[0]),getattr(row,columns[1])],color=color,alpha=.7,
                      marker='o',markersize=5,label=str(row.data_seed))
        axis.plot([-.10,.10],[rows[columns[0]].mean()]*2,color='black',linewidth=2.5)
        axis.plot([.90,1.10],[rows[columns[1]].mean()]*2,color='black',linewidth=2.5)
        axis.set_xticks([0,1],['RAW','SHIFT']);axis.set_ylabel(ylabel);axis.set_ylim(-.03,1.04);axis.grid(axis='y',alpha=.2)
    axes[0].legend(title='Data seed',loc='upper left',fontsize=8,title_fontsize=8)
    fig.suptitle('30% mutation, matched candidate budgets: each pair averages two CNN seeds',fontsize=11)
    fig.tight_layout(rect=(0,0,1,.94));fig.savefig(SUMMARY/'boundary_validation.png',dpi=170);plt.close(fig)


def main():
    cfg=read(PROTOCOL);manifest=read(SUMMARY/'manifest.json');verified=read(SUMMARY/'verification.json')
    qa=read(SUMMARY/'independent_qa.json');decision=read(SUMMARY/'decision.json')
    assert manifest['status']=='COMPLETED' and verified['status']==qa['status']=='PASSED'
    summary=pd.read_csv(SUMMARY/'test_summary.tsv',sep='\t');paired=pd.read_csv(SUMMARY/'paired_deltas.tsv',sep='\t')
    seeds=pd.read_csv(SUMMARY/'per_data_seed.tsv',sep='\t')
    filters=pd.read_csv(SUMMARY/'filter_decisions.tsv',sep='\t');classes=pd.read_csv(SUMMARY/'classification.tsv',sep='\t')
    figure(paired)
    chosen=paired[(paired.scenario=='mut30_random')&(paired.view=='MATCH_SHIFT')]
    main_summary=summary[(summary.scenario=='mut30_random')&(summary.view.isin(['FULL','CAP5','MATCH_SHIFT','MATCH_OH']))]
    weak_oh=summary[(summary.scenario=='mut30_random')&(summary.view=='MATCH_OH')].set_index('representation')
    mixed_oh=paired[(paired.scenario=='mixed_strand_mut15')&(paired.view=='MATCH_OH')]
    largest=mixed_oh.loc[mixed_oh.delta.abs().idxmax()];remaining=mixed_oh[mixed_oh.data_seed!=largest.data_seed]
    direction_rows=[]
    for scenario in cfg['scenarios']:
        cohort=seeds[(seeds.scenario==scenario['id'])&(seeds.view=='CAP5')].pivot(index='data_seed',columns='representation',values='f1')
        delta=cohort.TIED_SHIFT-cohort.FWD_SHIFT
        direction_rows.append({'scenario':scenario['id'],'mean_tied_minus_forward':delta.mean(),
            'tied_wins':int((delta>0).sum()),'tied_losses':int((delta<0).sum()),'ties':int((delta==0).sum())})
    cap=summary[summary.view=='CAP5'].pivot(index='scenario',columns='representation',values='f1_mean').reset_index()
    membership=[]
    for folder in sorted(ART.iterdir()):
        data=read(folder/'data.json')
        for parent in sorted(folder.glob('a*')):
            members={r['representation']:set(r['candidates']) for r in read(parent/'views.json') if r['view']=='MATCH_SHIFT'}
            membership.append({'scenario':data['scenario']['id'],'same_set':members['TIED_RAW']==members['TIED_SHIFT'],
                'overlap':len(members['TIED_RAW']&members['TIED_SHIFT'])})
    membership=pd.DataFrame(membership)
    weak_members=membership[membership.scenario=='mut30_random']
    by_seed=classes[classes.split=='test'].groupby(['scenario','data_seed','direction']).agg(auroc=('auroc','mean'),ap=('average_precision','mean')).reset_index()
    class_summary=by_seed.groupby(['scenario','direction']).agg(auroc=('auroc','mean'),average_precision=('ap','mean')).reset_index()
    kept=filters[filters.retained];shift_counts=kept.groupby(['scenario','direction','shift']).size().reset_index(name='retained_filters')
    edges=filters.groupby(['scenario','direction']).agg(activation_selected=('selected','sum'),edge_excluded=('excluded_edges','sum'),
        retained_filters=('retained','sum')).reset_index()
    edges['excluded_fraction']=edges.edge_excluded/edges.activation_selected.clip(lower=1)
    status='通过预设进一步独立复核门槛' if decision['advance_gate_passed'] else '未通过预设进一步独立复核门槛'
    guidance=('先增加模型初始化种子并在新的独立合成数据上复核平移规则，分析候选误报和方向结构；不要直接推断真实生物学效果。'
              if decision['advance_gate_passed'] else
              '不直接扩大到真实序列。结合逐数据差值、方向消融和偏移分布判断假设是否成立；如继续改发现机制，应另开协议并更换Test。')
    checks=pd.DataFrame([{'预设门槛':k,'结果':'通过' if v else '未通过'} for k,v in decision['gate_checks'].items()])
    text=f'''# CNN参考无关边界与方向消融

日期：2026-10-27。协议：{cfg['version']}。本轮为新的合成数据验证，完整结果与阴性结果均保留。

## 主要结果与指导

30%退化、匹配候选预算时，TIED_RAW平均定位F1为{chosen.f1_baseline.mean():.3f}，TIED_SHIFT为{chosen.f1_shift.mean():.3f}，配对差{decision['mean_delta']:+.3f}；5个新数据种子中SHIFT胜出{decision['wins']}个、落后{decision['losses']}个。候选集空白误报均值分别为{chosen.family_fpr_baseline.mean():.3f}与{chosen.family_fpr_shift.mean():.3f}。**{status}**。

{guidance}

预设门槛仅用于决定是否值得更大规模的独立合成复核，不证明真实启动子功能，也不能把CNN监督方案与OH48的差距直接归因于embedding。主要比较逐数据种子的效果、误报和预算见下表；算法种子不当作独立数据样本。

与OH48同预算比较，30%退化时SHIFT为{weak_oh.loc['TIED_SHIFT','f1_mean']:.3f}，OH48为{weak_oh.loc['OH48_BG','f1_mean']:.3f}，SHIFT在5个数据种子中胜出{int((paired[(paired.scenario=='mut30_random')&(paired.view=='MATCH_OH')].delta>0).sum())}个；两者使用的监督量不同。本次边界处理改善了CNN自身，尚未显示稳定超过OH48。

混合反向链对OH48的5组数据平均差为{mixed_oh.delta.mean():+.3f}，其中数据种子{int(largest.data_seed)}的差为{largest.delta:+.3f}，其余4组差范围为{remaining.delta.min():+.3f}到{remaining.delta.max():+.3f}。该均值受到单组影响明显，不能把它解读为普遍的大幅优势。

事后按10次算法运行拆解主要结果：4次由RAW的全宽参考门槛失败转为SHIFT通过，其余6次RAW已通过且F1保持相同。4次改善均出现在同一个算法种子20261112；另一种子20261113在5组数据中没有F1变化。此项拆解是对已固定结果的探索性解释，10次运行仍只来自5组数据和2个模型初始化，不能当作10次独立重复。收益主要表现为修复部分初始化下的全宽PWM门槛，而不是已通过门槛的候选定位进一步提升；后续应扩大初始化种子复核。

## 方法与公平比较

该方案由[上一轮的窗口偏移诊断](../cnn_features_20261026/REPORT.md)提出，协议在新Test评分前冻结。3种条件（强模式、30%退化、混合反向链15%退化）×5个新数据种子×2算法种子；Discovery植入/空白各256，Development和Test各128/128，另校准空白256。81bp序列和已知12bp连续人工模式不是实测启动子语法。同一种子下不同难度复用空白随机流，条件间不能算独立重复。

实际{manifest['cnn_trains']}次CNN训练、{manifest['cluster_fits']}次OH48聚类、{manifest['pwm_pools']}组PWM池与{manifest['evaluation_sets']}组评价集合。TIED与FWD同一数据/种子使用相同初始权重和批次流；TIED训练/提取共用正反向卷积，FWD只用输入正向，最终PWM评价均双链扫描。两者同时改变训练及提取方向，不能称单独只改变推理方向。

每个filter的Discovery正序列先按其Discovery空白最大激活95百分位固定筛选，每个序列最多一个实例。再排除任何±2bp候选窗口越界的实例：RAW与SHIFT使用**同一批**剩余实例，至少20条，配对避免由边界样本变化造成伪差异。RAW保留原窗口；SHIFT对该filter统一比较−2至+2bp五种窗口，选PWM相对于Discovery空白背景的信息量最大的偏移，平局优先绝对偏移较小、再较小的偏移。反向链移动方向按filter方向定义。该规则不读取参考motif、植入真值或Dev/Test，也没有逐实例自由调整；它可能因Discovery随机波动选择错误方向。

注意本轮RAW也剔除边缘不完整窗口，不能与前一轮未剔边缘的RAW直接当作完全相同流程。额外OH48是完整One-hot48 KMeans+背景拒绝，与CNN的序列级正负监督和反向共享归纳偏好不同。

随机权重和打乱标签控制已在上一轮运行，本轮没有重新训练这两类控制；本轮主要识别同一模型、同一保留实例在窗口处理前后的变化，不能把上一轮控制结果当作本轮独立重复。

所有PWM先写入MEME文本，Discovery排序和后续评分均从相同序列化PWM读取，避免上一轮排名与评分PWM精度不一致。排序是Discovery正减空白命中率及固定平局规则，无参考参与。FULL、CAP5，及按两方法共同候选数裁剪的MATCH_SHIFT/MATCH_OH均保留；m=0保留0分。全体模型、实例、平移选择、排名和预算在任何Dev/Test评分前锁定。

RAW和SHIFT分别对各自PWM做Discovery排序，匹配预算不一定产生相同的候选身份。30%退化条件下，{len(weak_members)}次算法运行中有{int(weak_members.same_set.sum())}次两者的候选集合相同，平均交集{weak_members.overlap.mean():.2f}个；其他条件的差异可能同时包含重排序影响。参考评价仍可能在同一候选集合中选中不同filter。

## 主要终点逐数据种子

{table(chosen[['data_seed','f1_baseline','f1_shift','delta','family_fpr_baseline','family_fpr_shift','budget_baseline','budget_shift','zero_budget_count']])}

![CNN boundary matched validation](summary/boundary_validation.png)

图中每条线对应一个数据种子，两端各自平均两个算法种子；黑线为数据均值。候选集合误报是Test空白中任一固定候选超阈值的比例，不是FIMO q或正式FDR。

{table(checks)}

配对符号翻转参照p={decision['sign_flip_p']:.4f}，32种符号分配；它依赖零假设下符号可交换，5组数据双侧p最小0.0625，不能宣称p<.05。CAP5/FULL、OH48、其他场景及方向对照均作探索性解释。

## 30%退化的全部预算和候选误报

{table(main_summary.drop(columns='scenario'))}

全12列参考cosine≥0.8才允许参考辅助F1，命中方向正确且起点误差≤2bp；参考仅在已锁定子集内选家族，无法找回预算外候选。各PWM阈值来自独立校准空白最大分95百分位；不以Test更改阈值、候选、平移或配置。

## 条件和方向对照：CAP5定位F1

{table(cap)}

RAW和SHIFT是同一个对应方向模型的两个发现输出，不是额外CNN训练。双向与单向比较是包含方向训练和提取的完整流程消融；OH48同样是完整流程对照。数据种子范围和各预算的详细误报、候选数、零预算在 `test_summary.tsv`、`per_data_seed.tsv`。

同一数据种子、CAP5下SHIFT的双向减单向配对比较如下；方向变化并未在所有场景表现出一致优势。

{table(pd.DataFrame(direction_rows))}

## 分类与实例边界

{table(class_summary)}

分类AUROC/AP仅描述模型对序列级标签的预测，不代表定位正确或功能motif发现。不同方向模型共享初始化与训练标签，但其预测不同；不根据Test分类表现选最终模型。

| 场景与方向 | 筛选后实例 | 越界剔除 | 比例 | 保留filter数 |
| --- | ---: | ---: | ---: | ---: |
{chr(10).join(f"| {r.scenario} / {r.direction} | {r.activation_selected} | {r.edge_excluded} | {r.excluded_fraction:.3f} | {r.retained_filters} |" for r in edges.itertuples())}

上表的剔除数按filter实例计，同一序列可被不同filter计入多次，不是独立生物样本。RAW/SHIFT对每个filter使用相同剩余实例；不足20条均不产生PWM。保留filter按Discovery选择的平移如下：

{table(shift_counts)}

平移数值的单位是filter方向的碱基；反链在原输入坐标中的位移符号相反。选择偏移不使用参考；某一偏移被选不等于找到了真实motif边界。

## 实际核验和限制

25项预拟合相关测试通过；{verified['instances_checked']}条实例和{verified['pwms_rebuilt']}个PWM从原始DNA与坐标重建；{verified['cnn_activation_maps_recomputed']}组CNN激活、{verified['filter_shift_decisions_recomputed']}项filter平移决策、{verified['rankings_recomputed']}组Discovery排序和{verified['evaluation_sets']}组冻结预算集合复算通过，{verified['historical_files_unchanged']}个历史保护文件哈希不变。

另核验{qa['paired_initial_states_recomputed']}组共同初始化、{qa['label_streams_and_weight_changes_checked']}组训练标签及权重变化，以及{qa['OH48_pools_reconstructed']}组OH48背景拒绝实例；各条件内/跨分区exact/RC重复为0、{verified['truth_sites_checked']}个植入位点正确。未做全量近同源审计或独立重新训练，真值仅用于评分/核验。

Python {manifest['python']}、Torch {manifest['dependencies']['torch']}、CPU 4线程，整批拟合连同事前历史文件哈希核验约{manifest['fit_seconds']:.1f}秒；不能把共享成本当成单个模型独立计时。真实Holdout未使用，未执行真实数据拟合、SAE或生物学功能验证。

## 复现与材料

```powershell
python -m experiments.pr01_02_embedding_v2.cnn_boundary_20261027.pipeline fit
python -m experiments.pr01_02_embedding_v2.cnn_boundary_20261027.pipeline score
python -m experiments.pr01_02_embedding_v2.cnn_boundary_20261027.qa
python -m experiments.pr01_02_embedding_v2.cnn_boundary_20261027.pipeline verify
python -m experiments.pr01_02_embedding_v2.cnn_boundary_20261027.report
```

拟合拒绝覆盖旧artifacts。人工序列、模型、实例与聚类对象留本地，公开协议、源码、分数据种子统计、shift/边界记录、图和哈希。严格字节哈希只适用于该环境下原始产物；重训需另开目录并报告环境差异。
'''
    (HERE/'REPORT.md').write_text(text,encoding='utf8');print('CNN boundary report and figure generated')


if __name__=='__main__':main()
