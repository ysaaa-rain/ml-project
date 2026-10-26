"""Budget-aware, data-seed-level report of fresh background rejection validation."""
import numpy as np
import pandas as pd

from .pipeline import HERE,PROTOCOL,SUMMARY,read


def table(frame):
    frame=frame.copy()
    for key in frame:
        if frame[key].dtype.kind=='f':frame[key]=frame[key].map(lambda v:f'{v:.3f}')
    return '\n'.join(['| '+' | '.join(map(str,frame.columns))+' |',
        '| '+' | '.join(['---']*len(frame.columns))+' |',
        *['| '+' | '.join(map(str,row))+' |' for row in frame.itertuples(index=False,name=None)]])


def figures(seeds):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    group=seeds[(seeds.scenario=='mut30_random')&(seeds.variant=='BG_ONLY')]
    fig,axes=plt.subplots(1,2,figsize=(12,4.8))
    views=['FULL','CAP5','MATCH_OH24','MATCH_OH48']
    styles=[('EMB128','#d06426',-.2),('OH24','#2873af',0),('OH48','#2c8c63',.2)]
    for axis,metric,label in [(axes[0],'f1','Reference-assisted localization F1 (±2bp)'),
                              (axes[1],'family_blank_fpr','Any-candidate blank FPR')]:
        for representation,color,offset in styles:
            for index,view in enumerate(views):
                values=group[(group.representation==representation)&(group.view==view)][metric].to_numpy()
                if not len(values):continue
                axis.scatter(index+offset+np.linspace(-.035,.035,len(values)),values,s=25,color=color,alpha=.85,
                             label=representation if index==0 else None)
                axis.plot([index+offset-.08,index+offset+.08],[values.mean()]*2,color=color,lw=2.5)
        axis.set_xticks(range(4),views,rotation=18,ha='right');axis.set(ylabel=label,ylim=(-.03,1.04));axis.grid(axis='y',alpha=.2)
    axes[0].legend(loc='upper right',fontsize=9)
    fig.suptitle('30% mutation, BG_ONLY: five data seeds; each dot averages three clustering seeds',fontsize=10)
    fig.tight_layout();fig.savefig(SUMMARY/'F1_budget_validation.png',dpi=170);plt.close(fig)


def main():
    cfg=read(PROTOCOL);manifest=read(SUMMARY/'manifest.json');verified=read(SUMMARY/'verification.json')
    qa=read(SUMMARY/'independent_qa.json');decision=read(SUMMARY/'primary_decision.json')
    assert manifest['status']=='COMPLETED' and verified['status']==qa['status']=='PASSED'
    summary=pd.read_csv(SUMMARY/'test_summary.tsv',sep='\t');seeds=pd.read_csv(SUMMARY/'per_data_seed.tsv',sep='\t')
    paired=pd.read_csv(SUMMARY/'paired_deltas.tsv',sep='\t');effects=pd.read_csv(SUMMARY/'paired_summary.tsv',sep='\t')
    figures(seeds)
    primary=paired[(paired.scenario=='mut30_random')&(paired.variant=='BG_ONLY')&(paired.view=='MATCH_OH48')&(paired.comparator=='OH48')]
    primary=primary[['data_seed','f1_emb','f1_baseline','delta','family_fpr_emb','family_fpr_baseline','budget_emb','budget_baseline','zero_budget_count']]
    comparison=[]
    for (variant,view,representation),group in summary.groupby(['variant','view','representation']):
        row={'variant':variant,'view':view,'representation':representation}
        for scenario in cfg['scenarios']:
            r=group[group.scenario==scenario['id']].iloc[0]
            row[scenario['id']]=f'{r.f1_mean:.3f} [{r.f1_data_min:.3f}, {r.f1_data_max:.3f}]'
        comparison.append(row)
    main_summary=summary[(summary.scenario=='mut30_random')&(summary.variant=='BG_ONLY')]
    gate_labels={'matched_mean_positive':'匹配预算的平均F1差为正',
                 'cap5_mean_positive':'CAP5平均F1差为正',
                 'at_least4_of5_wins':'至少4/5组数据胜出',
                 'no_zero_common_budget':'共同预算不为零',
                 'family_fpr_not_above_OH48':'候选集误报不高于One-hot48'}
    gate_table=pd.DataFrame([{'预设检查':gate_labels[name],'结果':'通过' if passed else '未通过'}
                             for name,passed in decision['gate_checks'].items()])
    status='通过预设进一步研究门槛' if decision['advance_gate_passed'] else '未通过预设进一步研究门槛'
    text=f'''# 背景拒绝的独立复核与候选预算匹配

日期：2026-10-26。协议：{cfg['version']}。本轮已完成合成验证，未启用真实Holdout或功能验证。

## 结论与后续指导

上一轮的小样本均值优势没有在本轮独立数据上复现。主要比较中双方每次均保留5个候选：BERT的定位F1为0.295，完整One-hot48为0.543，配对均值差−0.248；BERT在5/5组数据中均较低。候选集空白误报均值为0.163对0.146，BERT在5组数据中也均较高。FULL视图的平均F1差距仍存在，不能把结论归因于预算裁剪或零预算；主要比较没有候选不足。

背景拒绝确实改善BERT自身：CAP5下，30%退化从RAW的0.083升至0.295，混合反向链15%退化从0.385升至0.679；但后者仍低于One-hot24的0.861和One-hot48的0.841。One-hot48的混合链RAW为0.879，背景拒绝也不是所有表示的通用改进。强模式控制三种表示均约0.974，与加载和维度核验一起，未显示简单的模型失载或对照维度错误。

据此，建议收束当前“冻结ProkBERT-mini局部表示→PCA/K-Means→背景拒绝”的常规调参，将其保留为补充探索与失败机制材料；主研究继续依托简单基线和传统motif路线。只有提出与本轮明显不同、可证伪的表示或发现机制假设，再开新协议与新Test；本轮不直接启动更多编码器、SAE、微调或真实数据拟合。这里的负面结果只针对当前流程和合成任务，不能推广为所有BERT/DNA语言模型无效，也不构成传统方法已获胜的实测证明。

## 预设问题与完成范围

上一轮BERT BG_ONLY的30%退化均值约0.482，高于同处理One-hot约0.418，但仅1/3数据种子胜出，候选数与误报也不同。本轮检验这个优势在新数据、更充分算法重复及候选预算匹配下是否仍存在，不复用旧Test。

3类条件（强模式、30%退化、混合反向链15%退化）×5新数据种子×3聚类种子，比较BERT局部有序128维、One-hot PCA24及不降维One-hot48。RAW/BG_ONLY均运行。实际{manifest['cluster_fits']}次聚类、{manifest['pwm_pools']}组PWM池；从固定池形成{manifest['evaluation_sets']}组预算评价集合，不把子集重复计为新拟合。

每条件Discovery256、Development128、Test128条植入序列，各配同数空白，另有独立校准空白256条。81bp序列、连续12bp人工模式和12bp窗口；BERT权重冻结，局部上下文24bp。人工模式不代表真实−35/−10间距结构。未再跑普通15%退化，聚焦上一轮主要疑问及强模式控制。

数据种子202610263–267，聚类种子20261030/31/32；PCA种子固定。先在每个数据种子内平均3次聚类，再汇总5组数据。下表范围为5个数据均值min/max，不是置信区间，15次运行不是15个独立样本。

同一种子的不同难度条件复用空白随机流，条件之间也不是额外独立重复。exact/RC重复核查的范围为每个条件内部的各分区，不主张跨所有难度条件的每条序列都唯一。

## 排名与预算冻结

候选先按Discovery正类减空白命中率排序，分数阈值为对应PWM的Discovery空白序列最大分95百分位；平局按相对背景的信息量、支持序列数降序及ID升序。没有参考PWM、植入位置或Test参与排名。背景拒绝仍按RAW PWM的Discovery空白最大分95百分位严格筛选实例，至少20条重建一次；不叠加对齐。

- FULL：全部候选，显示完整流程能力。
- CAP5：各自最多5个；候选不足不补造PWM，标记underfill。
- MATCH_OH24 / MATCH_OH48：BERT与对应基线各取同样m个，m=min(5,双方可用数)。m=0保留0分，不删除或算胜出。

共同裁剪可能掩盖候选不足的方法，因此必须同时看FULL/CAP5。所有拟合、Discovery排名、候选成员名单和预算均在任何新Development/Test评分前锁定；本轮未按Test调整方法。

## 主要终点与预设门槛

主要终点为30%退化、BG_ONLY、MATCH_OH48的BERT减完整One-hot48定位F1。状态：**{status}**。均值差{decision['mean_delta']:.3f}，5组数据胜出{decision['wins']}组，CAP5均值差{decision['cap5_mean_delta']:.3f}。

{table(gate_table)}

门槛要求匹配和CAP5均值差都正、至少4/5数据种子胜出、共同预算不为零、匹配后候选集空白误报不高于One-hot48；它仅决定是否值得更大复核，不是总体/功能优势证明。

{table(primary)}

配对均值符号翻转参照p={decision['paired_sign_flip_p']:.4f}，枚举32种符号分配。其推断依赖零假设下配对差值符号可交换；5组数据的双侧p最小0.0625，因此不能宣称p<.05。其他条件、预算和OH24对照均为探索性。

## 30%退化的预算、误报和覆盖

{table(main_summary.drop(columns=['scenario','variant']))}

![Budget-controlled validation](summary/F1_budget_validation.png)

budget为每个数据种子内先平均3次算法运行的实际候选数，再跨数据平均。零预算和不足5个的运行没有被删掉。误报是Test空白中任一候选超过自身固定阈值的经验比例，不是正式FDR；数量匹配仍不保证候选相关结构或发现难度相同。

## 全部预设结果

{table(pd.DataFrame(comparison))}

评估沿用全12列参考cosine≥.8、方向正确、±2bp参考辅助定位F1。参考只在已冻结的候选子集内作家族匹配，不能找回预算排除的候选。阈值来自独立校准空白最大分95百分位，严格大于；不是FIMO位点p/q，也不是无监督绝对precision。

## 配对比较

{table(effects)}

完整逐数据种子差值见 `paired_deltas.tsv`，算法运行和分母见 `metrics.tsv`。不把同一数据的3个算法种子当额外生物学重复；不按Test选择最佳预算、方法或候选。更换测试版本和加入One-hot48后，不能与上一轮均值拼成连续提升曲线。

## 核验与资源

37项相关测试通过，覆盖共同预算、候选不足、禁止参考越过子集、误报不被参考门槛隐藏、完整评分流程与符号翻转边界。{verified['hashes_checked']}项产物哈希、{verified['instances_checked']}条实例、{verified['pwms_rebuilt']}个重建PWM、{verified['background_subset_pairs']}组BG子集关系及900个冻结预算集合核验通过；{verified['protected_files_unchanged']}个历史文件哈希不变。

额外重算{qa['discovery_rankings_recomputed']}组Discovery候选排名，核验{qa['truth_sites_checked']}个植入位置及{qa['synthetic_sequences_checked']}条序列记录，每个条件内部的分区内/间exact/RC重复为0；没有全量近同源审计。以上是产物一致性，不等于独立重新拟合模型。

保存的聚类模型输入维度核验为BERT128=128、One-hot24=24、One-hot48=48；后者确实没有PCA降维。图表已作视觉检查。

Python {manifest['python']}，Torch {manifest['dependencies']['torch']} CPU，Transformers {manifest['dependencies']['transformers']}；4线程，无GPU。共享拟合约{manifest['fit_elapsed_seconds']:.1f}秒，不把共享编码成本当单方法独立时间。固定模型版本/哈希和无缺失权重加载见 `encoder.json`。

STREME环境仍未补齐，不能宣称总体胜过传统方法。真实Holdout未用，SAE、微调、新编码器和真实数据拟合均未执行。若依据本轮Test改方法，另开版本并更换Test种子。

## 复现

```powershell
python -m experiments.pr01_02_embedding_v2.budget_validation_20261026.pipeline fit
python -m experiments.pr01_02_embedding_v2.budget_validation_20261026.pipeline score
python -m experiments.pr01_02_embedding_v2.budget_validation_20261026.qa
python -m experiments.pr01_02_embedding_v2.budget_validation_20261026.pipeline verify
python -m experiments.pr01_02_embedding_v2.budget_validation_20261026.report
```

拟合拒绝覆盖已有目录。权重、特征、人工序列、实例与拟合对象留本地；公开代码、冻结协议、汇总、预算/排名审计、图表及哈希，保留失败与零候选结果。
'''
    (HERE/'REPORT.md').write_text(text,encoding='utf-8');print('Budget validation report and figure generated')


if __name__=='__main__':main()
