"""Report seed-level experimental evidence without treating algorithm runs as samples."""
import numpy as np
import pandas as pd

from .pipeline import HERE,PROTOCOL,SUMMARY,read

METHODS=['OH_P24','LOCAL_P24','LOCAL_P64','LOCAL_P128','LOCAL_P24_CORE50','LOCAL_P128_CORE50']
LABELS=['One-hot24','Local24','Local64','Local128','Local24+Core50','Local128+Core50']


def table(frame,digits=3):
    frame=frame.copy()
    for column in frame:
        if frame[column].dtype.kind=='f': frame[column]=frame[column].map(lambda x:f'{x:.{digits}f}')
    header='| '+' | '.join(map(str,frame.columns))+' |'
    separator='| '+' | '.join(['---']*len(frame.columns))+' |'
    rows=['| '+' | '.join(map(str,row))+' |' for row in frame.itertuples(index=False,name=None)]
    return '\n'.join([header,separator,*rows])


def figures(per_seed,probes):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    scenarios=['strong_random','mut15_random','mut30_random']
    fig,axes=plt.subplots(1,3,figsize=(14,4.5),sharey=True)
    for ax,scenario in zip(axes,scenarios):
        group=per_seed[(per_seed.tier=='main')&(per_seed.scenario==scenario)]
        for index,method in enumerate(METHODS):
            values=group[group.method==method].f1.to_numpy()
            x=index+np.linspace(-.14,.14,len(values))
            ax.scatter(x,values,s=23,alpha=.8,color='#2563a6')
            ax.plot([index-.22,index+.22],[values.mean()]*2,color='#c34931',lw=2)
        ax.set_xticks(range(len(METHODS)),LABELS,rotation=42,ha='right')
        ax.set(title=scenario,ylim=(-.03,1.04)); ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Reference-assisted localization F1 (±2 bp)')
    fig.suptitle('Dots: five data seeds (three clustering seeds averaged); red: mean',fontsize=11)
    fig.tight_layout(); fig.savefig(SUMMARY/'F1_seed_validation.png',dpi=170);plt.close(fig)
    reps=['OH_RAW','LOCAL_RAW','LOCAL_P24','LOCAL_P64','LOCAL_P128']
    fig,ax=plt.subplots(figsize=(8.5,4.5))
    for index,representation in enumerate(reps):
        values=probes[probes.representation==representation].top1_within2.to_numpy()
        ax.scatter(index+np.linspace(-.12,.12,len(values)),values,s=35,color='#2563a6')
        ax.plot([index-.2,index+.2],[values.mean()]*2,color='#c34931',lw=2)
    ax.set_xticks(range(len(reps)),reps);ax.set(ylim=(-.03,1.04),ylabel='Test top1 within ±2bp',
        title='Supervised ranking probes on 30% mutation (not motif discovery)')
    ax.grid(axis='y',alpha=.2);fig.tight_layout();fig.savefig(SUMMARY/'F2_supervised_probes.png',dpi=170);plt.close(fig)


def main():
    cfg=read(PROTOCOL); manifest=read(SUMMARY/'manifest.json'); verification=read(SUMMARY/'verification.json')
    replication=read(SUMMARY/'replication_manifest.json')
    assert manifest['status']=='COMPLETED' and verification['status']=='PASSED' and replication['status']=='COMPLETED'
    summary=pd.read_csv(SUMMARY/'test_summary.tsv',sep='\t')
    seeds=pd.read_csv(SUMMARY/'per_data_seed.tsv',sep='\t')
    probes=pd.read_csv(SUMMARY/'supervised_probes.tsv',sep='\t')
    recipes=pd.read_csv(SUMMARY/'replication_summary.tsv',sep='\t')
    figures(seeds,probes)
    comparison=pd.DataFrame({'method':METHODS})
    for scenario in cfg['scenarios']:
        indexed=summary[(summary.tier=='main')&(summary.scenario==scenario['id'])].set_index('method')
        comparison[scenario['id']]=[f'{indexed.loc[m,"f1_mean"]:.3f} [{indexed.loc[m,"f1_data_min"]:.3f}, {indexed.loc[m,"f1_data_max"]:.3f}]' for m in METHODS]
    probe_summary=probes.groupby('representation').agg(top1_within2_mean=('top1_within2','mean'),
        top1_within2_min=('top1_within2','min'),top1_within2_max=('top1_within2','max'),
        exact_mean=('top1_exact','mean'),data_seeds=('data_seed','nunique'),max_iterations=('iterations','max')).reset_index()
    stability=pd.read_csv(SUMMARY/'algorithm_stability.tsv',sep='\t')
    stability=stability[(stability.tier=='main')&(stability.cosine>=.8)].groupby('method').agg(
        matched_pairs=('cosine','size'),pwm_cosine_median=('cosine','median'),
        instance_jaccard_median=('instance_jaccard','median')).reset_index()
    delta=pd.read_csv(SUMMARY/'paired_delta_vs_OH.tsv',sep='\t')
    delta_rows=[]
    for scenario in [s['id'] for s in cfg['scenarios']]:
        group=delta[(delta.tier=='main')&(delta.scenario==scenario)]
        for method in METHODS[1:]:
            values=group[method]
            delta_rows.append({'scenario':scenario,'method':method,'mean_delta_vs_OH':values.mean(),
                'data_seeds_better_than_OH':int((values>0).sum()),'data_seeds':len(values)})
    stress=summary[summary.tier=='stress'].drop(columns=['tier'])
    nulls=pd.read_csv(SUMMARY/'candidate_null_fpr.tsv',sep='\t')
    null_summary=nulls[nulls.split=='test'].groupby(['tier','method']).raw_blank_fpr.agg(['mean','min','max']).reset_index()
    text=f'''# 有序局部表示复核与PCA/实例瓶颈：阶段报告

日期：2026-10-25。协议：{cfg['version']}。状态：本轮合成复核、消融、压力测试及定位探针完成；未进入真实数据与功能验证。

## 实际执行与证据边界

主实验3难度×5全新数据种子×3聚类种子，6种发现配置；另3种压力条件各2数据种子×3聚类种子×3配置。主矩阵完成{manifest['cluster_fits']}次实际聚类、{len(manifest['runs'])}组PWM输出。CORE50复用聚类，不重复算拟合。另以原版fit函数完成90次原流程复核，合计324次实际聚类、414组发现输出。监督定位探针25个，单独统计。

先冻结协议和源码、拟合并锁定全部PWM，再读取Test评价；未依据本轮Test修改任何发现方法、门槛或候选。原流程复核也有独立协议、fit lock与源码哈希。上一轮Test不复用。真实Holdout未启用，未新增功能motif结论。

每条件Discovery256、Development128、Test128条植入序列，各配同数空白，另校准空白256条。序列81bp，人工连续模式和核心窗口均12bp；不代表真实−35/−10间距结构。压力条件是GC=.75强模式、一阶持续性=.7相关背景强模式、50%概率反向链的15%退化模式。生成时核验分区内/间无exact/RC重复，未执行全量近同源审计。

## 可比性与统计单位

局部编码仍是固定ProkBERT-mini的L6，24bp上下文、核心7个6-mer hidden states按位置拼接。每条序列70个窗口对应58个不同裁剪，复用边缘裁剪只减少重复计算，坐标核验不变。

PCA消融固定一个Discovery PCA128基底，取前24/64/128维，固定PCA种子；只改变KMeans种子。CORE50使用同一批聚类分配，在每簇按Discovery距离中位数保留窗口，再每序列取最近窗口，至少20条形成PWM，不使用植入真值。没有进行片段偏移/反向互补对齐。

另列上一轮原流程复核：直接调用未修改的 `bottleneck_20261025.pilot.fit`，每个算法种子重新拟合PCA24与KMeans；One-hot与有序局部表示均跑。这样避免把新的嵌套PCA计算路径误称为原流程逐项复现。裁剪batch96与历史64不同，不主张编码结果逐字节一致。

主统计**先平均同一数据种子内的3个聚类运行，再汇总5个数据种子**；不是15个独立样本。压力测试仅两个数据种子，只描述。以下范围是5个数据均值的min/max，不是置信区间。

## 原流程复核

{table(recipes)}

## PCA与实例筛选主实验

{table(comparison)}

![Seed-level validation](summary/F1_seed_validation.png)

每个数据种子相对One-hot的配对变化：

{table(pd.DataFrame(delta_rows))}

所有算法运行、阈值门槛、参考cosine与TP/预测分母保留在 `summary/metrics.tsv`；逐数据种子表见 `per_data_seed.tsv`。主要F1要求全宽参考cosine≥.8、方向正确、起点误差≤2bp；未过参考门槛记0，不能推断hidden state完全无信息。参考辅助选最佳家族，不是无人监督实际precision。

## 压力条件

{table(stress)}

混合反向链评价使用真实植入方向与PWM参考方向的组合，不能只按物理位置给TP。当前发现仍按正向窗口聚合，没有实例方向对齐，因此压力失败可指向流程短板，但不是单因素因果证明。

## 监督定位探针：30%退化

{table(probe_summary)}

![Supervised probes](summary/F2_supervised_probes.png)

探针只在Discovery已知植入位置训练；每条序列1个准确核心窗和4个无重叠背景窗，StandardScaler加balanced LogisticRegression(C=1,max_iter=1000)。固定训练设置，Test只在全部无监督拟合锁定后编码，枚举70窗取最高分位置。报告Test植入正序列内top1定位率，没有空白检测/FPR评价，**不是无监督发现、不是生物学证据，也不能与上表PWM定位F1等同**。原始高维与PCA探针还受监督正则化、维度和样本量影响，不能仅凭差异证明单一因果。

## 跨聚类种子的描述性稳定性

{table(stability)}

全宽cosine≥.8的一对一匹配家族，报告PWM相似与实例精确起点Jaccard；包含背景家族，坐标不作方向对齐。高PWM相似不等于实例稳定，算法配对不是独立数据重复，不据此计算生物学显著性。

## 全候选空白误报

{table(null_summary)}

阈值为独立合成校准空白最大分95百分位，严格score>threshold；上述FPR包含未通过参考匹配的所有候选，不把失败候选藏掉。每候选单独校准，均值不等于整个候选家族误报率；不是FIMO位点p/q。

## 核验与资源

20项相关测试通过。{verification['hashes_checked']}项产物哈希、{verification['instances_checked']}条实例、{verification['truth_sites_checked']}条植入真值、{verification['crop_core_checks']}次核心还原核查通过；{verification['protected_files_unchanged']}个历史保护文件不变。原流程复核输出另有独立 `replication_lock.json`，原fit源码哈希保留。哈希一致不等于独立重跑历史实验。

Python {manifest['python']}，Torch {manifest['dependencies']['torch']} CPU，Transformers {manifest['dependencies']['transformers']}，4线程，无GPU。主矩阵拟合约{manifest['fit_elapsed_seconds']:.1f}秒，后续评价、探针和复核分别执行，不将混合墙钟时间当单一方法成本。编码器revision/文件哈希与无缺失权重加载证据在 `summary/encoder.json`。

## 尚未闭合与真实数据门槛

当前环境无STREME、Docker或WSL发行版，见 `tool_availability.json`；没有用替代实现冒充STREME。传统基线未补齐，任何合成优势都不能直接写成总体胜过传统方法。片段偏移/方向对齐、可变宽度、更多人工模式、真实数据注释验证仍未执行。SAE、微调和更换编码器继续后置。

若依据本轮Test继续改方案，另开版本并更换Test种子。只有恢复能力、压力条件、实例稳定性及传统基线得到支持后，才考虑真实数据小试；真实Development已反复使用，项目Holdout曾被传统主线查看，均不能称为全新盲验证。

## 复现

```powershell
python -m experiments.pr01_02_embedding_v2.validation_20261025.pipeline fit
python -m experiments.pr01_02_embedding_v2.validation_20261025.replicate fit
python -m experiments.pr01_02_embedding_v2.validation_20261025.pipeline score
python -m experiments.pr01_02_embedding_v2.validation_20261025.replicate score
python -m experiments.pr01_02_embedding_v2.validation_20261025.pipeline probes
python -m experiments.pr01_02_embedding_v2.validation_20261025.pipeline verify
python -m experiments.pr01_02_embedding_v2.validation_20261025.report
```

拟合拒绝覆盖已有目录。重跑应使用新checkout或新版本。模型/特征/人工序列/实例/拟合对象留本地；公开配置、源码、汇总、图表、诊断与哈希。保护原始输入、旧结果和用户已有无关改动。
'''
    (HERE/'REPORT.md').write_text(text,encoding='utf-8')
    print('Validated report and figures generated')


if __name__=='__main__': main()
