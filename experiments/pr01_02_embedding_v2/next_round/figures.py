"""Plots only real tables and explicitly labeled algorithm diagrams."""
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from ..audit import BRANCH,ROOT,DATA,load,save,sha
from ..discovery import read_meme
from ..run import logo
from .common import NEXT,SOURCES,PRESENTATION,pool_weights

plt.rcParams.update({'font.family':'Microsoft YaHei','axes.unicode_minus':False,'font.size':11,
                     'axes.spines.top':False,'axes.spines.right':False,'savefig.facecolor':'white'})
NAMES=['枯草芽孢杆菌','鲍曼不动杆菌','慢生根瘤菌','白喉杆菌','大肠杆菌','葡萄球菌']
METHODS=['OH','A0','A1','A2']
COLORS=['#64748b','#334e88','#008e92','#e69138']
INDEX=[]


def export(fig,name,title,source,unit,limitation,kind='真实实验图'):
    dest=NEXT/'figures'/f'{name}.png'
    fig.tight_layout();fig.savefig(dest,dpi=190);plt.close(fig)
    INDEX.append({'figure_id':name,'title':title,'path':str(dest.relative_to(ROOT)),
                  'source_file':str(source.relative_to(ROOT)) if hasattr(source,'relative_to') else source,
                  'script':str((NEXT/'figures.py').relative_to(ROOT)),'manifest':'next_round/run_manifests/',
                  'statistical_unit':unit,'kind':kind,'limitation':limitation,'sha256':sha(dest)})


def main():
    (PRESENTATION/'figures').mkdir(parents=True,exist_ok=True)
    (PRESENTATION/'evidence').mkdir(exist_ok=True)
    # F01: true frozen data counts, metadata only for Holdout.
    source=BRANCH/'audit/counts_gc.tsv';d=pd.read_csv(source,sep='\t')
    fig,ax=plt.subplots(figsize=(10,4.5));bottom=np.zeros(6)
    for split,color in [('discovery','#008e92'),('development','#e69138'),('holdout','#b5c3d4')]:
        vals=[int(d[(d.source_tier.str.startswith(s))&(d.label==1)&(d.split==split)].n.iloc[0]) for s in SOURCES]
        ax.bar(NAMES,vals,bottom=bottom,label=split,color=color);bottom+=vals
    ax.set(ylabel='冻结正样本数量',title='六物种固定划分（仅展示样本数量）');ax.legend(frameon=False,ncol=3)
    export(fig,'F01_data','六物种数据组成',source,'sequence','Holdout只展示冻结数量，不参与本轮评价')
    # F02: definition matters; full-length one-to-one matched PWMs only.
    source=NEXT/'stability_diagnostics/tolerance_stability.tsv';d=pd.read_csv(source,sep='\t')
    filtered=d[(d.strategy=='full')&(d.pwm_similarity>=.8)&(d.method=='prokbert')]
    fig,axes=plt.subplots(1,2,figsize=(10,4.2),sharey=True)
    for ax,stage,title in zip(axes,['discovery_instances','development_scans'],['用于构建PWM的片段','Development独立扫描位点']):
        subset=filtered[filtered.stage==stage]
        g=subset.groupby('tolerance_bp')[['physical_position_jaccard','direction_aware_aligned_jaccard','uniform_null_jaccard_mean']].median()
        ax.plot(g.index,g.physical_position_jaccard,'o-',color='#008e92',label='只匹配位置')
        ax.plot(g.index,g.direction_aware_aligned_jaccard,'s--',color='#334e88',label='位置与方向')
        ax.plot(g.index,g.uniform_null_jaccard_mean,':',color='#a0a8b4',label='均匀随机位置')
        ax.set(title=title,xlabel='允许起点偏移（bp）',ylim=(0,.48));ax.set_xticks([0,1,2,3]);ax.grid(alpha=.15)
    axes[0].set_ylabel('位点集合 Jaccard 中位数');axes[1].legend(frameon=False,fontsize=9)
    export(fig,'F02_stability','精确与容差定位稳定性',source,'matched PWM seed-pair','cos≥0.8稳定家族；各阶段存在共有sequence_id的配对；不是单一生物学样本统计')
    # F03: all signed shifts, including large shifts.
    source=NEXT/'stability_diagnostics/site_offsets.tsv';d=pd.read_csv(source,sep='\t')
    fig,axes=plt.subplots(1,2,figsize=(10,3.7))
    for ax,stage in zip(axes,['discovery_instances','development_scans']):
        x=d[(d.method=='prokbert')&(d.stage==stage)].delta_bp
        ax.hist(x,bins=np.arange(-72.5,73.5,2),density=True,color='#008e92',alpha=.8)
        ax.set(title=stage,xlabel='跨种子起点偏移（bp）',ylabel='描述性密度')
    export(fig,'F03_offsets','定位偏移分布',source,'paired same-sequence instances','位点数有重复家族/种子对，不用于显著性检验')
    # F04: full Discovery/Development GC distributions by species and label.
    source=NEXT/'gc_diagnostics/gc_per_sequence.tsv';d=pd.read_csv(source,sep='\t')
    fig,axes=plt.subplots(2,6,figsize=(18,6.5),sharex=True,sharey=True)
    for row,split in enumerate(['discovery','development']):
        for col,(s,name) in enumerate(zip(SOURCES,NAMES)):
            ax=axes[row,col]
            for label,color in [(1,'#008e92'),(0,'#e69138')]:
                x=d[(d.source==s)&(d.split==split)&(d.label==label)].gc
                ax.hist(x,bins=np.linspace(0,1,18),density=True,histtype='step',lw=1.8,color=color,label=f'{label}, n={len(x)}')
            ax.set(title=f'{name}\n{split}',xlabel='GC比例');ax.legend(frameon=False,fontsize=8)
    axes[0,0].set_ylabel('密度');axes[1,0].set_ylabel('密度')
    export(fig,'F04_gc','六物种GC分布',source,'sequence','观察组成差异，不证明GC导致全部富集差异')
    # F05: effect size, not a ranking by p-values.
    source=NEXT/'gc_diagnostics/motif_gc_adjustment.tsv';a=pd.read_csv(source,sep='\t')
    b=pd.read_csv(NEXT/'gc_diagnostics/shuffle_discordance_diagnostic.tsv',sep='\t')
    d=a.merge(b,on=['source','method','motif'])
    fig,ax=plt.subplots(figsize=(7,4.5))
    for method,color in zip(['onehot','prokbert','tokenwindow','STREME'],COLORS):
        sub=d[d.method==method];ax.scatter(sub.natural_risk_difference,sub.paired_risk_difference,c=color,s=22,alpha=.55,label=method)
    ax.axhline(0,color='#b7bfca',lw=.8);ax.axvline(0,color='#b7bfca',lw=.8)
    ax.set(xlabel='天然对照：出现率差',ylabel='配对Shuffle：出现率差',title='同一固定PWM在两种背景下的效应量');ax.legend(frameon=False)
    export(fig,'F05_controls','天然背景与Shuffle效应量',source,'frozen leakage group','两类背景检验不同；不按显著候选数给方法排序')
    # F06: exact analytic pooling weights, explicitly a method diagram.
    fig,axes=plt.subplots(1,2,figsize=(10,3.7))
    for ax,pos,title in zip(axes,[25,0],['内部窗口 [25,35)：A0与A2重合','序列边缘窗口 [0,10)：权重有差异']):
        for method,color in zip(['A0','A1','A2'],COLORS[1:]):
            w=pool_weights(10,method)[pos];ax.plot(np.arange(76),w,label=method,color=color,lw=2,alpha=.85,ls='--' if method=='A2' else '-')
        ax.set(xlim=(max(0,pos-7),pos+13),title=title,xlabel='6-mer token起点（0-based）',ylabel='归一化聚合权重')
        ax.legend(frameon=False)
    export(fig,'F06_pooling','聚合权重公式核验',NEXT/'pooling_ablation/weight_equivalence.tsv','analytic token weight','不是实验发现位点；hidden state仍含全序列上下文','原理示意')
    # F07: the actual strict, independent synthetic test matrix.
    source=NEXT/'synthetic_benchmark/strict_metrics.tsv';d=pd.read_csv(source,sep='\t')
    x=d[(d.split=='test')&(d.tolerance_bp==2)&(((d.width==10)&~d.scenario.str.startswith(('S5','S6')))|((d.width==6)&d.scenario.str.startswith(('S5','S6'))))]
    means=x.groupby(['scenario','method']).f1.mean().unstack().reindex(columns=METHODS)
    fig,ax=plt.subplots(figsize=(8.4,6));im=ax.imshow(means.values,vmin=0,vmax=1,cmap='YlGnBu')
    ax.set_xticks(range(4),['OH','A0','A1','A2']);ax.set_yticks(range(len(means)),means.index)
    ax.set_xlabel('OH: One-hot   A0: 原聚合   A1: 窗口内Token   A2: 覆盖加权',fontsize=9)
    for (i,j),v in np.ndenumerate(means.values):ax.text(j,i,f'{v:.3f}',ha='center',va='center',color='white' if v>.65 else '#182b4d',fontsize=11)
    ax.set_title('通过参考PWM门槛后的 ±2bp 定位 F1\n独立合成Test，三个种子均值');fig.colorbar(im,ax=ax,fraction=.04,label='F1')
    export(fig,'F07_synthetic','合成难度与方法效果矩阵',source,'synthetic seed × motif family','0包含未过PWM匹配门槛；不能推断hidden state完全无信息；不是分类准确率')
    # F08: window width dependence, error bar is observed seed range.
    sub=d[(d.scenario=='S3_strong_random')&(d.split=='test')&(d.tolerance_bp==2)]
    fig,ax=plt.subplots(figsize=(8,4));xx=np.arange(3)
    for j,(method,color) in enumerate(zip(METHODS,COLORS)):
        g=sub[sub.method==method].groupby('width').f1.agg(['mean','min','max']).reindex([6,10,16])
        ax.bar(xx+(j-1.5)*.18,g['mean'],width=.18,color=color,label=method,
               yerr=np.vstack([g['mean']-g['min'],g['max']-g['mean']]),capsize=2)
    ax.set_xticks(xx,[6,10,16]);ax.set(xlabel='窗口宽度（bp）',ylabel='±2bp定位F1',ylim=(0,1.17),title='随机位置强模式的窗口消融');ax.legend(frameon=False,ncol=4,loc='upper center')
    export(fig,'F08_width','窗口宽度消融',source,'synthetic seed','误差线是三个种子的实际范围，不是置信区间')
    # F09: synthetic pair grammar. Display true paired observations only.
    fig,axes=plt.subplots(1,2,figsize=(10,4));gap_rows=[]
    for ax,scenario in zip(axes,['S5_pair_fixed','S6_pair_variable']):
        for method,color in zip(METHODS,COLORS):
            paths=list((NEXT/'synthetic_benchmark').glob(f'{scenario}*/*__w6/strict/grammar.tsv'))
            chunks=[]
            for p in paths:
                if p.parent.parent.name.startswith(method+'__'):chunks.append(pd.read_csv(p,sep='\t'))
            if not chunks:continue
            g=pd.concat(chunks);g=g[g.split=='test'];valid=g[g.both_detected]
            ax.scatter(valid.true_gap,valid.predicted_gap,s=7,alpha=.15,color=color,label=f'{method} n={len(valid)}')
            gap_rows.append({'scenario':scenario,'method':method,'test_pairs':len(g),'both_detected':len(valid),
              'correct_order_rate_all_pairs':float(g.correct_order.mean()),
              'gap_within3_rate_all_pairs':float((g.gap_error.abs()<=3).mean())})
        ax.plot([0,30],[0,30],color='#a2acb9',ls=':');ax.set(title=scenario,xlabel='真值间距（bp）',ylabel='预测间距（bp）');ax.legend(fontsize=8,frameon=False)
    pd.DataFrame(gap_rows).to_csv(NEXT/'synthetic_benchmark/grammar_summary.tsv',sep='\t',index=False)
    export(fig,'F09_pair','双Motif间距恢复',NEXT/'synthetic_benchmark/grammar_summary.tsv','independent synthetic test promoter','仅展示双检测者；总体分母与漏检率另表保留，不声称真实生物学Grammar')
    # F10: true traditional-PWM Development positions and gaps, not strict FIMO results.
    source=NEXT/'figures/traditional_development_positions.tsv';h=pd.read_csv(source,sep='\t')
    gaps=pd.read_csv(NEXT/'figures/traditional_development_pairs.tsv',sep='\t')
    fig,axes=plt.subplots(1,2,figsize=(10,4))
    s=SOURCES[0];pos=h[h.source==s]
    for name,g in pos.groupby('motif'):axes[0].hist(g.center_tss,bins=np.arange(-60,22,4),histtype='step',lw=1.8,label=name.removeprefix('STREME_'))
    axes[0].axvline(0,color='#a0a8b4',ls=':');axes[0].set(xlabel='Motif中心相对TSS（bp）',ylabel='命中序列数',title='枯草芽孢杆菌固定STREME PWM');axes[0].legend(fontsize=8,frameon=False)
    axes[1].hist(gaps[gaps.source==s].gap,bins=np.arange(-20,72,4),color='#008e92',alpha=.8)
    axes[1].set(xlabel='同序列不同主位点 edge gap（bp）',ylabel='配对数',title='同序列不同Motif的描述性间距')
    export(fig,'F10_grammar','传统PWM位置与间距描述',source,'Development sequence / motif pair','自定义最大分扫描，非严格FIMO q；不解释为−35/−10机制或显著组合')
    # F11: one real, deterministically selected instance with unstable extraction.
    source=NEXT/'stability_diagnostics/site_offsets.tsv';offset=pd.read_csv(source,sep='\t')
    sub=offset[(offset.source==SOURCES[0])&(offset.method=='prokbert')&(offset.stage=='discovery_instances')&(offset.delta_bp.abs()>10)]
    if len(sub):
        r=sub.sort_values(['sequence_id','seed_a','seed_b','motif_a']).iloc[0]
        frame=load(SOURCES[0],'discovery',1).set_index('sequence_id');sequence=frame.loc[r.sequence_id].sequence
        fig,ax=plt.subplots(figsize=(14,3))
        for i,c in enumerate(sequence):ax.text(i,1.0,c,ha='center',va='center',fontfamily='Consolas',fontsize=9)
        for y,pos,label,color in [(0.4,r.start_a0,f'种子 {r.seed_a} / {r.motif_a}',COLORS[1]),(-.2,r.start_b0,f'种子 {r.seed_b} / {r.motif_b}',COLORS[2])]:
            ax.plot([pos,pos+9],[y,y],lw=10,color=color,solid_capstyle='butt');ax.text(0,y-.2,label,fontsize=9)
        ax.axvline(60,color='#e69138',ls='--');ax.text(60,1.5,'TSS（索引60）',ha='center',fontsize=10)
        ax.set(xlim=(-2,83),ylim=(-.7,1.8),xlabel='输入序列索引（0-based）',title=str(r.sequence_id));ax.set_yticks([])
        export(fig,'F11_instance','真实构建实例对照',source,'same frozen promoter','按字典序取首个偏移>10bp的示例，属于失败案例，不代表总体分布')
        save(NEXT/'figures/instance_example.json',r.to_dict())
    # F12: a genuine original STREME logo.
    p=ROOT/'results/motif/tju_streme_v2_20261008'/SOURCES[0]/'selected_main.meme'
    matrices=read_meme(p);name=next(iter(matrices));dest=NEXT/'figures/F12_streme_logo.png';logo(dest,matrices[name],name)
    INDEX.append({'figure_id':'F12_streme_logo','title':'原固定STREME候选Logo','path':str(dest.relative_to(ROOT)),
                  'source_file':str(p.relative_to(ROOT)),'script':'first-round run.logo','statistical_unit':'discovery candidate PWM',
                  'kind':'真实实验图','limitation':'候选及历史注释相似，不是功能确认','sha256':sha(dest)})
    pd.DataFrame(INDEX).to_csv(PRESENTATION/'evidence/figure_index.csv',index=False,encoding='utf-8-sig')
    print(f'Exported {len(INDEX)} figures',flush=True)


if __name__=='__main__':main()
