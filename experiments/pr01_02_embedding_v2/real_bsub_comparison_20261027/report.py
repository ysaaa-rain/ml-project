"""Produce figures and summary for the same-data comparison."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from .run import HERE, OUT, sha

def main():
    f=pd.read_csv(OUT/'family.tsv',sep='\t',dtype={'seed':str});p=pd.read_csv(OUT/'cnn_discovery_positions.tsv',sep='\t',dtype={'seed':str})
    s=pd.read_csv(OUT/'cnn_streme_similarity.tsv',sep='\t',dtype={'seed':str});rank=pd.read_csv(OUT/'discovery_rankings.tsv',sep='\t',dtype={'seed':str})
    assert len(f)==14 and len(s)==10
    names=['STREME','CNN_SHIFT','ONEHOT','PROKBERT']; colors=['#447ea5','#58a27c','#d69b43','#a06c96']
    fig,axes=plt.subplots(1,3,figsize=(13,4),sharey=True)
    for ax,col,title in zip(axes,['positive_rate','natural_rate','shuffle_rate'],['Positive groups','Natural negative groups','Paired shuffle groups']):
        for i,(name,color) in enumerate(zip(names,colors)):
            x=f[(f.view=='TOP5')&(f.method==name)][col].to_numpy();ax.bar(i,float(np.mean(x)),color=color);ax.scatter(np.full(len(x),i),x,color='#20252b',s=22,zorder=3)
        ax.set_title(title);ax.set_ylim(0,1);ax.set_xticks(range(4),names,rotation=30,ha='right',fontsize=8);ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Fraction with any selected candidate hit')
    fig.suptitle('Same B. subtilis data and common top-5 scan');fig.tight_layout(rect=(0,0,1,.95));fig.savefig(OUT/'comparison.png',dpi=160);plt.close(fig)
    fig,axes=plt.subplots(2,3,figsize=(14,6),sharex=True,sharey=True)
    for i,seed in enumerate(['20261005','20261006']):
        selected=rank[(rank.method=='CNN_SHIFT')&(rank.seed==seed)].sort_values('rank').head(3)
        for j,row in enumerate(selected.itertuples()):
            pp=p[(p.seed==seed)&(p.candidate==row.candidate)];ax=axes[i,j]
            for strand,color in [('+','#3983a8'),('-','#d48643')]:
                q=pp[pp.strand==strand];ax.bar(q.center_tss,q.instance_n,width=.9,alpha=.8,color=color,label=strand)
            ax.axvline(0,color='#555',ls='--',lw=1);ax.set_xlim(-60,20);ax.set_title(f'{seed} / {row.candidate}');ax.set_xlabel('PWM center relative to TSS (bp)');ax.set_ylabel('Discovery instances')
    fig.legend(handles=axes[0,0].containers[:2],labels=['+','−'],loc='upper center',bbox_to_anchor=(.5,.95),ncol=2,frameon=False)
    fig.suptitle('Discovery candidate locations; strand shown relative to input',y=.995)
    fig.tight_layout(rect=(0,0,1,.88));fig.savefig(OUT/'cnn_positions.png',dpi=160);plt.close(fig)
    lines=['# 枯草芽孢杆菌真实数据：同一规则的候选比较','',
      '基于 MEME/STREME 使用的同一 TJU 冻结数据和划分，比较已有 STREME 自然对照运行的完整10个候选、CNN 平移候选、One-hot 与 ProkBERT 候选。CNN、One-hot、BERT 各2个种子，STREME为1次既有运行；不重新训练。','',
      '所有候选按 Discovery 正类与天然负类的同源组命中差排序，使用共同的天然对照90%分位阈值和双链PWM扫描；下表各取前5个候选，以候选集合至少命中一次计组。天然负类与配对 Shuffle 分开；正类66组、天然对照77组，Shuffle与66组正类配对。图中柱表示已有种子的均值、黑点表示每次运行；这两个种子共享数据，不是两份独立数据。','',
      '![同数据统一扫描](comparison.png)','','| 方法 | 种子 | 正类 | 天然对照 | 差值 | 配对 Shuffle |','| --- | --- | ---: | ---: | ---: | ---: |']
    for r in f[f.view=='TOP5'].itertuples():
        lines.append(f'| {r.method} | {r.seed} | {r.positive_rate:.3f} | {r.natural_rate:.3f} | {r.natural_difference:+.3f} | {r.shuffle_rate:.3f} |')
    lines += ['', '各方法均有较高的 Shuffle 命中，合并多个候选也会提高对照命中，因此需同时看正类与背景。完整候选池和配对正类/Shuffle仅一方命中的组数见 [family.tsv](family.tsv)；单个候选见 [candidates.tsv](candidates.tsv)。CNN使用序列级监督，其他候选产生过程不同；既有 Development 已观察，此表是探索性同数据比较。MEME完整 PWM 当前缺失，尚未纳入。','',
      '## CNN 候选位置与 STREME 局部相似','',
      '![CNN Discovery 实例位置](cnn_positions.png)','','图中为各 CNN Discovery 排名前3个候选的训练实例中心，横轴相对 TSS；正负方向相对已标准化输入。它描述拟合时选中的位置，不是独立定位准确率。','',
      '每种子前5个 CNN 候选与同批 STREME 全部10个候选比较，允许正反向和偏移，最短重叠为较短 PWM 的60%且至少5 bp。最佳相似度和重叠长度见 [cnn_streme_similarity.tsv](cnn_streme_similarity.tsv)，这些 STREME 候选来自同一数据，仅用于描述性参照。','',
      '[Discovery 排序](discovery_rankings.tsv)｜[位置计数](cnn_discovery_positions.tsv)｜[来源与输出哈希](manifest.json)｜[方案](../README.md)。未使用真实 Holdout；未做新的生物功能验证。']
    (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    m=json.loads((OUT/'manifest.json').read_text(encoding='utf8'))
    for name in ['REPORT.md','comparison.png','cnn_positions.png']:m['outputs_sha256'][name]=sha(OUT/name)
    m['report_source_sha256']=sha(Path(__file__))
    (OUT/'manifest.json').write_text(json.dumps(m,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print('Report complete')

if __name__=='__main__':main()
