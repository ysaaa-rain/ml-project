"""Describe the uniform scan of real-data PWM candidates."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/same_batch_real_equal_20261027'

def main():
    runs=pd.read_csv(OUT/'runs.tsv',sep='\t');s=pd.read_csv(OUT/'summary.tsv',sep='\t')
    assert len(runs)==42 and len(s)==18
    assert (runs.candidate_count==runs.common_budget).all()
    for _,g in runs.groupby('source'):assert g.common_budget.nunique()==1
    labels={'tjupan_bacillus_subtilis':'B. subtilis','tjupan_baumannii':'A. baumannii','tjupan_bradyrhizobium':'Bradyrhizobium',
            'tjupan_diphtheria':'C. diphtheriae','tjupan_escherichia_coli':'E. coli','tjupan_staphylococcus':'Staphylococcus'}
    methods=['STREME','onehot','prokbert']; colors=['#477eab','#dd9e48','#4c9b70']
    fig,axes=plt.subplots(1,2,figsize=(14,5),sharey=True)
    for ax,col,title in zip(axes,['positive_rate','control_rate'],['Positive sequence hit rate','Control sequence hit rate']):
        for i,(method,color) in enumerate(zip(methods,colors)):
            a=s[s.method==method].set_index('source').loc[list(labels)]
            ax.bar(np.arange(6)+(i-1)*.24,a[col],width=.23,label=method,color=color)
        ax.set_title(title);ax.set_xticks(range(6),list(labels.values()),rotation=30,ha='right');ax.set_ylim(0,1);ax.grid(axis='y',alpha=.2)
    axes[0].set_ylabel('Fraction with at least one selected PWM hit');axes[1].legend()
    fig.suptitle('Same species, same split and common candidate budget');fig.tight_layout();fig.savefig(OUT/'comparison.png',dpi=160);plt.close(fig)
    lines=['# 六物种真实数据：同批候选统一扫描','',
      '比较已有STREME主候选、One-hot和ProkBERT候选，使用同一物种、同一Discovery/Development划分和相同天然对照。每个物种按三种方法所有运行中共同可用的候选数取相同数量，最多5个；因此不同物种的共同数量可以不同。STREME为1次既有候选集，One-hot与ProkBERT各为3个种子、L6W10。','',
      '统一在Discovery的天然对照上确定每个PWM的95百分位扫描阈值，双向扫描；候选按Discovery正类命中率减对照命中率排序，信息量和名称用于平局排序。对Development序列，统计候选集合至少命中一次的正类及对照比例。图中对照命中不是可以忽略的噪声：应与正类命中一起阅读。','',
      '这是对已有候选及已观察过的Development数据进行的事后探索。STREME候选是原有已筛选主候选，Embedding候选来自各运行的完整PWM池；候选产生与初筛方式不同。结果不等于准确率或功能验证。CNN和MEME尚没有在这批真实数据上按此方案生成可比较的候选，故不在表中虚填。','',
      '![真实数据命中比较](comparison.png)','','| 物种 | 方法 | 共同候选数 | 正类命中率 | 对照命中率 | 差值 |','| --- | --- | ---: | ---: | ---: | ---: |']
    for source in labels:
        for method in methods:
            r=s[(s.source==source)&(s.method==method)].iloc[0]
            lines.append(f'| {labels[source]} | {method} | {r.candidate_count:.0f} | {r.positive_rate:.3f} | {r.control_rate:.3f} | {r.rate_difference:+.3f} |')
    lines += ['', '[运行明细](runs.tsv)｜[选中候选](candidates.tsv)｜[汇总](summary.tsv)｜[来源与校验](manifest.json)。']
    (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    m=json.loads((OUT/'manifest.json').read_text(encoding='utf8'))
    for name in ['REPORT.md','comparison.png']:
        m['output_sha256'][name]=hashlib.sha256((OUT/name).read_bytes()).hexdigest()
    m['report_script_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (OUT/'manifest.json').write_text(json.dumps(m,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print('real report complete')

if __name__=='__main__':main()
