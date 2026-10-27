"""Produce a readable comparison from the scored same-batch output."""
import hashlib
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'reports/same_batch_motif_equal_20261027'
METHODS=['EMB128','OH24','OH48','TIED_RAW','TIED_SHIFT']
LABELS={'EMB128':'BERT','OH24':'One-hot PCA24','OH48':'One-hot48','TIED_RAW':'CNN raw','TIED_SHIFT':'CNN shift'}

def main():
    per=pd.read_csv(OUT/'per_data_seed.tsv',sep='\t')
    run=pd.read_csv(OUT/'runs.tsv',sep='\t')
    assert len(run)==450 and len(per)==150 and not run.duplicated(['scenario','data_seed','algorithm_seed','method','split']).any()
    assert run.f1.between(0,1).all() and run.family_blank_fpr.between(0,1).all()
    summary=pd.read_csv(OUT/'summary.tsv',sep='\t');test=summary[summary['split']=='test']
    assert (test.n_data==5).all() and len(test)==15
    scenarios=['strong_random','mut30_random','mixed_strand_mut15']
    fig,ax=plt.subplots(figsize=(11,5));w=.15
    for i,method in enumerate(METHODS):
        vals=[test[(test.scenario==s)&(test.method==method)].f1_mean.iloc[0] for s in scenarios]
        ax.bar(np.arange(3)+(i-2)*w,vals,width=w,label=LABELS[method])
    ax.set_xticks(range(3),['Strong pattern','30% mutation','Mixed strand, 15% mutation'])
    ax.set_ylim(0,1);ax.set_ylabel('Reference-assisted partial localization F1')
    ax.set_title('Identical synthetic data, candidate budget and scoring rule')
    ax.grid(axis='y',alpha=.2);ax.legend(ncol=3,fontsize=8);fig.tight_layout()
    fig.savefig(OUT/'comparison.png',dpi=160);plt.close(fig)
    lines=['# 同一批合成数据：BERT、One-hot 与 CNN 的比较','',
      '本轮使用原BERT实验的同一批输入：3种场景×5个数据种子，每个数据种子下3个算法种子。BERT、One-hot沿用已保存的Discovery候选；CNN在同一批Discovery正类和空白序列上新训练，共45个模型，原窗口和平移共用模型。没有使用真实数据或重新训练BERT。','',
      '统一按Discovery正类与空白命中率差排序，每组数据和算法种子取各方法共同可用的候选数量，最多5个。合成评价使用已知参考，在候选集合内选最佳局部相似项：重叠至少8/12 bp，正反向及偏移均考虑，位置容差±3 bp；空白序列预测计入精确率，同时单列候选集合空白误报。算法种子先在每组数据内求平均，再汇总5组数据。','',
      '该批BERT和One-hot的Test结果在此前已观察，本轮CNN是在旧数据上追加的探索性同批对比；参考辅助选择候选，不等于参考无关的自动发现。CNN有序列级标签监督，BERT及One-hot聚类没有同样的分类训练。','',
      '![同批合成数据比较](comparison.png)','','| 场景 | 方法 | 平均F1 | 数据种子范围 | 空白误报 |','| --- | --- | ---: | ---: | ---: |']
    for scenario in scenarios:
        for method in METHODS:
            row=test[(test.scenario==scenario)&(test.method==method)].iloc[0]
            lines.append(f'| {scenario} | {LABELS[method]} | {row.f1_mean:.3f} | {row.f1_min:.3f}–{row.f1_max:.3f} | {row.blank_fpr:.3f} |')
    lines += ['', '逐数据种子配对差值（CNN平移−BERT）用于显示数据间变化：','', '| 场景 | 平均差值 | 正/负/平 |', '| --- | ---: | --- |']
    for scenario in scenarios:
        d=per[(per.scenario==scenario)&(per['split']=='test')].pivot(index='data_seed',columns='method',values='f1')
        diff=d.TIED_SHIFT-d.EMB128
        lines.append(f'| {scenario} | {diff.mean():+.3f} | {(diff>1e-12).sum()}/{(diff< -1e-12).sum()}/{(abs(diff)<=1e-12).sum()} |')
    lines += ['', '[逐次运行](runs.tsv)｜[逐数据种子](per_data_seed.tsv)｜[汇总](summary.tsv)｜[来源与校验](manifest.json)。候选不足5个及低分均保留。',
      '', '传统MEME/STREME没有在这批合成数据上运行，本机当前未找到可用的MEME Suite，因此没有填入或借用其他数据的分数。真实数据的统一扫描另待完成。']
    (OUT/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    manifest=json.loads((OUT/'manifest.json').read_text(encoding='utf8'))
    manifest['outputs_sha256'].update({p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [OUT/'REPORT.md',OUT/'comparison.png']})
    manifest['report_source_sha256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(test.to_string(index=False))

if __name__=='__main__':main()
