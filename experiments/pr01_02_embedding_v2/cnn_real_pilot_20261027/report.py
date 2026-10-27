"""Generate course-facing figures and a summary from locked real-data candidates."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from .pipeline import HERE, ART, SUMMARY, STREME, read, read_meme, sha

def main():
    c=pd.read_csv(SUMMARY/'candidates.tsv',sep='\t',dtype={'seed':str});f=pd.read_csv(SUMMARY/'family.tsv',sep='\t',dtype={'seed':str})
    checks=read(SUMMARY/'verification.json');assert checks['status']=='PASSED'
    seeds=[20261005,20261006];refs=read_meme(STREME);fig,axes=plt.subplots(2,3,figsize=(14,6));colors=['#338a62','#397eae','#e8a444','#c95850']
    selected=[]
    for i,seed in enumerate(seeds):
        d=ART/f'seed_{seed}'/'SHIFT';pwms=read_meme(d/'motifs.meme');names=[x['candidate'] for x in read(d/'ranking.json')[:3]]
        for j,name in enumerate(names):
            ax=axes[i,j];pwm=pwms[name];bottom=np.zeros(len(pwm))
            for k,base in enumerate('ACGT'):
                ax.bar(np.arange(1,len(pwm)+1),pwm[:,k],bottom=bottom,color=colors[k],width=.85,label=base);bottom+=pwm[:,k]
            ax.set_ylim(0,1);ax.set_title(f'Seed {seed} / {name}');ax.set_xlabel('PWM position');ax.set_ylabel('Base probability')
            row=c[(c.method=='CNN_SHIFT')&(c.seed==str(seed))&(c.candidate==name)].iloc[0]
            ref=refs[row.nearest_streme_id];off=int(row.streme_offset);overlap=max(0,min(len(pwm),off+len(ref))-max(0,off))
            selected.append({'seed':seed,'candidate':name,'width':len(pwm),'reference_id':row.nearest_streme_id,'reference_width':len(ref),
                'overlap_bp':overlap,'offset':off,'orientation':row.streme_orientation,'cosine':row.streme_cosine,
                'positive_rate':row.dev_positive_rate,'natural_rate':row.dev_natural_rate,'shuffle_rate':row.dev_shuffle_rate})
    fig.legend(handles=axes[0,0].containers,labels=list('ACGT'),loc='upper center',bbox_to_anchor=(.5,.95),ncol=4,frameon=False)
    fig.suptitle('CNN candidates on real B. subtilis promoters',y=.995);fig.tight_layout(rect=(0,0,1,.89));fig.savefig(SUMMARY/'candidate_profiles.png',dpi=160,bbox_inches='tight');plt.close(fig)
    selected=pd.DataFrame(selected);selected.to_csv(SUMMARY/'selected_candidates.tsv',sep='\t',index=False,lineterminator='\n')
    lines=['# CNN 真实数据 Motif 候选小试','',
      '使用与传统 MEME/STREME 相同的 TJU B. subtilis 冻结数据：Discovery 正类 473、天然对照 559；Development 正类 67、天然对照 77。共训练2个 CNN，使用序列级正负标签；只用 Discovery 拟合和排序，Development 在过去已被观察，本轮属于探索。未使用 Holdout。','',
      '每个滤波器按 Discovery 天然对照激活的90%分位选择支持实例，至少10条，并比较原窗口与 Discovery 信息量选择的 ±2bp 统一平移。按 Discovery 命中差排序；各方法候选按天然对照90%分位阈值进行同一双链扫描。Development 天然对照比较按同源组计数：正类66组、天然负类77组；Shuffle 按正类同源组配对。不同候选的阈值按同一规则分别计算，并非 FIMO p/q。','',
      '## 前3个候选的序列组命中','','| 方法 | 种子 | 正类命中 | 天然负类命中 | 正类减天然 | 配对 Shuffle 命中 |','| --- | --- | ---: | ---: | ---: | ---: |']
    for row in f[f.view=='top3'].itertuples():
        lines.append(f'| {row.method} | {row.seed} | {row.positive_rate:.3f} | {row.natural_rate:.3f} | {row.natural_difference:+.3f} | {row.shuffle_rate:.3f} |')
    lines += ['', '对照为原有3个已筛选 STREME 主候选，候选来源与 CNN 不同；MEME 和全部 STREME 候选尚未纳入，因此本表是初步参考，不是完整方法排名。所有候选集合的结果另见[family.tsv](family.tsv)。','',
      '## CNN 平移候选展示','','![CNN候选碱基概率](candidate_profiles.png)','','上图为 PWM 每个位置的碱基概率，分别取两次运行 Discovery 排名前3；不是功能强度。下表比较既有同数据 STREME 主候选的局部相似，该参考并非独立来源。','',
      '| 种子 | 候选 | 最相似 STREME | 重叠长度 | 中心化余弦 | 正类命中 | 天然对照 | Shuffle |','| --- | --- | --- | ---: | ---: | ---: | ---: | ---: |']
    for r in selected.itertuples():
        lines.append(f'| {r.seed} | {r.candidate} | {r.reference_id} | {r.overlap_bp}/{r.reference_width} | {r.cosine:.3f} | {r.positive_rate:.3f} | {r.natural_rate:.3f} | {r.shuffle_rate:.3f} |')
    lines += ['', '两个种子各有1个保留滤波器的最佳平移不是0；前3候选的 RAW 与 SHIFT 命中相同，不能据此声称平移带来改善。更多候选和各运行结果见 [candidates.tsv](candidates.tsv)、[selected_candidates.tsv](selected_candidates.tsv) 和 [family.tsv](family.tsv)。',
      '', '本次候选是实际 CNN 输出，不等于验证的功能 motif；CNN 使用序列标签监督，打乱背景命中提示需要继续检查序列组成及候选特异性。下一步可先核对两种子候选的位置和实例，再按同一数据与口径补齐完整 MEME/STREME、One-hot 和 BERT 对照。','',
      '[拟合锁与来源](fit_lock.json)｜[评价哈希](score_manifest.json)｜[核验](verification.json)｜[方案](../README.md)。']
    (SUMMARY/'REPORT.md').write_text('\n'.join(lines)+'\n',encoding='utf8')
    save={'report_source_sha256':sha(Path(__file__)),'outputs_sha256':{x:sha(SUMMARY/x) for x in ['REPORT.md','candidate_profiles.png','selected_candidates.tsv']}}
    (SUMMARY/'report_manifest.json').write_text(json.dumps(save,indent=2)+'\n',encoding='utf8')
    print('Report generated:',len(selected),'profiles')

if __name__=='__main__':main()
