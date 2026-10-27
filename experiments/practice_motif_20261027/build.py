"""Build a reproducible practice-oriented view of existing motif evidence."""
import html
import json
import platform
import subprocess
import time
from pathlib import Path

import numpy as np
import pandas as pd

from ..pr01_02_embedding_v2.budget_validation_20261026.pipeline import candidate_cache
from ..pr01_02_embedding_v2.discovery import ALPHABET,read_meme,similarity
from ..pr01_02_embedding_v2.refinement_20261026.pipeline import sequence_windows

HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
OUT=ROOT/'reports/practice_motif_20261027';BRANCH=ROOT/'experiments/pr01_02_embedding_v2'
PROTOCOL=HERE/'protocol.json'
SOURCE_TAGS={'tjupan_bacillus_subtilis':'B. subtilis','tjupan_baumannii':'A. baumannii',
    'tjupan_bradyrhizobium':'Bradyrhizobium','tjupan_diphtheria':'C. diphtheriae',
    'tjupan_escherichia_coli':'E. coli','tjupan_staphylococcus':'Staphylococcus'}


def read(path):return json.loads(path.read_text(encoding='utf8'))


def sha(path):
    import hashlib
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


def save(path,value):path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n',encoding='utf8')


def project_reference_start(starts,forward,offset):
    """For equal widths, project reference start after forward/reverse PWM scan."""
    return np.asarray(starts)+np.where(np.asarray(forward),offset,-offset)


def evaluate_partial(cache,matches,names,data,split,tolerance):
    family=np.zeros(len(data[split+'_control']),dtype=bool)
    for name in names:family|=cache[name][split]['blank']
    common={'split':split,'candidate_count':len(names),'family_blank_fpr':float(family.mean()),
            'truth_n':len(data[split]),'position_tolerance_bp':tolerance}
    winner=max(names,key=lambda name:(matches[name][0],name)) if names else None
    if winner is None:
        return {**common,'candidate':'','partial_cosine':0.,'offset':0,'overlap_bp':0,'orientation':'',
                'precision':0.,'recall':0.,'f1':0.,'tp':0,'prediction_count':0}
    sim,offset,orientation=matches[winner];scans=cache[winner][split]
    actual=np.array([r['start0'] for r in data[split+'_truth']])
    actual_forward=np.array([r['strand']=='+' for r in data[split+'_truth']])
    predicted_forward=scans['forward'] if orientation=='+' else ~scans['forward']
    predicted_start=project_reference_start(scans['start'],scans['forward'],offset)
    tp=int((scans['present']&(predicted_forward==actual_forward)&(np.abs(predicted_start-actual)<=tolerance)).sum())
    predictions=int(scans['present'].sum()+scans['blank'].sum())
    precision=tp/predictions if predictions else 0.;recall=tp/len(actual) if len(actual) else 0.
    return {**common,'candidate':winner,'partial_cosine':float(sim),'offset':int(offset),
        'overlap_bp':12-abs(int(offset)),
        'orientation':orientation,'precision':precision,'recall':recall,
        'f1':2*precision*recall/(precision+recall) if precision+recall else 0.,'tp':tp,'prediction_count':predictions}


def traditional(inputs):
    folder=ROOT/'results/motif/m3_expanded_20261009';manifest=read(folder/'publication_sha256.json')
    def load(name):
        p=folder/name;assert sha(p)==manifest[name];inputs[p.relative_to(ROOT).as_posix()]=sha(p)
        return pd.read_csv(p,sep='\t')
    candidates=load('candidates.tsv');cases=load('cases.tsv');scans=load('strict_FIMO_summary.tsv')
    stats=load('cluster_tests.tsv');similar=load('PWM_similarity.tsv')
    result=candidates.merge(cases[['key','source','group','control','status']],left_on='case',right_on='key',validate='many_to_one').drop(columns='key')
    selected_path=ROOT/'results/motif/tju_streme_v2_20261008/selected_main_motifs.tsv'
    inputs[selected_path.relative_to(ROOT).as_posix()]=sha(selected_path)
    selected=pd.read_csv(selected_path,sep='\t');original={(r.species,r.motif_id) for r in selected.itertuples()}
    result['original_main17']=[r.method=='streme' and r.control=='natural_control' and (r.source,r.native_motif_id) in original for r in result.itertuples()]
    assert result.original_main17.sum()==17 and len(result)==540
    for split in ['discovery','development']:
        rows=scans[scans.split==split][['motif_id','positive_count','control_count','strict_positive_sequence_hits','strict_control_sequence_hits','strict_site_count']].copy()
        rows=rows.rename(columns={c:f'{split}_{c}' for c in rows.columns if c!='motif_id'})
        result=result.merge(rows,on='motif_id',how='left',validate='one_to_one')
    for family in ['natural_presence','paired_null_presence','sigma_specificity']:
        rows=stats[(stats.split=='development')&(stats.family==family)][['motif_id','p','BH_q','rate_difference','positive_clusters','control_clusters']].copy()
        rows=rows.rename(columns={c:f'dev_{family}_{c}' for c in rows.columns if c!='motif_id'})
        result=result.merge(rows,on='motif_id',how='left',validate='one_to_one')
    ref=similar[similar.purpose=='known_reference_descriptive'].sort_values(['Query_ID','q-value','p-value','Target_ID'],na_position='last').drop_duplicates('Query_ID')
    ref=ref[['Query_ID','Target_ID','q-value','p-value','Optimal_offset','Overlap','Orientation']].rename(columns={'Query_ID':'motif_id','Target_ID':'reference_id',
        'q-value':'reference_q','p-value':'reference_p','Optimal_offset':'reference_offset','Overlap':'reference_overlap','Orientation':'reference_orientation'})
    result=result.merge(ref,on='motif_id',how='left',validate='one_to_one')
    return result,selected


def embedding(inputs):
    path=BRANCH/'summary/all_development_enrichment.tsv';inputs[path.relative_to(ROOT).as_posix()]=sha(path)
    frame=pd.read_csv(path,sep='\t');keys=['run','motif','method','source','layer','width','seed']
    assert not frame.duplicated(keys+['control']).any()
    columns=['positive_groups','control_groups','positive_hits','control_hits','odds_ratio','p','q','q_six_species','threshold']
    natural=frame[frame.control=='natural'][keys+columns].rename(columns={c:'natural_'+c for c in columns})
    shuffle=frame[frame.control=='shuffle'][keys+columns].rename(columns={c:'shuffle_'+c for c in columns})
    result=natural.merge(shuffle,on=keys,how='outer',validate='one_to_one');assert len(result)==1289
    result['natural_positive_rate']=result.natural_positive_hits/result.natural_positive_groups.replace(0,np.nan)
    result['natural_control_rate']=result.natural_control_hits/result.natural_control_groups.replace(0,np.nan)
    result['shuffle_positive_rate']=result.shuffle_positive_hits/result.shuffle_positive_groups.replace(0,np.nan)
    result['shuffle_control_rate']=result.shuffle_control_hits/result.shuffle_control_groups.replace(0,np.nan)
    metadata=[];manifests={}
    for run in result.run.unique():
        folder=BRANCH/'runs'/run;p=folder/'candidates.tsv';manifest_path=folder/'manifest.json'
        manifest=read(manifest_path);assert sha(p)==manifest['output_sha256']['candidates.tsv']
        manifests[run]=manifest
        inputs[p.relative_to(ROOT).as_posix()]=sha(p);inputs[manifest_path.relative_to(ROOT).as_posix()]=sha(manifest_path)
        f=pd.read_csv(p,sep='\t');f['run']=run;metadata.append(f)
    meta=pd.concat(metadata,ignore_index=True)[['run','motif','consensus','promoters','information_bits']]
    result=result.merge(meta,on=['run','motif'],how='left',validate='one_to_one')
    result['consensus']=result['consensus'].fillna(result.motif.str.replace(r'^STREME_[0-9]+-','',regex=True))
    def logo(r):
        p=BRANCH/'runs'/r.run/(str(r.motif)+'.png')
        if not p.is_file():return ''
        digest=sha(p);assert manifests[r.run]['output_sha256'][p.name]==digest
        inputs[p.relative_to(ROOT).as_posix()]=digest
        return p.relative_to(ROOT).as_posix()
    result['logo_path']=[logo(r) for r in result.itertuples()]
    return result


def synthetic(inputs,cfg):
    rows=[];view_count=0
    for version in cfg['synthetic_versions']:
        base=BRANCH/version;version_cfg=read(base/'protocol.json');lock=read(base/'summary/fit_lock.json')
        inputs[(base/'protocol.json').relative_to(ROOT).as_posix()]=sha(base/'protocol.json')
        inputs[(base/'summary/fit_lock.json').relative_to(ROOT).as_posix()]=sha(base/'summary/fit_lock.json')
        for source,digest in lock['source_sha256'].items():assert sha(ROOT/source)==digest,source
        def checked(path):
            rel=path.relative_to(base).as_posix();digest=sha(path)
            assert lock['outputs'].get(rel)==digest,rel
            inputs[path.relative_to(ROOT).as_posix()]=digest
            return path
        for folder in sorted((base/'artifacts').iterdir()):
            data=read(checked(folder/'data.json'));prepared={k:sequence_windows(data[k],version_cfg['width']) for k in
                ['calibration_control','development','development_control','test','test_control']}
            if version=='budget_validation_20261026':
                groups=[]
                for seed in version_cfg['algorithm_seeds']:
                    for variant in version_cfg['variants']:
                        view_file=checked(folder/f'views__a{seed}__{variant}.json')
                        groups.append((seed,variant,read(view_file),lambda method,s=seed,v=variant:folder/f'{method}__a{s}'/v/'motifs.meme'))
            else:
                groups=[]
                for parent in sorted(folder.glob('a*')):
                    seed=int(parent.name[1:]);view_file=checked(parent/'views.json')
                    groups.append((seed,'standard',read(view_file),lambda method,p=parent:p/method/'motifs.meme'))
            for seed,variant,views,pwm_path in groups:
                names={v['representation'] for v in views};caches={};matches={}
                reference=np.full((version_cfg['width'],4),data['scenario']['mutation']/3)
                for i,c in enumerate(version_cfg['synthetic_pattern']):reference[i,ALPHABET.index(c)]=1-data['scenario']['mutation']
                for method in names:
                    motifs=read_meme(checked(pwm_path(method)));caches[method],_=candidate_cache(motifs,data,version_cfg,prepared)
                    matches[method]={name:similarity(pwm,reference,minimum_overlap=cfg['partial_alignment']['min_overlap']) for name,pwm in motifs.items()}
                for view in views:
                    common={'version':version,'scenario':data['scenario']['id'],'data_seed':data['seed'],'algorithm_seed':seed,
                        'variant':variant,**{k:v for k,v in view.items() if k!='candidates'}}
                    for split in ['development','test']:
                        rows.append({**common,**evaluate_partial(caches[view['representation']],matches[view['representation']],view['candidates'],
                            data,split,cfg['partial_alignment']['position_tolerance_bp'])})
                    view_count+=1
            print(f'partial {version} {folder.name}: {view_count} views',flush=True)
    frame=pd.DataFrame(rows);assert view_count==cfg['expected']['synthetic_views'] and len(frame)==cfg['expected']['synthetic_split_rows']
    key=['version','scenario','variant','view','representation','split','data_seed']
    seeds=frame.groupby(key).agg(f1=('f1','mean'),partial_cosine=('partial_cosine','mean'),overlap_bp=('overlap_bp','mean'),
        family_blank_fpr=('family_blank_fpr','mean'),budget=('budget','mean'),zero_budget_runs=('zero_budget','sum')).reset_index()
    summary=seeds.groupby(key[:-1]).agg(f1_mean=('f1','mean'),f1_min=('f1','min'),f1_max=('f1','max'),
        partial_cosine=('partial_cosine','mean'),overlap_bp=('overlap_bp','mean'),family_blank_fpr=('family_blank_fpr','mean'),
        budget=('budget','mean'),zero_budget_runs=('zero_budget_runs','sum'),data_seeds=('data_seed','nunique')).reset_index()
    return frame,seeds,summary


def figures(selected,seeds,inputs):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    colors=['#3a8b62','#3579ac','#e5a13d','#c45d53']
    fig,axes=plt.subplots(2,3,figsize=(15,8));axes=axes.ravel()
    chosen=selected.sort_values(['species','rank']).groupby('species',sort=True).first().reset_index()
    assert len(chosen)==6
    for axis,row in zip(axes,chosen.itertuples()):
        p=ROOT/'results/motif/tju_streme_v2_20261008'/row.species/'selected_main.meme'
        inputs[p.relative_to(ROOT).as_posix()]=sha(p);pwm=read_meme(p)[row.motif_id];bottom=np.zeros(len(pwm))
        for j,base in enumerate(ALPHABET):
            axis.bar(np.arange(1,len(pwm)+1),pwm[:,j],bottom=bottom,color=colors[j],width=.88,label=base)
            bottom+=pwm[:,j]
        axis.set_title(SOURCE_TAGS[row.species]+'  |  '+row.motif_id,fontsize=10)
        axis.set(xlabel='Position in motif',ylabel='Base probability',ylim=(0,1));axis.set_xticks(range(1,len(pwm)+1,2))
    fig.legend(handles=axes[0].containers,labels=list(ALPHABET),loc='upper center',ncol=4,frameon=False,bbox_to_anchor=(.5,.96))
    fig.suptitle('Six existing discovery candidates: base-probability profiles',fontsize=13)
    fig.tight_layout(rect=(0,0,1,.91));fig.savefig(OUT/'six_species_profiles.png',dpi=160);plt.close(fig)
    recent=seeds[(seeds.version=='cnn_seed_validation_20261027')&(seeds.split=='test')&(seeds.view=='CAP5')]
    fig,axis=plt.subplots(figsize=(11,5));scenarios=['strong_random','mut30_random','mixed_strand_mut15']
    styles=[('TIED_RAW','#c26736',-.22),('TIED_SHIFT','#3a8b62',0),('OH48_BG','#3579ac',.22)]
    for name,color,offset in styles:
        for i,scenario in enumerate(scenarios):
            vals=recent[(recent.scenario==scenario)&(recent.representation==name)].f1.to_numpy()
            axis.scatter(i+offset+np.linspace(-.04,.04,len(vals)),vals,color=color,s=28,alpha=.85,label=name if i==0 else None)
            axis.plot([i+offset-.08,i+offset+.08],[vals.mean()]*2,color=color,lw=2.5)
    axis.set_xticks(range(3),['Strong','30% mutation','Mixed strand, 15%']);axis.set_ylim(-.03,1.04)
    axis.set_ylabel('Reference-assisted partial-pattern localization F1');axis.grid(axis='y',alpha=.2)
    axis.legend(loc='lower left',fontsize=9);axis.set_title('Synthetic exploration: overlap ≥8 bp, position tolerance ±3 bp',fontsize=11)
    fig.tight_layout();fig.savefig(OUT/'synthetic_partial_comparison.png',dpi=160);plt.close(fig)


def browser(traditional_frame,embedding_frame):
    # Inline aggregate metadata only; user input is rendered as text, not HTML.
    rows=[]
    def clean(v):
        if pd.isna(v):return None
        return v.item() if hasattr(v,'item') else v
    for r in traditional_frame.itertuples():
        rows.append({'kind':'传统方法','source':r.source,'method':r.method,'group':r.group,'background':r.control,
            'motif':r.motif_id,'consensus':r.consensus,'width':r.width,'positive_hits':r.development_strict_positive_sequence_hits,
            'positive_n':r.development_positive_count,'control_hits':r.development_strict_control_sequence_hits,
            'control_n':r.development_control_count,'natural_q':r.dev_natural_presence_BH_q,
            'shuffle_q':r.dev_paired_null_presence_BH_q,'reference':r.reference_id,'logo':''})
    for r in embedding_frame.itertuples():
        rows.append({'kind':'Embedding历史','source':r.source,'method':r.method,
            'group':f'L{r.layer} / W{r.width} / s{r.seed}','background':'natural + shuffle',
            'motif':r.motif,'consensus':r.consensus,'width':r.width,'positive_hits':r.natural_positive_hits,
            'positive_n':r.natural_positive_groups,'control_hits':r.natural_control_hits,'control_n':r.natural_control_groups,
            'natural_q':r.natural_q,'shuffle_q':r.shuffle_q,'reference':'','logo':r.logo_path})
    rows=[{k:clean(v) for k,v in row.items()} for row in rows]
    payload=json.dumps(rows,ensure_ascii=False,allow_nan=False).replace('</','<\\/')
    template='''<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Motif 候选浏览</title><style>
body{font-family:"Microsoft YaHei",system-ui,sans-serif;margin:0;color:#20314a;background:#f4f7fa}main{max-width:1480px;margin:auto;padding:28px}
h1{margin:0 0 8px;font-size:29px}p{line-height:1.6;color:#4d5e73}.panel{background:white;border:1px solid #dce4ed;border-radius:12px;padding:18px;margin:16px 0}
.filters{display:flex;gap:12px;flex-wrap:wrap;align-items:end}label{font-size:13px;color:#52647b;display:flex;flex-direction:column;gap:5px}
input,select,button{font:inherit;border:1px solid #bccbd9;border-radius:6px;padding:8px;background:white}input{width:260px}button{cursor:pointer;color:#155a82}button:disabled{opacity:.4}
.scroll{overflow:auto}table{border-collapse:collapse;width:100%;font-size:13px;white-space:nowrap}th,td{padding:9px;border-bottom:1px solid #e9eff5;text-align:left}th{background:#edf4f9;position:sticky;top:0}td.motif{max-width:310px;overflow:hidden;text-overflow:ellipsis}
.muted{font-size:13px}.pages{display:flex;gap:10px;align-items:center;margin-top:14px}a{color:#17678d}code{font-size:12px}</style></head><body><main>
<h1>Motif 候选浏览</h1><p>浏览已有真实数据运行记录，按来源、方法或关键词筛选。每一行是一条运行级候选，跨运行重复不算独立发现。</p>
<div class="panel"><div class="filters"><label>结果来源<select id="kind"><option value="">全部</option><option>传统方法</option><option>Embedding历史</option></select></label>
<label>方法<select id="method"><option value="">全部</option></select></label><label>搜索来源、序列或候选名<input id="search" type="search" placeholder="例如 ecoli、TATAAT"></label><label>每页<select id="size"><option>25</option><option selected>50</option><option>100</option></select></label></div>
<p class="muted">正/对照命中均取既有 Development 汇总；传统方法和 Embedding 的扫描、q 值计算口径不同，不宜直接比较。参考名称仅表示相似对照。表格可左右滚动。</p>
<div class="scroll"><table><thead><tr><th>来源</th><th>方法</th><th>组/背景</th><th>候选</th><th>共识</th><th>宽度</th><th>正类命中</th><th>对照命中</th><th>天然 q</th><th>shuffle q</th><th>参考/图</th></tr></thead><tbody id="body"></tbody></table></div>
<div class="pages"><button id="prev">上一页</button><span id="count" aria-live="polite"></span><button id="next">下一页</button></div></div>
<p class="muted">完整字段与来源见 <a href="REPORT.md">汇总说明</a>、<a href="traditional_candidates.tsv">传统候选表</a> 和 <a href="embedding_candidates.tsv">Embedding候选表</a>。</p>
</main><script>const rows=__DATA__;const $=id=>document.getElementById(id);let page=0;
const esc=x=>String(x??'—').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const fmt=x=>x==null?'—':typeof x==='number'?(x>0&&x<.001?x.toExponential(2):Number(x.toPrecision(4)).toString()):String(x);
const methods=[...new Set(rows.map(x=>x.method))].sort();for(const m of methods){const o=document.createElement('option');o.value=m;o.textContent=m;$('method').append(o)}
function filtered(){const q=$('search').value.trim().toLowerCase();return rows.filter(r=>(!$('kind').value||r.kind===$('kind').value)&&(!$('method').value||r.method===$('method').value)&&(!q||Object.values(r).join(' ').toLowerCase().includes(q)))}
function draw(){const list=filtered(),size=Number($('size').value),last=Math.max(0,Math.ceil(list.length/size)-1);page=Math.min(page,last);const start=page*size;
$('body').innerHTML=list.slice(start,start+size).map(r=>{const detail=r.logo?'<a href="../../'+esc(r.logo)+'">Logo</a>':esc(r.reference);
return '<tr><td>'+esc(r.source)+'</td><td>'+esc(r.method)+'</td><td>'+esc(r.group)+' / '+esc(r.background)+'</td><td class="motif" title="'+esc(r.motif)+'">'+esc(r.motif)+'</td><td>'+esc(r.consensus)+'</td><td>'+esc(r.width)+'</td><td>'+esc(fmt(r.positive_hits))+'/'+esc(fmt(r.positive_n))+'</td><td>'+esc(fmt(r.control_hits))+'/'+esc(fmt(r.control_n))+'</td><td>'+esc(fmt(r.natural_q))+'</td><td>'+esc(fmt(r.shuffle_q))+'</td><td>'+detail+'</td></tr>'}).join('');
$('count').textContent=list.length?'第 '+(start+1)+'–'+Math.min(start+size,list.length)+' 条，共 '+list.length+' 条':'没有匹配结果';$('prev').disabled=page===0;$('next').disabled=page===last}
for(const id of ['kind','method','search','size'])$(id).addEventListener('input',()=>{page=0;draw()});$('prev').onclick=()=>{page--;draw()};$('next').onclick=()=>{page++;draw()};draw();</script></body></html>'''
    (OUT/'index.html').write_text(template.replace('__DATA__',payload),encoding='utf8')


def table(frame):
    f=frame.copy()
    for key in f:
        if f[key].dtype.kind=='f':f[key]=f[key].map(lambda x:'—' if pd.isna(x) else f'{x:.3f}')
    return '\n'.join(['| '+' | '.join(map(str,f.columns))+' |','| '+' | '.join(['---']*len(f.columns))+' |',
        *['| '+' | '.join(map(str,row))+' |' for row in f.itertuples(index=False,name=None)]])


def report(trad,emb,summary,cfg):
    latest=summary[(summary.version=='cnn_seed_validation_20261027')&(summary.variant=='standard')&(summary.view=='CAP5')&(summary.split=='test')]
    latest=latest[['scenario','representation','f1_mean','partial_cosine','overlap_bp','family_blank_fpr','budget','data_seeds']]
    counts=trad.groupby(['source','method']).size().reset_index(name='运行级候选数')
    methods=emb.groupby('method').size().reset_index(name='运行级记录数')
    text=f'''# 启动子 Motif 课程实践：候选与结果总览

整理日期：2026-10-27。当前目录把已有真实数据结果与合成方法实验集中展示，便于从来源、候选形状、对照和位置几个角度阅读。

## 浏览入口

- [交互式候选浏览](index.html)：按来源、方法或关键词筛选，查看 Development 正类/对照命中、既有 q 值与已保存 Logo。
- [传统方法候选表](traditional_candidates.tsv)：{len(trad)}条 MEME/STREME 运行级候选，含来源、组、背景、发现结果、Development命中及参考相似记录。
- [Embedding历史候选表](embedding_candidates.tsv)：{len(emb)}条运行级候选/对照记录，分别保留天然和Shuffle结果；STREME对照在不同运行中重复出现，不算新的独立motif。
- [合成片段评价汇总](synthetic_partial_summary.tsv)｜[逐数据种子](synthetic_partial_per_data_seed.tsv)｜[全部评价记录](synthetic_partial_metrics.tsv)。

候选数是运行记录数，不是去重后的生物学功能motif数。两条真实数据路线的扫描方法、对照、q值家族不同，不把它们的q值当作相同指标。表中空值表示没有可用记录，不填零；已知参考相似仅作注释。

## 真实数据：候选分布

传统发现中实际记录540条PWM，覆盖TJU六物种、RegulonDB和DBTBS各可运行组、MEME与STREME；部分小样本组或辅助运行有既存限制，见[扩展结果原报告](../M3_扩展双方法与sigma实验_20261009.md)。原有17个自然对照主候选在表中有标记；这里另展示全部运行级候选。DBTBS与TJU有高度重叠，不计作独立外部验证。

{table(counts)}

![六物种候选碱基概率](six_species_profiles.png)

每张图使用各物种Discovery阶段既有主候选中的第1个PWM；柱高是每个位置的A/C/G/T概率，并非功能强度。每个候选、完整PWM与原运行依据见传统方法表和历史结果。

## 真实数据：Embedding及对照

{table(methods)}

53个历史运行提供天然对照和配对Shuffle的Development聚合记录；每条记录展示正类/对照命中数、原始p、各自q及运行参数。天然对照和Shuffle含义不同，分别阅读；重复的STREME对照记录并非新发现。详情见[Embedding原报告](../../experiments/pr01_02_embedding_v2/REPORT.md)和[事实汇总](../../experiments/pr01_02_embedding_v2/facts_20261010/实验事实汇总_GitHub.md)。

## 合成实验：局部模式与位置

基于四轮已保存候选和锁定的预算集合，重新计算局部模式评价，共{cfg['expected']['synthetic_views']}个候选集合、Development/Test各一行。已知人工参考只用于评价：在每个已冻结集合内寻找至少8/12bp重叠的最佳参考相似，考虑正反向和偏移，将参考起点投影后按±3bp统计定位F1；不按相似度额外筛掉候选。表中同时给出平均相似度、重叠长度、候选数量及整个候选集的空白误报；空集合保留0分。

这是一项使用已知参考和已看过的数据的事后探索：低相似度也有记录，F1不能单独视作motif真实性、真实数据性能或新的独立验证。各数据种子先平均同一数据的算法运行，再汇总数据种子；CNN使用序列级监督，One-hot聚类没有相同分类训练。

### 最新CNN合成轮次：CAP5

{table(latest)}

![合成局部模式定位](synthetic_partial_comparison.png)

每点表示1组数据中各算法种子的平均，横线为数据均值。更多轮次、预算和Development/Test结果可在TSV中按字段筛选。

## 使用说明

1. 先看候选的来源、背景和序列形状，再结合命中数、对照和已保存的统计值解释。
2. 真实序列候选、人工植入模式和参考相似度分别理解；不能把合成参考辅助定位结果写成真实启动子发现。
3. 本次没有重训模型、重跑FIMO或调整真实留出数据；如果后续根据这份探索结果修改算法，应另设运行和验证数据。既有科学结果、负面观察及工具失败记录仍可从原报告追溯。

所有输入和输出SHA256、处理配置及核验计数见 [manifest.json](manifest.json)。
'''
    (OUT/'REPORT.md').write_text(text,encoding='utf8')


def main():
    cfg=read(PROTOCOL);OUT.mkdir(exist_ok=False,parents=True);began=time.time();inputs={}
    inputs[PROTOCOL.relative_to(ROOT).as_posix()]=sha(PROTOCOL)
    for source in [HERE/'build.py',HERE/'__init__.py',HERE/'test_overview.py']:
        inputs[source.relative_to(ROOT).as_posix()]=sha(source)
    trad,selected=traditional(inputs);emb=embedding(inputs);frame,seeds,summary=synthetic(inputs,cfg)
    trad.to_csv(OUT/'traditional_candidates.tsv',sep='\t',index=False)
    emb.to_csv(OUT/'embedding_candidates.tsv',sep='\t',index=False)
    frame.to_csv(OUT/'synthetic_partial_metrics.tsv',sep='\t',index=False)
    seeds.to_csv(OUT/'synthetic_partial_per_data_seed.tsv',sep='\t',index=False)
    summary.to_csv(OUT/'synthetic_partial_summary.tsv',sep='\t',index=False)
    figures(selected,seeds,inputs);browser(trad,emb);report(trad,emb,summary,cfg)
    outputs={p.name:sha(p) for p in OUT.iterdir() if p.is_file()}
    manifest={'status':'COMPLETED','version':cfg['version'],'started_from_git':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
        'python':platform.python_version(),'traditional_candidate_rows':len(trad),'embedding_candidate_rows':len(emb),
        'synthetic_view_sets':len(frame)//2,'synthetic_split_rows':len(frame),'source_hashes':inputs,'output_sha256':outputs,
        'analysis_status':'post-hoc exploratory, no retraining, no new real Holdout access','elapsed_seconds':time.time()-began}
    save(OUT/'manifest.json',manifest)
    print(json.dumps({k:v for k,v in manifest.items() if k not in ['source_hashes','output_sha256']},ensure_ascii=False),flush=True)


if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--refresh-browser',action='store_true')
    if parser.parse_args().refresh_browser:
        assert (OUT/'manifest.json').is_file()
        browser(pd.read_csv(OUT/'traditional_candidates.tsv',sep='\t'),pd.read_csv(OUT/'embedding_candidates.tsv',sep='\t'))
        manifest=read(OUT/'manifest.json')
        manifest['source_hashes'][(HERE/'build.py').relative_to(ROOT).as_posix()]=sha(HERE/'build.py')
        manifest['output_sha256']['index.html']=sha(OUT/'index.html')
        manifest['browser_revision']='Updated only the browser labels; scientific tables and figures are unchanged.'
        save(OUT/'manifest.json',manifest);print('Candidate browser refreshed')
    else:main()
