"""Aggregate expanded M3 without publishing sequences or per-record hits."""
from pathlib import Path
import csv,json,collections,math,sys
from scipy.stats import fisher_exact,binomtest
from .follow_expanded_m3 import ROOT,BASE,BIN,OUT as FOLLOW,readtable,readfa,table,execute,sha
OUT=ROOT/'results/motif/m3_expanded_20261009'
def bh(rows):
 groups=collections.defaultdict(list)
 for r in rows:
  if r['p']!='NA':groups[(r['split'],r['family'])].append(r)
 for group in groups.values():
  order=sorted(group,key=lambda r:r['p']);q=1.;n=len(order)
  for i in range(n-1,-1,-1):
   q=min(q,order[i]['p']*n/(i+1));order[i]['BH_q']=q;order[i]['family_size']=n
 return rows
def representatives(mapping,master):
 clusters=collections.defaultdict(list)
 for sid,label in mapping.items():clusters[master[sid]['leakage_group']].append((sid,label))
 reps={};mixed=0
 for group,items in clusters.items():
  if len({label for _,label in items})>1:mixed+=1;continue
  sid,label=min(items);reps[sid]=label
 return reps,mixed
def compare(mapping,master,hits,paired=False):
 reps,mixed=representatives(mapping,master)
 if paired:
  reps={sid:label for sid,label in reps.items() if label==1}
  a=sum(sid in hits for sid in reps);b=sum(sid+'__dinucleotide_null' in hits for sid in reps);n1=n0=len(reps)
  x=sum(sid in hits and sid+'__dinucleotide_null' not in hits for sid in reps);y=sum(sid not in hits and sid+'__dinucleotide_null' in hits for sid in reps);p=binomtest(x,x+y,.5).pvalue if x+y else 1.;test='paired_exact_binomial'
 else:
  n1=sum(label==1 for label in reps.values());n0=len(reps)-n1;a=sum(sid in hits and label==1 for sid,label in reps.items());b=sum(sid in hits and label==0 for sid,label in reps.items());p=float(fisher_exact([[a,n1-a],[b,n0-b]]).pvalue) if n1 and n0 else 'NA';x=y='NA';test='two_sided_fisher'
 return dict(positive_clusters=n1,control_clusters=n0,positive_hits=a,control_hits=b,rate_difference=a/n1-b/n0 if n1 and n0 else 'NA',p=p,BH_q='NA',mixed_clusters_excluded=mixed,discordant_positive_only=x,discordant_control_only=y,test=test,inference='exploratory_not_new_holdout_confirmation')
def main():
 m=json.loads((FOLLOW/'manifest.json').read_text());assert m['status']!='running' or '--partial' in sys.argv;native=json.loads((ROOT/'tmp/m3_expanded_20261009/manifest.json').read_text());OUT.mkdir(exist_ok=True);stats=[];sites=[];cases=[];similarity=[];method_runs=[];candidate=readtable(FOLLOW/'candidates.tsv')
 for case in m['cases']:
  cases.append({k:v for k,v in case.items() if k not in ['FIMO','comparisons']})
  if case.get('status') in ['discovery_failed_no_pwm','no_pwm_tool_found_no_motif']:continue
  key=case['key'];folder=FOLLOW/key;source=case['source'];group=case['group'];cohort=ROOT/case['cohort'];parent=BASE/source/('main' if source.startswith('tjupan') else 'core');master={r['sequence_id']:r for r in readtable(parent/'sequence_master.tsv')};mids=[r['motif_id'] for r in candidate if r['case']==key]
  for split in ['discovery','development','holdout']:
   p=readfa(cohort/f'{split}.positive.fasta');n=readfa(cohort/f'{split}.{case["control"]}.fasta');rows=readtable(folder/(split+'_strict')/'fimo.tsv');by=collections.defaultdict(set)
   for r in rows:by[r['motif_id']].add(r['sequence_name'])
   for mid in mids:
    if p and n:
     mapping={sid:1 for sid,_ in p}
     if case['control']=='natural_control':mapping.update({sid:0 for sid,_ in n})
     stats.append(dict(case=key,motif_id=mid,split=split,family='natural_presence' if case['control']=='natural_control' else 'paired_null_presence',**compare(mapping,master,by[mid],case['control']=='dinucleotide_null')))
    positions=[float(r['TSS_relative_center']) for r in readtable(folder/(split+'_strict_primary_sites.tsv')) if r['motif_id']==mid and r['sequence_name'] in {sid for sid,_ in p}]
    positions.sort();sites.append(dict(case=key,split=split,motif_id=mid,strict_site_count=sum(r['motif_id']==mid for r in rows),strict_positive_sequence_hits=sum(sid in by[mid] for sid,_ in p),strict_control_sequence_hits=sum(sid in by[mid] for sid,_ in n),positive_count=len(p),control_count=len(n),position_n=len(positions),TSS_center_min=min(positions) if positions else 'NA',TSS_center_median=(positions[(len(positions)-1)//2]+positions[len(positions)//2])/2 if positions else 'NA',TSS_center_max=max(positions) if positions else 'NA'))
   if group!='overall':
    pool=readtable(folder/f'{split}.sigma_pool_mapping.tsv');mapping={r['sequence_id']:int(r['target_label']) for r in pool};sigmahits=collections.defaultdict(set)
    for r in readtable(folder/(split+'_sigma_pool')/'fimo.tsv'):sigmahits[r['motif_id']].add(r['sequence_name'])
    for mid in mids:
     if mapping:stats.append(dict(case=key,motif_id=mid,split=split,family='sigma_specificity',**compare(mapping,master,sigmahits[mid])))
  # Record all applicable known-reference matches, no candidate deletion.
  for r in readtable(folder/'known_tomtom/tomtom.tsv'):similarity.append(dict(case=key,purpose='known_reference_descriptive',**r))
 # Tool-to-tool PWM comparison: each source/group/control separately.
 index={(c['source'],c['group'],c.get('control'),c['method']):c for c in m['cases'] if c.get('motif_count',0)}
 for (source,group,control,method),c in index.items():
  if method!='meme':continue
  other=index.get((source,group,control,'streme'))
  if not other:continue
  dest=FOLLOW/(c['key']+'__method_tomtom')
  if not dest.exists():
   r=execute([BIN/'tomtom','--o',dest,'--dist','pearson','--min-overlap','5','--thresh','1',FOLLOW/c['key']/'all_reported.meme',FOLLOW/other['key']/'all_reported.meme'],dest.with_suffix('.log'));assert r['exit_code']==0,r
   method_runs.append(r)
  for r in readtable(dest/'tomtom.tsv'):similarity.append(dict(case=c['key'],purpose='meme_vs_streme_exploratory',**r))
 table(OUT/'candidates.tsv',candidate);table(OUT/'cases.tsv',cases);table(OUT/'strict_FIMO_summary.tsv',sites);table(OUT/'cluster_tests.tsv',bh(stats));table(OUT/'PWM_similarity.tsv',similarity)
 evidence={'status':'partial_running' if m['status']=='running' else 'completed','native_status':native['status'],'versions':native['versions'],'tool_sha256':native['tool_sha256'],'config':json.loads((ROOT/'tmp/m3_expanded_20261009/config.json').read_text()),'jobs':native['jobs'],'reuse':native['reuse'],'skipped':native['skipped'],'followup_cases':m['cases'],'original17_hash_unchanged':sha(ROOT/'results/motif/tju_streme_v2_20261008/selected_main_motifs.tsv')==json.loads((ROOT/'tmp/m3_expanded_20261009/config.json').read_text())['original_selection_sha256'],'native_manifest_sha256':sha(ROOT/'tmp/m3_expanded_20261009/manifest.json'),'followup_manifest_sha256':sha(FOLLOW/'manifest.json'),'method_comparison_commands':method_runs,'statistics':'BH per split/family across all emitted candidates; independent pure-label cluster representatives; paired null exact binomial; every new holdout result exploratory','public_output_sha256':{p.name:sha(p) for p in OUT.glob('*.tsv')}}
 (OUT/'run_evidence.json').write_text(json.dumps(evidence,indent=2)+'\n')
 # Original FIMO distribution diagnostics are aggregates only.
 diagnostic=ROOT/'tmp/fimo_diagnostic_20261009';dm=json.loads((diagnostic/'manifest.json').read_text());table(OUT/'original17_diagnostic_distribution.tsv',readtable(diagnostic/'site_distribution_summary.tsv'));sea=[]
 for s in dm['sea']:
  for r in readtable(diagnostic/s['species']/(s['split']+'_sea/sea.tsv')):sea.append(dict(species=s['species'],split=s['split'],**{k:v for k,v in r.items() if k!='DB'}))
 table(OUT/'original17_SEA.tsv',sea);(OUT/'original17_diagnostic_evidence.json').write_text(json.dumps({k:v for k,v in dm.items() if k not in ['output_sha256','summary']},indent=2)+'\n')
 write_report(evidence,cases,sites,stats,candidate)
 print('cases',len(cases),'candidates',len(candidate),'tests',len(stats),'BH significant',sum(r['BH_q']!='NA' and r['BH_q']<=.05 for r in stats),flush=True)
def write_report(evidence,cases,sites,stats,candidate):
 partial=evidence['status']=='partial_running';jobs=evidence['jobs'];done=sum(r['status']=='completed' for r in jobs);failed=sum(r['status']=='failed' for r in jobs)
 text='# M3 扩展双方法与 σ 实验（2026-10-09）\n\n'
 text+=f"状态：{'**仍在运行，以下为阶段性结果**' if partial else '**本轮已完成（工具失败如实保留）**'}。新增发现完成 {done}/42，失败 {failed}；下游 case {len(cases)}/54（含已有 12 组 STREME）。\n\n"
 text+='原 natural STREME 的 17 主候选及原严格 holdout 结论不变。新增分析是探索性扩展；不据候选匹配或 holdout 结果重新选择。DBTBS 自身 σ 分析不等于独立外部验证。\n\n'
 if partial:text+='完整家族尚未汇齐，当前 BH 值仅阶段性计算，不能作为最终显著性结论。\n\n'
 text+='## 完整发现与扫描记录\n\n|来源/组|对照|方法|发现显著数|discovery 严格位点|development 严格位点|holdout 严格位点|状态|\n|---|---|---|---:|---:|---:|---:|---|\n'
 for c in cases:
  counts=[sum(int(r['strict_site_count']) for r in sites if r['case']==c['key'] and r['split']==split) for split in ['discovery','development','holdout']]
  text+='|'+c['source']+'/'+c['group']+'|'+c.get('control','NA')+'|'+c['method']+'|'+str(c.get('native_significant','NA'))+'|'+ '|'.join(map(str,counts))+'|'+c['status']+'|\n'
 text+='\n位点数量合计正样本与对应对照，并非独立启动子数；详细正/负序列命中数见 strict_FIMO_summary.tsv。发现显著不等于位点显著，也不等于 σ 特异。空分区登记跳过；未完成 case 未进入上表。\n\n'
 text+='## 样本不足与方法限制\n\n'
 text+='；'.join(f"{r['source']}/{r['group']} n={r['n']}，两方法跳过" for r in evidence['skipped'])+'。10–19 条组照预设执行，结果功效有限。STREME 内部评估正负总数各不足 5 的原生 p 不标为可靠发现显著。\n\n'
 text+='MEME DE/ZOOPS、零阶背景、内部评估比例 0.5；STREME 二阶背景、内部比例 0.1；宽度 5–15，最多 10 PWM，seed=20261005。内部拆分由工具实现，不保证簇级或原序列/shuffle 成对隔离，原生显著性只描述，下游以独立簇统计。天然对照 GC 差异和 shuffle 成对依赖仍是限制。\n\n'
 if not partial:
  text+='## 簇级统计汇总\n\n|检验家族|分区|有效检验数|BH q≤0.05|\n|---|---|---:|---:|\n'
  for family in sorted({r['family'] for r in stats}):
   for split in ['discovery','development','holdout']:
    rs=[r for r in stats if r['family']==family and r['split']==split and r['p']!='NA'];sig=sum(r['BH_q']!='NA' and r['BH_q']<=.05 for r in rs);text+=f'|{family}|{split}|{len(rs)}|{sig}|\n'
  text+='\n天然正负及实际 σ 对照用纯标签同源簇代表的双侧 Fisher；shuffle 用同一正序列与其 null 的成对精确二项检验。BH 按分区/比较目的跨所有输出 PWM 校正；混标签簇排除。完整效应量与计数见 cluster_tests.tsv；新增 holdout 均探索性，不能当新的确认性重复。位置只描述，不改变候选。\n\n'
 text+='## 证据与复现\n\n'
 text+='[首次 FIMO 诊断](M3_首次FIMO诊断_20261009.md)已完成：未发现明显输入/截断错误，1494 个 raw p<1e-4 位点与严格 site q 的稀疏结果并存。\n\n'
 text+='公共证据：results/motif/m3_expanded_20261009/{run_evidence.json,cases.tsv,candidates.tsv,strict_FIMO_summary.tsv,cluster_tests.tsv,PWM_similarity.tsv}。完整命令、版本、输入哈希与本地原生输出证据留在 manifest；原始序列、逐位点匹配序列、受限来源 PWM 均不公开。\n\n'
 text+='运行：`python -m experiments.run_expanded_m3` → `python -m experiments.follow_expanded_m3` → `python -m experiments.summarize_expanded_m3`。阶段性汇总使用 `--partial`，不能据其临时 BH 宣称最终显著性。测试：155 passed（2026-10-09）；原 17 selection SHA256 保持一致。\n'
 (ROOT/'reports/M3_扩展双方法与sigma实验_20261009.md').write_text(text)
if __name__=='__main__':main()
