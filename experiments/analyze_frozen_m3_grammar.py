"""Four prespecified holdout test families; cluster representatives precede hits."""
from pathlib import Path
from itertools import combinations
import json,csv,hashlib,math
import numpy as np
from scipy.stats import fisher_exact,ks_2samp
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'data/processed/pr01_02_data_v2';OUT=ROOT/'tmp/m3_followup_20261008'
def read(p):return list(csv.DictReader((s for s in p.read_text().splitlines() if s and not s.startswith('#')),delimiter='\t'))
def table(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=sorted(set(k for r in rows for k in r)),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def choose_sites(rows):
 out={}
 for r in rows:
  if float(r['q-value'])>.05:continue
  key=(r['sequence_name'],r['motif_id']);order=lambda s:(float(s['q-value']),-float(s['score']),int(s['start']),int(s['stop']),s['strand'])
  if key not in out or order(r)<order(out[key]):out[key]=r
 return out
def gap(a,b):
 a,b=sorted([a,b],key=lambda s:(int(s['start']),int(s['stop']),s['strand']));return int(b['start'])-int(a['stop'])-1
def bh(ps):
 if not ps:return []
 ix=np.argsort(ps);adj=np.minimum.accumulate((np.array(ps)[ix]*len(ps)/np.arange(1,len(ps)+1))[::-1])[::-1];out=np.empty(len(ps));out[ix]=np.minimum(adj,1);return out.tolist()
def fisher(a,b,c,d):
 result=fisher_exact([[a,b],[c,d]],alternative='two-sided');o=result.statistic
 if a+c==0 or b+d==0:return {'raw_p':float(result.pvalue),'odds_ratio':'NA','OR_CI_low':'NA','OR_CI_high':'NA','CI_method':'unidentifiable_all_hit_or_all_absent'}
 cells=np.array([a,b,c,d],float);correction=bool((cells==0).any());cells+=.5 if correction else 0;a1,b1,c1,d1=cells;ratio=a1*d1/(b1*c1);se=np.sqrt(np.sum(1/cells));return {'raw_p':float(result.pvalue),'odds_ratio':o if np.isfinite(o) else 'inf','OR_CI_low':float(np.exp(np.log(ratio)-1.96*se)),'OR_CI_high':float(np.exp(np.log(ratio)+1.96*se)),'CI_method':'log_Wald_Haldane0.5' if correction else 'log_Wald'}
def permute_ks(x,y,seed=20261005):
 observed=ks_2samp(x,y).statistic;pool=np.array(x+y);rng=np.random.default_rng(seed);count=0
 for i in range(5000):
  rng.shuffle(pool);count+=ks_2samp(pool[:len(x)],pool[len(x):]).statistic>=observed-1e-12
 rng=np.random.default_rng(seed);d=[np.median(rng.choice(x,len(x),replace=True))-np.median(rng.choice(y,len(y),replace=True)) for _ in range(2000)];ci=np.quantile(d,[.025,.975]);return {'raw_p':(count+1)/5001,'KS_D':float(observed),'median_difference':float(np.median(x)-np.median(y)),'CI_low':float(ci[0]),'CI_high':float(ci[1])}
def main():
 config=json.loads((OUT/'config.json').read_text());manifest=json.loads((OUT/'manifest.json').read_text());assert manifest['status']=='completed';stats=[];description=[];exclusions=[];units=[]
 for species in sorted(OUT.glob('tjupan_*')):
  name=species.name;ids=[x['motif_id'] for x in config['main_motifs'] if x['species']==name];rows=[r for r in read(BASE/name/'main/sequence_master.tsv') if r['split']=='holdout'];clusters={}
  for r in rows:clusters.setdefault(r['leakage_group'],[]).append(r)
  pos=[];neg=[]
  for group,rs in sorted(clusters.items()):
   if len({r['label'] for r in rs})!=1:exclusions.append({'species':name,'cluster':group,'reason':'mixed_labels_in_cluster_confirmation_excluded','records':len(rs)});continue
   chosen=min(rs,key=lambda r:r['sequence_id']);(pos if chosen['label']=='1' else neg).append(chosen['sequence_id']);units.append({'species':name,'cluster':group,'representative':chosen['sequence_id'],'label':chosen['label'],'members':len(rs)})
  hits=read(species/'holdout_fimo/fimo.tsv');sites=choose_sites(hits)
  for mid in ids:
   a=sum((sid,mid) in sites for sid in pos);c=sum((sid,mid) in sites for sid in neg);row={'species':name,'family':'presence','motif1':mid,'motif2':'','positive_hits':a,'positive_clusters':len(pos),'negative_hits':c,'negative_clusters':len(neg),'status':'tested',**fisher(a,len(pos)-a,c,len(neg)-c)};stats.append(row)
   values={}
   for label,sids in [('positive',pos),('negative',neg)]:
    values[label]=[(int(sites[(sid,mid)]['start'])+int(sites[(sid,mid)]['stop']))/2-61 for sid in sids if (sid,mid) in sites]
    vs=values[label];description.append({'species':name,'kind':'position','motif1':mid,'motif2':'','label':label,'independent_hit_clusters':len(vs),'median':float(np.median(vs)) if vs else 'NA','min':min(vs) if vs else 'NA','max':max(vs) if vs else 'NA'})
   rs={'species':name,'family':'position','motif1':mid,'motif2':'','positive_hits':len(values['positive']),'negative_hits':len(values['negative'])}
   if min(map(len,values.values()))<20:rs.update(status='descriptive_only_insufficient_20_hit_clusters_per_arm',raw_p='NA',BH_q='NA')
   else:rs.update(status='tested',**permute_ks(values['positive'],values['negative']))
   stats.append(rs)
  for x,y in combinations(ids,2):
   pairsets={label:[sid for sid in sids if (sid,x) in sites and (sid,y) in sites] for label,sids in [('positive',pos),('negative',neg)]};a=len(pairsets['positive']);c=len(pairsets['negative']);stats.append({'species':name,'family':'cooccurrence','motif1':x,'motif2':y,'positive_hits':a,'positive_clusters':len(pos),'negative_hits':c,'negative_clusters':len(neg),'status':'tested',**fisher(a,len(pos)-a,c,len(neg)-c)})
   vals={label:[gap(sites[(sid,x)],sites[(sid,y)]) for sid in sids] for label,sids in pairsets.items()}
   for label,vs in vals.items():description.append({'species':name,'kind':'gap','motif1':x,'motif2':y,'label':label,'independent_hit_clusters':len(vs),'median':float(np.median(vs)) if vs else 'NA','min':min(vs) if vs else 'NA','max':max(vs) if vs else 'NA'})
   rs={'species':name,'family':'gap','motif1':x,'motif2':y,'positive_hits':a,'negative_hits':c}
   if min(a,c)<20:rs.update(status='descriptive_only_insufficient_20_double_hit_clusters_per_arm',raw_p='NA',BH_q='NA')
   else:rs.update(status='tested',**permute_ks(vals['positive'],vals['negative']))
   stats.append(rs)
 summary=[]
 for family in ['presence','cooccurrence','position','gap']:
  rows=[r for r in stats if r['family']==family];tested=[r for r in rows if r['status']=='tested']
  for r,q in zip(tested,bh([r['raw_p'] for r in tested])):r['BH_q']=q
  summary.append({'family':family,'planned_hypotheses':len(rows),'tested':len(tested),'BH_q_le_0_05':sum(r['BH_q']<=.05 for r in tested),'descriptive_only':len(rows)-len(tested)})
 table(OUT/'grammar_tests.tsv',stats);table(OUT/'grammar_family_summary.tsv',summary);table(OUT/'position_gap_descriptive.tsv',description);table(OUT/'confirmation_units.tsv',units);table(OUT/'confirmation_exclusions.tsv',exclusions)
 manifest['grammar_summary']=summary;manifest['confirmation_units']=len(units);manifest['mixed_label_clusters_excluded']=len(exclusions);manifest['output_sha256']={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json'};(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(summary,indent=2));print('mixed label clusters excluded:',len(exclusions))
if __name__=='__main__':main()
