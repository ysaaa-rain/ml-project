"""Frozen reference comparisons and native FIMO q-values, discovery first."""
from pathlib import Path
import csv,json,hashlib,subprocess,os,re,collections
ROOT=Path(__file__).resolve().parents[1];TOOLS=ROOT/'tmp/meme-suite-5.5.9/bin';M3=ROOT/'results/motif/tju_streme_v2_20261008';BASE=ROOT/'data/processed/pr01_02_data_v2';OUT=ROOT/'tmp/m3_followup_20261008'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(path,rows):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['status'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def run(cmd,log):
 env=dict(os.environ);env['PATH']=str(TOOLS)+os.pathsep+env.get('PATH','')
 with log.open('w') as f:r=subprocess.run([str(x) for x in cmd],stdout=f,stderr=f,cwd=ROOT,env=env)
 assert r.returncode==0,(cmd,log)
def main():
 OUT.mkdir(exist_ok=False)
 selection=json.loads((M3/'manifest.json').read_text());v2=json.loads((BASE/'manifest.json').read_text());reference=ROOT/'tmp/reference_library_v1_20261008/reference_library_v1.meme'
 for f,h in v2['output_sha256'].items():assert sha(BASE/f)==h,f
 config={'main_candidates_sha256':sha(M3/'selected_main_motifs.tsv'),'main_motifs':selection['main_candidates'],'reference_sha256':sha(reference),'site_q':.05,'fimo_q_family':'per motif across all positions and strands in combined positive/natural-negative sequences of one species and one split','fimo_native_q_method':'native BH/Storey pi0 implementation, reservoir sampling per native defaults; no truncated-table recalculation','max_stored_scores':5000000,'motif_pseudocount':.1,'background':'zero-order, estimated once from discovery natural negatives; complement symmetric; total uniform pseudocount1','primary_site_sort':['q ascending','score descending','start ascending','stop ascending','strand lexical'],'cluster_unit':'holdout homogeneous-label clusters only; one lexicographically smallest sequence_id chosen per cluster BEFORE site inspection; mixed-label clusters excluded from confirmation and logged','tests':'Fisher presence/cooccurrence; two-sided KS with5000 whole-cluster permutations positions/gap; each distribution arm >=20 hit clusters','BH':'four separate families across all six species; NA not included; total hypothesis count retained','CI':'OR log-Wald with Haldane0.5 only for zero cells; cluster bootstrap 2000 for median differences','seed':20261005,'status':'frozen_before_any_project_holdout_scan'}
 (OUT/'config.json').write_text(json.dumps(config,indent=2)+'\n');jobs=[];comparisons=[];text=reference.read_text();parts=re.split(r'(?m)^MOTIF ',text)
 for species in sorted(M3.glob('tjupan_*')):
  name=species.name;folder=OUT/name;folder.mkdir();cohort=BASE/name/'main';motif=species/'selected_main.meme';assert sha(motif)==selection['output_sha256'][str(motif.relative_to(M3))]
  prefix='ecoli_' if name=='tjupan_escherichia_coli' else 'bsub_' if name=='tjupan_bacillus_subtilis' else None
  if prefix:
   ref=folder/'reference.meme';ref.write_text(parts[0]+''.join('MOTIF '+s for s in parts[1:] if s.split()[0].startswith(prefix)))
   cmd=[TOOLS/'tomtom','--o',folder/'known_tomtom','--dist','pearson','--min-overlap','5','--thresh','1',motif,ref];run(cmd,folder/'known_tomtom.log');comparisons.append({'species':name,'command':list(map(str,cmd)),'status':'completed_small_target_library_descriptive_similarity_only'})
  else:comparisons.append({'species':name,'status':'no_qualified_reference_all_candidates_retained_unannotated'})
  counts=collections.Counter(''.join(s.strip() for s in (cohort/'discovery.natural_control.fasta').read_text().splitlines() if not s.startswith('>')));n=sum(counts.values());freq={b:((counts[b]+counts[{'A':'T','T':'A','C':'G','G':'C'}[b]])/2+.25)/(n+1) for b in 'ACGT'};background=folder/'discovery_natural_zero_order.bg';background.write_text('# Frozen discovery natural control; strand-symmetric order0, total pseudocount1\n'+'\n'.join(f'{b} {freq[b]:.12f}' for b in 'ACGT')+'\n')
  for split in ['discovery','development','holdout']:
   p=cohort/f'{split}.positive.fasta';neg=cohort/f'{split}.natural_control.fasta';input=folder/f'{split}.combined.fasta';input.write_bytes(p.read_bytes()+neg.read_bytes());dest=folder/f'{split}_fimo';cmd=[TOOLS/'fimo','--o',dest,'--bfile',background,'--motif-pseudo','0.1','--max-stored-scores','5000000','--qv-thresh','--thresh','0.05',motif,input];jobs.append({'species':name,'split':split,'command':list(map(str,cmd)),'input':str(input),'input_sha256':sha(input),'motif_sha256':sha(motif),'background_sha256':sha(background),'positive_sha256':sha(p),'negative_sha256':sha(neg)})
 manifest={'config_sha256':sha(OUT/'config.json'),'tools':{t:subprocess.check_output([str(TOOLS/t),'--version'],text=True).strip() for t in ['fimo','tomtom']},'tool_sha256':{t:sha(TOOLS/t) for t in ['fimo','tomtom']},'comparisons':comparisons,'planned_fimo':jobs,'runs':[]};(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 # Complete all discovery/development and their review BEFORE external holdout.
 for splitset in [('discovery','development'),('holdout',)]:
  summaries=[]
  for j in jobs:
   if j['split'] not in splitset:continue
   folder=OUT/j['species'];log=folder/(j['split']+'_fimo.log');run(j['command'],log);content=log.read_text()
   assert not re.search(r'discard|truncat|max.*stored.*reached',content,re.I),(j,content)
   tsv=folder/(j['split']+'_fimo')/'fimo.tsv';rows=list(csv.DictReader((s for s in tsv.read_text().splitlines() if s and not s.startswith('#')),delimiter='\t'))
   assert all(float(r['q-value'])<=.05 for r in rows)
   j['strict_site_count']=len(rows);j['exit_code']=0;manifest['runs'].append(j);summaries.append({'species':j['species'],'split':j['split'],'q_le_0_05_sites':len(rows),'sequence_hit_count':len({r['sequence_name'] for r in rows})});print(j['species'],j['split'],len(rows),'strict sites',flush=True)
  if splitset[0]=='discovery':table(OUT/'discovery_development_preflight.tsv',summaries);manifest['preflight_passed_before_holdout']=True
  else:table(OUT/'holdout_scan_summary.tsv',summaries)
  (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 manifest['status']='completed';manifest['output_sha256']={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json'};(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
if __name__=='__main__':main()
