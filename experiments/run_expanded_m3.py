"""User-authorized 2026-10-09 matrix; immutable original17, all new PWM kept."""
from pathlib import Path
import concurrent.futures as cf, csv,json,hashlib,subprocess,os,time,collections
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'data/processed/pr01_02_data_v2';BIN=ROOT/'tmp/meme-suite-5.5.9/bin';OUT=ROOT/'tmp/m3_expanded_20261009'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def fa(p):return p.read_text().count('>') if p.exists() else 0
def background(negative,path,order):
 # MEME and FIMO use explicit order0 discovery-control composition;
 # STREME estimates its own order2 background, including auxiliary branches.
 c=collections.Counter(''.join(s.strip() for s in negative.read_text().splitlines() if not s.startswith('>')));n=sum(c.values());freq={b:((c[b]+c[{'A':'T','T':'A','C':'G','G':'C'}[b]])/2+.25)/(n+1) for b in 'ACGT'};path.write_text('# Discovery control only; complement-symmetric zero-order\n'+'\n'.join(f'{b} {freq[b]:.12f}' for b in 'ACGT')+'\n')
def main():
 OUT.mkdir(exist_ok=False)
 frozen=json.loads((BASE/'manifest.json').read_text())
 for f,h in frozen['output_sha256'].items():assert sha(BASE/f)==h,f
 jobs=[];skips=[];reuse=[]
 def add(cohort,source,group,control,method):
  p=cohort/'discovery.positive.fasta';n=cohort/f'discovery.{control}.fasta';key=f'{source}__{group}__{control}__{method}';directory=OUT/key;directory.mkdir();b=directory/'discovery_control_zero_order.bg';background(n,b,0)
  if method=='meme':cmd=[str(BIN/'meme'),str(p),'-dna','-revcomp','-mod','zoops','-objfun','de','-test','mhg','-neg',str(n),'-bfile',str(b),'-markov_order','0','-hsfrac','0.5','-searchsize','100000','-nmotifs','10','-minw','5','-maxw','15','-seed','20261005','-maxsize','10000000','-o',str(directory/'discovery')]
  else:cmd=[str(BIN/'streme'),'--dna','--p',str(p),'--n',str(n),'--order','2','--minw','5','--maxw','15','--thresh','0.05','--nmotifs','10','--seed','20261005','--o',str(directory/'discovery')]
  # STREME remains order2, as the original user requirement. MEME DE uses
  # explicit order0; this is documented method difference, not hidden tuning.
  jobs.append({'key':key,'source':source,'group':group,'method':method,'control':control,'cohort':str(cohort.relative_to(ROOT)),'positive':str(p.relative_to(ROOT)),'negative':str(n.relative_to(ROOT)),'positive_count':fa(p),'negative_count':fa(n),'positive_sha256':sha(p),'negative_sha256':sha(n),'background_sha256':sha(b),'command':cmd,'output':str((directory/'discovery').relative_to(ROOT)),'status':'planned','inference':'new exploratory analysis; original holdout already inspected'})
 for cohort in sorted(BASE.glob('tjupan_*/main')):
  for control in ['natural_control','dinucleotide_null']:
   add(cohort,cohort.parent.name,'overall',control,'meme');reuse.append({'source':cohort.parent.name,'group':'overall','control':control,'method':'streme','native_output':'results/motif/tju_streme_v2_20261008/'+cohort.parent.name+('/natural' if control=='natural_control' else '/dinucleotide')})
 for source in ['regulondb_ecoli','dbtbs_bsub']:
  core=BASE/source/'core'
  for method in ['meme','streme']:add(core,source,'overall','dinucleotide_null',method)
  for group in sorted((core/'sigma').iterdir()):
   if not group.is_dir():continue
   count=fa(group/'discovery.positive.fasta')
   if count<10:skips.append({'source':source,'group':group.name,'n':count,'methods':['meme','streme'],'reason':'discovery_n_lt_10; no formal discovery; not fabricated significance'});continue
   for method in ['meme','streme']:add(group,source,group.name,'dinucleotide_null',method)
 config={'date':'2026-10-09','width':[5,15],'nmotifs':10,'seed':20261005,'meme_objfun':'de','meme_model':'zoops','meme_hsfrac':.5,'meme_searchsize':100000,'meme_background_order':0,'streme_background_order':2,'streme_hofract':.1,'significance':'MEME native E<=.05; STREME internal-test p<=.05 only if test positive and negative counts>=5; all reported candidates go to FIMO','q_vs_p':'strict q<=.05; diagnostic p<1e-4 discovery/development only','small_n':'10-19 attempt both; failure/low power kept; <10 skip both','original17_unchanged':True,'comparison_limits':'internal splits/background/model differ by tool; newly inspected holdout is descriptive/exploratory, never model selection','original_selection_sha256':sha(ROOT/'results/motif/tju_streme_v2_20261008/selected_main_motifs.tsv')};(OUT/'config.json').write_text(json.dumps(config,indent=2)+'\n')
 manifest={'status':'running','git_commit':subprocess.check_output(['git','rev-parse','HEAD'],text=True,cwd=ROOT).strip(),'config_sha256':sha(OUT/'config.json'),'versions':{t:subprocess.check_output([str(BIN/t),'-version' if t=='meme' else '--version'],text=True).strip() for t in ['meme','streme','fimo','tomtom','sea']},'tool_sha256':{t:sha(BIN/t) for t in ['meme','streme','fimo','tomtom','sea']},'jobs':jobs,'reuse':reuse,'skipped':skips}
 def save():(OUT/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 save();env=dict(os.environ);env['PATH']=str(BIN)+os.pathsep+env.get('PATH','')
 def worker(j):
  start=time.monotonic();folder=OUT/j['key']
  with (folder/'stdout.log').open('w') as stdout,(folder/'stderr.log').open('w') as stderr:r=subprocess.run(j['command'],cwd=ROOT,env=env,stdout=stdout,stderr=stderr)
  j.update(status='completed' if r.returncode==0 else 'failed',exit_code=r.returncode,elapsed_seconds=round(time.monotonic()-start,2));print(j['key'],j['status'],j['elapsed_seconds'],'s',flush=True);return j
 # Independent native processes; no timeout or significance-dependent retry.
 with cf.ThreadPoolExecutor(max_workers=3) as pool:
  fs=[pool.submit(worker,j) for j in jobs]
  for f in cf.as_completed(fs):f.result();save()
 manifest['status']='completed_with_recorded_tool_failures' if any(j['status']=='failed' for j in jobs) else 'completed';manifest['output_sha256']={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json'};save()
if __name__=='__main__':main()
