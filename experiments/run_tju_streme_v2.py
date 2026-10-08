"""Frozen discovery-only STREME runs; no external holdout is read."""
from pathlib import Path
import concurrent.futures as cf
import hashlib,json,subprocess,time,datetime,os,shutil
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 base=ROOT/'data/processed/pr01_02_data_v2';frozen=json.loads((base/'manifest.json').read_text())
 for f,h in frozen['raw_sha256'].items():assert sha(ROOT/f)==h,f
 for f,h in frozen['output_sha256'].items():assert sha(base/f)==h,f
 out=ROOT/'results/motif/tju_streme_v2_20261008';out.mkdir(exist_ok=False,parents=True)
 binary=ROOT/'tmp/meme-suite-5.5.9/bin/streme';version=subprocess.check_output([str(binary),'--version'],text=True).strip()
 protocol=ROOT/'docs/PR01-02_motif_reference_and_grammar_protocol_20261008.md';shutil.copy2(protocol,out/'protocol_frozen.md');shutil.copy2(ROOT/'tmp/tju_positive_subset_verification.json',out/'positive_subset_verification.json')
 selection={'significance':'STREME internal discovery holdout p <= 0.05 (not FIMO q or project holdout)','sort':'p ascending, numeric motif id ascending','max_per_species':5,'duplicate_rule':'exact PWM matrix or exact reverse-complement PWM only; no similarity-family collapsing','background':'natural only; robustness never changes selection','order':2,'internal_hofract':0.1,'nmotifs_overrides_threshold_stopping':True}
 manifest={'status':'running','date':'2026-10-08','git_commit':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),'tool_version':version,'tool_sha256':sha(binary),'input_manifest_sha256':sha(base/'manifest.json'),'protocol_sha256':sha(protocol),'selection_rule':selection,'external_holdout_used':False,'runs':[]}
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 jobs=[]
 for cohort in sorted(base.glob('tjupan_*/main')):
  for name,control in [('natural','natural_control'),('dinucleotide','dinucleotide_null')]:
   p=cohort/'discovery.positive.fasta';n=cohort/f'discovery.{control}.fasta'
   for f in [p,n]:
    lines=f.read_text().splitlines();assert all(len(s)==81 and set(s)<=set('ACGT') for s in lines if not s.startswith('>'))
   dest=out/cohort.parent.name/name;dest.parent.mkdir(exist_ok=True)
   cmd=[str(binary),'--dna','--p',str(p),'--n',str(n),'--order','2','--minw','5','--maxw','15','--thresh','0.05','--nmotifs','10','--seed','20261005','--o',str(dest)]
   jobs.append({'species':cohort.parent.name,'background':name,'positive':str(p.relative_to(ROOT)),'negative':str(n.relative_to(ROOT)),'positive_sha256':sha(p),'negative_sha256':sha(n),'positive_count':p.read_text().count('>'),'negative_count':n.read_text().count('>'),'command':cmd,'output':str(dest.relative_to(ROOT))})
 manifest['planned_runs']=jobs;(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 def run(job):
  dest=ROOT/job['output'];start=time.monotonic();env=dict(os.environ);env['PATH']=str(binary.parent)+os.pathsep+env.get('PATH','')
  with (dest.parent/(job['background']+'.stdout.log')).open('w') as stdout,(dest.parent/(job['background']+'.stderr.log')).open('w') as stderr:
   result=subprocess.run(job['command'],cwd=ROOT,env=env,stdout=stdout,stderr=stderr)
  job.update(exit_code=result.returncode,elapsed_seconds=round(time.monotonic()-start,3));print(job['species'],job['background'],'exit',result.returncode,flush=True);return job
 with cf.ThreadPoolExecutor(max_workers=2) as pool:
  for fut in cf.as_completed([pool.submit(run,j) for j in jobs]):
   manifest['runs'].append(fut.result());(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 manifest['status']='completed' if all(j['exit_code']==0 for j in manifest['runs']) else 'failed'
 manifest['output_sha256']={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'}
 (out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 if manifest['status']!='completed':raise RuntimeError('STREME failed; consult retained logs')
if __name__=='__main__':main()
