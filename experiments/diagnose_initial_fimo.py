"""Read-only diagnostic of original 17 PWM, discovery/development only."""
from pathlib import Path
from collections import defaultdict
from array import array
import csv,json,hashlib,subprocess,os,re
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OLD=ROOT/'tmp/m3_followup_20261008';M3=ROOT/'results/motif/tju_streme_v2_20261008';OUT=ROOT/'tmp/fimo_diagnostic_20261009';BIN=ROOT/'tmp/meme-suite-5.5.9/bin'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 OUT.mkdir(exist_ok=False);summary=[];commands=[];bestrows=[];sea_runs=[]
 env=dict(os.environ);env['PATH']=str(BIN)+os.pathsep+env.get('PATH','')
 for species in sorted(M3.glob('tjupan_*')):
  name=species.name;folder=OUT/name;folder.mkdir();motif=species/'selected_main.meme';bg=OLD/name/'discovery_natural_zero_order.bg';cohort=ROOT/'data/processed/pr01_02_data_v2'/name/'main'
  for split in ['discovery','development']:
   inp=OLD/name/f'{split}.combined.fasta';cmd=[str(BIN/'fimo'),'--o',str(folder/split),'--bfile',str(bg),'--motif-pseudo','0.1','--max-stored-scores','5000000','--thresh','1',str(motif),str(inp)]
   log=folder/(split+'.log')
   with log.open('w') as f:r=subprocess.run(cmd,stdout=f,stderr=f,env=env,cwd=ROOT)
   assert r.returncode==0,(cmd,log);assert not re.search(r'discard|truncat',log.read_text(),re.I)
   commands.append({'species':name,'split':split,'command':cmd,'input_sha256':sha(inp),'motif_sha256':sha(motif),'background_sha256':sha(bg)})
   vals=defaultdict(lambda: {'p':array('d'),'q':array('d'),'score':array('d')});best={}
   with (folder/split/'fimo.tsv').open() as f:
    reader=csv.DictReader((s for s in f if s.strip() and not s.startswith('#')),delimiter='\t')
    for row in reader:
     p=float(row['p-value']);q=float(row['q-value']);score=float(row['score']);mid=row['motif_id'];v=vals[mid];v['p'].append(p);v['q'].append(q);v['score'].append(score);key=(row['sequence_name'],mid)
     if key not in best or score>float(best[key]['score']):best[key]=row
   bestrows.extend({'species':name,'split':split,**r} for r in best.values())
   for mid,v in vals.items():
    ps=np.array(v['p']);qs=np.array(v['q']);scores=np.array(v['score']);summary.append({'species':name,'split':split,'motif_id':mid,'reported_positions_p_lt_1':len(ps),'p_lt_1e_4':int((ps<1e-4).sum()),'p_lt_0_05':int((ps<.05).sum()),'q_le_0_05':int((qs<=.05).sum()),'minimum_p':float(ps.min()),'minimum_q':float(qs.min()),'maximum_score':float(scores.max()),'p_quantiles':np.quantile(ps,[0,.01,.1,.5,.9,1]).tolist(),'q_quantiles':np.quantile(qs,[0,.01,.1,.5,.9,1]).tolist()})
   old=list(csv.DictReader((s for s in (OLD/name/(split+'_fimo')/'fimo.tsv').read_text().splitlines() if s and not s.startswith('#')),delimiter='\t'))
   fresh=sum(x['q_le_0_05'] for x in summary if x['species']==name and x['split']==split)
   # Native pi0 estimator may differ with different stored p threshold; retain
   # discrepancies rather than replacing original strict outputs.
   commands[-1].update(original_strict_sites=len(old),diagnostic_full_table_strict_sites=fresh)
   seadir=folder/(split+'_sea');sea=[str(BIN/'sea'),'--p',str(cohort/f'{split}.positive.fasta'),'--n',str(cohort/f'{split}.natural_control.fasta'),'--m',str(motif),'--bfile',str(bg),'--order','0','--hofract','0.1','--seed','20261005','--thresh','1000000','--noseqs','--o',str(seadir)]
   with (folder/(split+'_sea.log')).open('w') as f:r=subprocess.run(sea,stdout=f,stderr=f,env=env,cwd=ROOT)
   sea_runs.append({'species':name,'split':split,'command':sea,'exit_code':r.returncode,'scope':'diagnostic_sequence_enrichment_only_no_holdout'})
   print(name,split,'p<1e-4',sum(x['p_lt_1e_4'] for x in summary if x['species']==name and x['split']==split),'original strict',len(old),'fulltable strict',fresh,flush=True)
 for name,rows in [('site_distribution_summary',summary),('best_sites',bestrows)]:
  with (OUT/(name+'.tsv')).open('w') as f:
   w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
 m={'status':'completed','scope':'original17; discovery/development only; no original files changed','background':'same original frozen zero-order natural controls','reported_distribution_scope':'p<1 positions, excludes p=1 sites','commands':commands,'sea':sea_runs,'no_truncation_warnings':True,'summary':summary,'output_sha256':{str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file()}};(OUT/'manifest.json').write_text(json.dumps(m,indent=2)+'\n')
if __name__=='__main__':main()
