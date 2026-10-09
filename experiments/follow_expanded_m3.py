"""Follow completed discovery jobs; every emitted PWM gets independent FIMO.

All per-sequence data stay in tmp. Existing original17/holdout outputs untouched.
"""
from pathlib import Path
import json,csv,hashlib,os,time,subprocess,re,collections,xml.etree.ElementTree as ET
from .run_expanded_m3 import background,fa
ROOT=Path(__file__).resolve().parents[1];NATIVE=ROOT/'tmp/m3_expanded_20261009';OUT=ROOT/'tmp/m3_expanded_followup_20261009';BASE=ROOT/'data/processed/pr01_02_data_v2';BIN=ROOT/'tmp/meme-suite-5.5.9/bin'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def readfa(path):
 rows=[];sid=None;s=''
 for line in path.read_text().splitlines():
  if line.startswith('>'):
   if sid:rows.append((sid,s))
   sid=line[1:].split()[0];s=''
  else:s+=line.strip()
 if sid:rows.append((sid,s))
 return rows
def readtable(p):
 if not p.exists():return []
 return list(csv.DictReader((s for s in p.read_text().splitlines() if s and not s.startswith('#')),delimiter='\t'))
def table(p,rows):
 keys=sorted(set(k for r in rows for k in r)) if rows else ['status']
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def parse_candidates(native,method,key):
 root=ET.parse(native/('meme.xml' if method=='meme' else 'streme.xml')).getroot();rows=[];matrices=[]
 for i,m in enumerate(root.findall('./motifs/motif'),1):
  a=m.attrib;global_id=f'{key}__M{i:02d}';width=int(a['width']);valid=True
  if method=='meme':
   mx=[{v.attrib['letter_id']:float(v.text) for v in pos.findall('value')} for pos in m.findall('./probabilities/alphabet_matrix/alphabet_array')];p=a.get('p_value','NA');e=a.get('e_value','NA');consensus=a.get('name','');native_id=a['id'];sites=a.get('sites','2');valid=float(p)>=0
  else:
   mx=[{b:float(pos.attrib[b]) for b in 'ACGT'} for pos in m.findall('pos')];p=a.get('test_pvalue','NA');e=a.get('test_evalue','NA');consensus=a['id'].split('-',1)[-1];native_id=a['id'];sites=a.get('npassing','2');valid=int(root.find('.//test_positives').attrib['count'])>=5 and int(root.find('.//test_negatives').attrib['count'])>=5 and 0<=float(p)<=1
  assert len(mx)==width and all(abs(sum(v.values())-1)<1e-5 for v in mx)
  rows.append({'case':key,'method':method,'motif_id':global_id,'native_motif_id':native_id,'consensus':consensus,'width':width,'native_p':p,'native_E':e,'internal_test_valid':valid,'native_significant':valid and float(e if method=='meme' else p)<=.05,'inference':'exploratory_expansion_not_replacement_of_original17'})
  matrices.append((global_id,native_id,width,sites,mx))
 return rows,matrices
def write_pwm(path,matrices):
 text='MEME version 5\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n'
 for mid,orig,w,n,mx in matrices:
  text+=f'MOTIF {mid} {orig}\nletter-probability matrix: alength= 4 w= {w} nsites= {n} E= 0\n'+'\n'.join(' '.join(f'{row[b]:.8f}' for b in 'ACGT') for row in mx)+'\n\n'
 path.write_text(text)
def execute(cmd,log):
 env=dict(os.environ);env['PATH']=str(BIN)+os.pathsep+env.get('PATH','')
 with log.open('w') as f:r=subprocess.run(list(map(str,cmd)),stdout=f,stderr=f,cwd=ROOT,env=env)
 return {'command':list(map(str,cmd)),'exit_code':r.returncode,'log':str(log.relative_to(ROOT))}
def scan(case,pwm,input,bg,dest,mode,split):
 cmd=[BIN/'fimo','--o',dest,'--bfile',bg,'--motif-pseudo','0.1','--max-stored-scores','5000000']
 if mode in ('strict','strict_sigma_pool'):cmd+=['--qv-thresh','--thresh','0.05']
 else:cmd+=['--thresh','0.0001']
 cmd +=[pwm,input];log=dest.parent/(dest.name+'.log');result=execute(cmd,log);result.update(case=case,split=split,mode=mode,input_sha256=sha(input),pwm_sha256=sha(pwm),background_sha256=sha(bg))
 if result['exit_code']==0:
  text=log.read_text();result['truncation_warning']=bool(re.search(r'discard|truncat',text,re.I));assert not result['truncation_warning']
  rows=readtable(dest/'fimo.tsv');result.update(site_count=len(rows),sequence_hit_count=len({r['sequence_name'] for r in rows}));table(dest.parent/(dest.name+'_primary_sites.tsv'),primary_sites(rows))
 return result
def primary_sites(rows):
 best={}
 for r in rows:
  key=(r['sequence_name'],r['motif_id']);order=lambda r:(float(r.get('q-value',1)),-float(r['score']),int(r['start']),int(r['stop']),r['strand'])
  if key not in best or order(r)<order(best[key]):best[key]=r
 return [{**r,'TSS_relative_center':(int(r['start'])+int(r['stop']))/2-61} for r in best.values()]
def main():
 OUT.mkdir(exist_ok=True);state={'status':'running','cases':[],'tools':{t:sha(BIN/t) for t in ['fimo','tomtom','sea']},'analysis_boundary':'new discovery and new inspections of existing holdout are exploratory; no threshold optimization; diagnostic raw-p only discovery/development','multiple_testing':'FIMO native site q per motif per case/split; downstream BH separate per split and comparison purpose across all newly/reused expanded cases; original17 untouched','raw_p_cutoff':1e-4,'q_cutoff':.05,'sigma_specificity':'same source other single-sigma groups; no shuffle-based specificity claim'};candidates=[];processed=set()
 if (OUT/'manifest.json').exists():
  state=json.loads((OUT/'manifest.json').read_text());candidates=readtable(OUT/'candidates.tsv');processed={c['key'] for c in state['cases']}
  for r in candidates:
   for k in ['internal_test_valid','native_significant']:r[k]=str(r[k]).lower()=='true'
 def save():
  table(OUT/'candidates.tsv',candidates);(OUT/'manifest.json').write_text(json.dumps(state,indent=2)+'\n')
 def follow(job,reused=False):
  source=job['source'];group=job['group'];method=job['method'];control=job['control'];key=job.get('key',f'{source}__{group}__{control}__{method}');folder=OUT/key;folder.mkdir();native=ROOT/(job['native_output'] if reused else job['output']);cohort=BASE/source/'main' if source.startswith('tjupan') else BASE/source/'core' if group=='overall' else BASE/source/'core/sigma'/group
  rows,mats=parse_candidates(native,method,key);candidates.extend(rows);pwm=folder/'all_reported.meme';write_pwm(pwm,mats);case={'key':key,'source':source,'group':group,'method':method,'control':control,'native_output':str(native.relative_to(ROOT)),'cohort':str(cohort.relative_to(ROOT)),'reused_discovery':reused,'motif_count':len(rows),'native_significant':sum(r['native_significant'] for r in rows),'pwm_sha256':sha(pwm),'native_xml_sha256':sha(native/('meme.xml' if method=='meme' else 'streme.xml')),'FIMO':[],'comparisons':[]}
  if not mats:case['status']='no_pwm_tool_found_no_motif';state['cases'].append(case);save();return
  bg=folder/'discovery_control_zero_order.bg';background(cohort/f'discovery.{control}.fasta',bg,0)
  for split in ['discovery','development','holdout']:
   p=cohort/f'{split}.positive.fasta';n=cohort/f'{split}.{control}.fasta'
   if not fa(p) or not fa(n):case['FIMO'].append({'split':split,'status':'skipped_empty_partition','positive_count':fa(p),'negative_count':fa(n)});continue
   input=folder/f'{split}.combined.fasta';input.write_bytes(p.read_bytes()+n.read_bytes());case['FIMO'].append(scan(key,pwm,input,bg,folder/(split+'_strict'),'strict',split))
   if split!='holdout':
    case['FIMO'].append(scan(key,pwm,input,bg,folder/(split+'_raw_p'),'diagnostic_raw_p',split))
    sea=[BIN/'sea','--p',p,'--n',n,'--m',pwm,'--bfile',bg,'--order','0','--seed','20261005','--hofract','0.1','--thresh','1000000','--noseqs','--o',folder/(split+'_sea')];result=execute(sea,folder/(split+'_sea.log'));result.update(split=split,purpose='sequence_level_enrichment_diagnostic_only');case['comparisons'].append(result)
  # Only appropriate species reference; no sigma70 architecture for sigma54.
  ref=ROOT/'tmp/reference_library_v1_20261008/reference_library_v1.meme';prefix='ecoli_' if source in ['regulondb_ecoli','tjupan_escherichia_coli'] and group in ['overall','Sigma70'] else 'bsub_' if source in ['dbtbs_bsub','tjupan_bacillus_subtilis'] and group in ['overall','SigA'] else None
  if ref.exists() and prefix:
   parts=re.split(r'(?m)^MOTIF ',ref.read_text());target=folder/'known_reference.meme';target.write_text(parts[0]+''.join('MOTIF '+s for s in parts[1:] if s.split()[0].startswith(prefix)));cmd=[BIN/'tomtom','--o',folder/'known_tomtom','--dist','pearson','--min-overlap','5','--thresh','1',pwm,target];result=execute(cmd,folder/'known_tomtom.log');result['purpose']='curated_reference_similarity_small_library_descriptive';case['comparisons'].append(result)
  else:case['reference_status']='no_applicable_qualified_reference_all_motifs_retained'
  # Sigma comparison input contains actual other sigma positives, no null DNA.
  if group!='overall':
   for split in ['discovery','development','holdout']:
    pool=[];maprows=[]
    for sg in sorted((BASE/source/'core/sigma').iterdir()):
     p=sg/f'{split}.positive.fasta'
     if p.exists():
      for sid,seq in readfa(p):pool.append((sid,seq));maprows.append({'sequence_id':sid,'sigma':sg.name,'target_label':int(sg.name==group)})
    if not pool:continue
    inp=folder/f'{split}.sigma_pool.fasta';inp.write_text(''.join(f'>{sid}\n{seq}\n' for sid,seq in pool));table(folder/f'{split}.sigma_pool_mapping.tsv',maprows);r=scan(key,pwm,inp,bg,folder/(split+'_sigma_pool'),'strict_sigma_pool',split);case['FIMO'].append(r)
  case['status']='completed' if all(r.get('exit_code',0)==0 for r in case['FIMO']+case['comparisons']) else 'completed_with_failures';state['cases'].append(case);save();print(key,'FIMO done; motifs',len(rows),'native significant',case['native_significant'],flush=True)
 while True:
  try:m=json.loads((NATIVE/'manifest.json').read_text())
  except (FileNotFoundError,json.JSONDecodeError):time.sleep(3);continue
  # Start existing STREME branches while MEME is still computing.
  for j in m['reuse']:
   key=f'{j["source"]}__{j["group"]}__{j["control"]}__{j["method"]}'
   if key not in processed:follow(j,True);processed.add(key)
  for j in m['jobs']:
   if j['key'] in processed or j['status']=='planned':continue
   if j['status']=='completed':follow(j)
   else:state['cases'].append({'key':j['key'],'status':'discovery_failed_no_pwm','exit_code':j['exit_code'],'source':j['source'],'group':j['group'],'method':j['method']});save()
   processed.add(j['key'])
  if m['status']!='running' and len(processed)==len(m['jobs'])+len(m['reuse']):break
  time.sleep(10)
 state['status']='completed_with_recorded_failures' if any(c['status'] not in ['completed','no_pwm_tool_found_no_motif'] for c in state['cases']) else 'completed';state['output_sha256']={str(p.relative_to(OUT)):sha(p) for p in OUT.rglob('*') if p.is_file() and p.name!='manifest.json'};save()
if __name__=='__main__':main()
