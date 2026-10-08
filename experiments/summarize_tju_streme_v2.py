"""Summarize discovery-only significance and PWM robustness; freeze main list."""
from pathlib import Path
import xml.etree.ElementTree as ET
import csv,json,hashlib,subprocess,os,re
ROOT=Path(__file__).resolve().parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def table(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]) if rows else ['motif_id'],delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def matrix(m):return tuple(tuple(float(p.attrib[k]) for k in 'ACGT') for p in m.findall('pos'))
def rc(mx):return tuple(tuple(reversed(r)) for r in reversed(mx))
def main():
 out=ROOT/'results/motif/tju_streme_v2_20261008';manifest=json.loads((out/'manifest.json').read_text());assert manifest['status']=='completed'
 summaries=[];motifs=[];chosen=[];comparisons=[];tool=ROOT/'tmp/meme-suite-5.5.9/bin/tomtom'
 for species in sorted(out.glob('tjupan_*')):
  lists={}
  for bg in ['natural','dinucleotide']:
   x=ET.parse(species/bg/'streme.xml').getroot();ms=x.findall('./motifs/motif');assert len(ms)<=10;lists[bg]=ms
   assert all(5<=int(m.attrib['width'])<=15 for m in ms)
   for m in ms:
    a=m.attrib;motifs.append({'species':species.name,'background':bg,'motif_id':a['id'],'width':a['width'],'internal_test_p':a['test_pvalue'],'internal_test_E':a['test_evalue'],'internal_pos_hits':a['test_pos_count'],'internal_neg_hits':a['test_neg_count'],'p_le_0_05':float(a['test_pvalue'])<=.05})
   summaries.append({'species':species.name,'background':bg,'reported_motifs':len(ms),'p_le_0_05':sum(float(m.attrib['test_pvalue'])<=.05 for m in ms),'E_le_0_05':sum(float(m.attrib['test_evalue'])<=.05 for m in ms)})
  selected=[];seen=set()
  for m in sorted(lists['natural'],key=lambda m:(float(m.attrib['test_pvalue']),int(m.attrib['id'].split('-')[0]))):
   if float(m.attrib['test_pvalue'])>.05:continue
   mx=matrix(m);key=min(mx,rc(mx))
   if key in seen:continue
   seen.add(key);selected.append(m)
   if len(selected)==5:break
  ids={m.attrib['id'] for m in selected}
  for rank,m in enumerate(selected,1):chosen.append({'species':species.name,'rank':rank,'motif_id':m.attrib['id'],'width':m.attrib['width'],'internal_test_p':m.attrib['test_pvalue'],'internal_test_E':m.attrib['test_evalue'],'selection_source':'natural_only'})
  text=(species/'natural/streme.txt').read_text();chunks=re.split(r'(?m)^MOTIF ',text);kept=[c for c in chunks[1:] if c.split()[0] in ids];(species/'selected_main.meme').write_text(chunks[0]+''.join('MOTIF '+c for c in kept));assert len(kept)==len(ids)
  if selected and lists['dinucleotide']:
   dest=species/'robustness_tomtom';cmd=[str(tool),'--oc',str(dest),'--dist','pearson','--min-overlap','5','--thresh','1',str(species/'selected_main.meme'),str(species/'dinucleotide/streme.txt')]
   env=dict(os.environ);env['PATH']=str(tool.parent)+os.pathsep+env.get('PATH','')
   with (species/'tomtom.stdout.log').open('w') as a,(species/'tomtom.stderr.log').open('w') as b:result=subprocess.run(cmd,cwd=ROOT,env=env,stdout=a,stderr=b)
   assert result.returncode==0,species
   comparisons.append({'species':species.name,'command':cmd,'exit_code':result.returncode,'interpretation':'PWM similarity only; small target library; no biological validation; absence not proof of non-enrichment'})
 table(out/'motif_summary.tsv',motifs);table(out/'run_summary.tsv',summaries);table(out/'selected_main_motifs.tsv',chosen)
 manifest['main_candidates']=chosen;manifest['run_summary']=summaries;manifest['robustness_comparisons']=comparisons;manifest['tomtom_version']=subprocess.check_output([str(tool),'--version'],text=True).strip();manifest['tomtom_sha256']=sha(tool);manifest['output_sha256']={str(p.relative_to(out)):sha(p) for p in sorted(out.rglob('*')) if p.is_file() and p.name!='manifest.json'};(out/'manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
 print(json.dumps(summaries,indent=2));print('Selected natural candidates:',len(chosen))
if __name__=='__main__':main()
