"""Build short-element references from database annotations, never full promoters.

Source snapshots and sequence-level provenance stay in ignored tmp/. Element
classification rules are frozen before comparisons; no discovery logo is used.
"""
from pathlib import Path
import json,re,csv,hashlib,collections
from .rebuild_audited_data import PageParser,dbtbs_records,genome as read_genome,rc
ROOT=Path(__file__).resolve().parents[1]
BASE=ROOT/'tmp/reference_library_v1_20261008'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
class RedParser(PageParser):
 def __init__(self):super().__init__();self.red=False
 def handle_starttag(self,tag,attrs):
  super().handle_starttag(tag,attrs)
  if tag in ('td','th') and self.cell is not None:self.cell['red']=[];self.cell['links']=[];self.red=False
  d=dict(attrs)
  if tag=='font' and d.get('color','').lower() in ('red','#ff0000'):self.red=True
  if tag=='a' and self.cell is not None:self.cell['links'].append(d.get('href',''))
 def handle_endtag(self,tag):
  if tag=='font':self.red=False
  super().handle_endtag(tag)
 def handle_data(self,data):
  if self.cell is not None and self.red:
   offset=len(re.sub(r'\s+','',self.cell['text']));self.cell.setdefault('red',[]).extend(range(offset,offset+len(re.sub(r'\s+','',data))))
  super().handle_data(data)
def write_table(p,rows):
 keys=sorted(set(k for r in rows for k in r))
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=keys,delimiter='\t',lineterminator='\n');w.writeheader()
  for r in rows:w.writerow({k:json.dumps(v) if isinstance(v,(dict,list)) else v for k,v in r.items()})
def main():
 rules={'scope':['E.coli Sigma70','B.subtilis SigA'],'minimum_independent_sites':10,'independence':'one site per genomic interval, identical short sequences at distinct coordinates retained','alignment':'source annotated 6bp boxes; DBTBS contiguous red 6bp spans only, minus10 start -15..-7, minus35 start -38..-30','pseudocount':'total 1 uniform (0.25 per base per column)','orientation':'source transcription direction verified against exact genome interval','evidence':'S/C RegulonDB with experimental PMID and no COMP promoter citations; DBTBS red cis element with source PMID and experimental code; promoter evidence does not prove individual nucleotide perturbation','no_consensus_or_whole_promoter':True}
 (BASE/'build_rules.json').write_text(json.dumps(rules,indent=2))
 sites=[];excluded=[];g=read_genome(ROOT/'data/raw/regulondb_e_coli_k12_20260916.fna')
 for r in json.loads((BASE/'regulondb_promoters.json').read_text()):
  sigma=(r.get('bindsSigmaFactor') or {}).get('abbreviatedName');cites=r.get('citations') or [];experimental=[c for c in cites if (c.get('evidence') or {}).get('type') in ('S','C') and (c.get('evidence') or {}).get('code','').startswith(('EXP-','EV-EXP-')) and (c.get('publication') or {}).get('pmid')];inferred=any((c.get('evidence') or {}).get('code','').startswith('COMP') for c in cites)
  for b in r.get('boxes') or []:
   s=(b.get('sequence') or '').upper();a=b.get('leftEndPosition');z=b.get('rightEndPosition');strand=r.get('operon_strand');reason=''
   if sigma!='sigma70' or b.get('type') not in ('minus10','minus35'):reason='outside_initial_sigma70_scope'
   elif r.get('confidenceLevel') not in ('S','C') or not experimental or inferred:reason='insufficient_or_mixed_computational_evidence'
   elif len(s)!=6 or set(s)-set('ACGT') or not a or not z or z-a+1!=6:reason='invalid_short_box'
   elif strand not in ('forward','reverse') or (g[a-1:z] if strand=='forward' else rc(g[a-1:z]))!=s:reason='coordinate_sequence_mismatch'
   if reason:excluded.append({'source':'RegulonDB','record_id':r['_id'],'box_type':b.get('type'),'reason':reason});continue
   sites.append({'species':'ecoli','sigma':'Sigma70','element_type':b['type'],'source':'RegulonDB','release':'14.5.0 API snapshot 2026-10-08','record_id':r['_id'],'promoter':r['name'],'sequence':s,'assembly':'NC_000913.3','left':a,'right':z,'strand':strand,'PMIDs':[c['publication']['pmid'] for c in experimental],'evidence':experimental,'evidence_tier':'curated_short_box_in_experimentally_supported_promoter_no_COMP_evidence','classification':'explicit_database_box_type'})
 g=read_genome(ROOT/'data/raw/dbtbs_reference_20261005/NC_000964.2.fasta')
 for p in sorted((ROOT/'data/raw/dbtbs_v4.1_20261005/pages').glob('*.html')):
  core={r['record_id']:r for r in dbtbs_records(p,g) if not r.get('reason')};parser=RedParser();parser.feed(p.read_text(errors='replace'));promoter_table=False
  for idx,cells in enumerate(parser.rows):
   ts=[re.sub(r'\s+',' ',c['text']).strip() for c in cells]
   if any('Binding' in t for t in ts) and any('Regulation' in t for t in ts):promoter_table=True;continue
   if not promoter_table or len(ts)<6 or ts[:2]!=['SigA','Promoter']:continue
   record=f'{p.name}:{idx}';r=core.get(record)
   if not r:continue
   c=cells[4];s=re.sub(r'\s+','',c['text']).upper();pmids=re.findall(r'(?:list_uids=\+?|pubmed/)(\d+)',' '.join(cells[5].get('links',[])));method=ts[5]
   if not pmids or not re.search(r'\b(PE|RG|RO|SDM|FT|S1|DP)\b',method) or re.search(r'\bHM\b',method):continue
   a=int(ts[2].split(':')[0]);z=lambda x:x if x<0 else x-1
   groups=[]
   for j in sorted(set(c.get('red',[]))):
    if not groups or j!=groups[-1][-1]+1:groups.append([j])
    else:groups[-1].append(j)
   for span in groups:
    if len(span)!=6:continue
    rel=z(a)+span[0];kind='minus10' if -15<=rel<=-7 else 'minus35' if -38<=rel<=-30 else None
    if not kind:continue
    lo=int(r['absolute_position'].split('..')[0]);hi=int(r['absolute_position'].split('..')[1]);left=lo+span[0] if r['strand']=='+' else hi-span[-1];right=lo+span[-1] if r['strand']=='+' else hi-span[0]
    sites.append({'species':'bsub','sigma':'SigA','element_type':kind,'source':'DBTBS','release':'4.1','record_id':record,'promoter':p.stem,'sequence':s[span[0]:span[-1]+1],'assembly':'NC_000964.2','left':left,'right':right,'strand':r['strand'],'PMIDs':pmids,'evidence':method,'evidence_tier':'literature_red_cis_element_in_experimentally_supported_promoter','classification':'coordinate_classification_of_source_red_span','relative_start':rel})
 write_table(BASE/'reference_sites.tsv',sites);write_table(BASE/'reference_exclusions.tsv',excluded)
 groups=collections.defaultdict(dict)
 for s in sites:groups[(s['species'],s['sigma'],s['element_type'])][(s['assembly'],s['left'],s['right'],s['strand'])]=s
 motifs=[];txt='MEME version 5\n\nALPHABET= ACGT\n\nstrands: + -\n\nBackground letter frequencies\nA 0.25 C 0.25 G 0.25 T 0.25\n\n'
 for (sp,sigma,kind),records in sorted(groups.items()):
  ss=list(records.values());n=len(ss);mid=f'{sp}_{sigma}_{kind}';tier='eligible_curated_reference' if n>=10 else 'descriptive_small_n'
  motifs.append({'motif_id':mid,'species':sp,'sigma_family':sigma,'element_type':kind,'n_sites':n,'source_database':ss[0]['source'],'evidence_tier':tier,'short_site_evidence_caveat':'source annotation plus promoter experiment, not guaranteed base-specific mutagenesis for every box','PMIDs':sorted(set(p for s in ss for p in s['PMIDs'])),'source_records':[s['record_id'] for s in ss]})
  txt+=f'MOTIF {mid}\nletter-probability matrix: alength= 4 w= 6 nsites= {n} E= 0\n'
  for pos in range(6):txt+=' '.join(f'{(sum(s["sequence"][pos]==b for s in ss)+.25)/(n+1):.8f}' for b in 'ACGT')+'\n'
  txt+='\n'
 (BASE/'reference_library_v1.meme').write_text(txt);write_table(BASE/'reference_library_v1.tsv',motifs)
 m={'status':'built_and_sequence_verified','rules':rules,'motifs':motifs,'unavailable':['Ecoli UP element: no verified short-site asset','other four species: no reference admitted in v1'],'source_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(BASE.glob('regulondb*.json'))},'output_sha256':{p.name:sha(p) for p in [BASE/'reference_sites.tsv',BASE/'reference_exclusions.tsv',BASE/'reference_library_v1.meme',BASE/'reference_library_v1.tsv',BASE/'build_rules.json']},'reuse_limit':'local only; database redistribution terms unresolved'};(BASE/'manifest.json').write_text(json.dumps(m,indent=2)+'\n');print(json.dumps(motifs,indent=2)[:2000])
if __name__=='__main__':main()
