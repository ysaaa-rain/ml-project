#!/usr/bin/env python3
"""Read-only supplemental M3 RQ1 audit; public output contains aggregates only."""
from __future__ import annotations
import csv, json, hashlib, math, re
from collections import defaultdict, Counter
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT=Path(__file__).resolve().parents[1]
V2=ROOT/'data/processed/pr01_02_data_v2'
REF=ROOT/'tmp/reference_library_v1_20261008'
FOLLOW=ROOT/'tmp/m3_expanded_followup_20261009'
PUBLIC=ROOT/'results/motif/m3_expanded_20261009'
OUT=ROOT/'results/motif/m3_rq1_annotation_audit_20261010'
LOCAL=ROOT/'tmp/m3_rq1_annotation_audit_20261010'
OUT.mkdir(parents=True,exist_ok=True); LOCAL.mkdir(parents=True,exist_ok=True)

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(1<<20),b''): h.update(b)
 return h.hexdigest()
def readtsv(p):
 with open(p,encoding='utf-8',newline='') as f:return list(csv.DictReader(f,delimiter='\t'))
def writetsv(p,rows,fields=None):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 if fields is None: fields=list(rows[0]) if rows else []
 with open(p,'w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader();w.writerows(rows)
def fasta(p):
 d={};name=None;seq=[]
 if not Path(p).exists():return d
 for ln in open(p,encoding='utf-8'):
  ln=ln.strip()
  if ln.startswith('>'):
   if name is not None:d[name]=''.join(seq).upper()
   name=ln[1:].split()[0];seq=[]
  elif ln:seq.append(ln)
 if name is not None:d[name]=''.join(seq).upper()
 return d
def json_arr(s):
 try:return json.loads(s)
 except:return []
def wc(path):return sha(path)
def plot_save(fig,name):
 fig.tight_layout();fig.savefig(OUT/(name+'.png'),dpi=180,bbox_inches='tight');fig.savefig(OUT/(name+'.svg'),bbox_inches='tight');plt.close(fig)

# Frozen source inventory; do not modify any input or prior receipt.
manifest=json.load(open(FOLLOW/'manifest.json'))
cases=manifest['cases']; cand=readtsv(PUBLIC/'candidates.tsv'); similarity=readtsv(PUBLIC/'PWM_similarity.tsv'); clusters=readtsv(PUBLIC/'cluster_tests.tsv')
refsites=readtsv(REF/'reference_sites.tsv')
regapi=json.load(open(REF/'regulondb_promoters.json'))
regapi_byname=defaultdict(list)
for r in regapi: regapi_byname[r['name']].append(r)
# Master and source-record indexes. Sequences are used only to verify annotated element intervals after ID/TSS mapping.
masters={}; srcmaps={}; refmap={}; sigma_map=defaultdict(set)
for source,sub in [('RegulonDB','regulondb_ecoli'),('DBTBS','dbtbs_bsub')]:
 base=V2/sub/'core'
 for row in readtsv(base/'sequence_master.tsv'): masters[row['sequence_id']]=row
 maps=readtsv(base/'source_record_mapping.tsv'); srcmaps[source]=maps
 for row in readtsv(base/'sigma_associations.tsv'):sigma_map[row['sequence_id']].add(row['sigma_raw'])
 # native record IDs are authoritative; no DNA matching is used to create this map.
 for m in maps:
  if source=='RegulonDB':
   name=re.sub(r'#\d+$','',m['record_id'])
   refmap[(source,name)]=m
  else: refmap[(source,m['record_id'])]=m

# Reference annotation projection: establish eligibility by provenance ID, assembly, single TSS/strand, and exact element subsequence.
apiid={r['_id']:r for r in regapi}
projected=[]; projection_reasons=Counter()
for r in refsites:
 source=r['source']; key=(source,r['promoter'] if source=='RegulonDB' else r['record_id'])
 m=refmap.get(key)
 if not m: projection_reasons[(source,'source_record_id_missing')]+=1;continue
 master=masters.get(m['sequence_id'])
 if master is None:projection_reasons[(source,'master_missing')]+=1;continue
 try:
  tss=int(m['tss']); left=int(r['left']); right=int(r['right'])
 except:projection_reasons[(source,'coordinate_missing')]+=1;continue
 if m.get('assembly')!=r.get('assembly') or master.get('assembly')!=r.get('assembly'):
  projection_reasons[(source,'assembly_mismatch')]+=1;continue
 if source=='RegulonDB':
  api=apiid.get(r['record_id']); at=api.get('transcriptionStartSite',{}) if api else {}
  if not api or not at or at.get('leftEndPosition')!=at.get('rightEndPosition') or int(at['leftEndPosition'])!=tss:
   projection_reasons[(source,'API_TSS_mismatch')]+=1;continue
 strand=m['strand']; refstrand=r['strand']
 if strand not in ['+','-'] or (source=='RegulonDB' and {'+':'forward','-':'reverse'}[strand]!=refstrand) or (source=='DBTBS' and strand!=refstrand):
  projection_reasons[(source,'strand_mismatch')]+=1;continue
 tssvals=set(json_arr(master['genomic_tss_values'])); strands=set(json_arr(master['genomic_strands']))
 if tssvals!={tss} or strands!={strand}:
  projection_reasons[(source,'ambiguous_master_TSS_or_strand')]+=1;continue
 if strand=='+': start,leftlocal=left-tss+61,right-tss+61
 else: start,leftlocal=tss-right+61,tss-left+61
 start=int(start); stop=int(leftlocal)
 if start<1 or stop>81 or stop<start:
  projection_reasons[(source,'outside_81bp_window')]+=1;continue
 if master['sequence'][start-1:stop]!=r['sequence'].upper():
  projection_reasons[(source,'element_sequence_mismatch')]+=1;continue
 if source=='DBTBS' and r.get('relative_start'):
  try:
   if int(r['relative_start'])!=start-61:projection_reasons[(source,'relative_location_mismatch')]+=1;continue
  except:projection_reasons[(source,'relative_location_unparseable')]+=1;continue
 projected.append({'source':source,'sequence_id':m['sequence_id'],'sigma':r['sigma'],'element_type':r['element_type'],'left':left,'right':right,'strand':strand,'start':start,'stop':stop,'record_id':r['record_id'],'promoter':r['promoter'],'assembly':r['assembly'],'PMIDs':r['PMIDs'],'seq':master['sequence']})

# Map all positive cohorts directly through input FASTA names -> frozen master IDs; coordinate transfer only source-matched records.
source_for=lambda c: 'RegulonDB' if c['source']=='regulondb_ecoli' else ('DBTBS' if c['source']=='dbtbs_bsub' else 'TJU')
cohort_data={}; cohort_master={}
for c in cases:
 cohort=ROOT/c['cohort']; key=c['key']
 # discovery FASTA consists only positive cohort sequences; its IDs also key to sequence_master IDs.
 pos=fasta(cohort/'discovery.positive.fasta'); cohort_data[key]=pos
 cm={}
 if c['source'] in ('regulondb_ecoli','dbtbs_bsub'):
  cm={sid:masters[sid] for sid in pos if sid in masters}
 cohort_master[key]=cm

# candidate significance keyed by run and motif
cand_by=defaultdict(dict)
for x in cand: cand_by[x['case']][x['motif_id']]=x

# Mapping coverage per run/group; based on FASTA positive IDs and source provenance, never TJU exact-sequence lookup.
mapping_rows=[]
for c in cases:
 key=c['key']; src=source_for(c); pos=cohort_data[key]
 if src=='TJU':
  mapping_rows.append(dict(case=key,source=src,sigma_group=c['group'],positive_sequences=len(pos),coordinate_unique_sequences=0,qualified_element_annotated_sequences=0,coordinate_sequences_without_qualified_annotation=0,excluded_no_genomic_provenance=len(pos),excluded_ambiguous_context=0,excluded_sigma_not_applicable=0,coverage='NA(no per-sequence genomic coordinate; PromLoop standardized TSS-relative sequences)'))
  continue
 group=c['group']; applicable='Sigma70' if src=='RegulonDB' else 'SigA'
 coords=set(); annseq=set();amb=set(); noctx=set(); nonsigma=set()
 for sid in pos:
  m=cohort_master[key].get(sid)
  if not m: noctx.add(sid);continue
  tv=set(json_arr(m['genomic_tss_values'])); st=set(json_arr(m['genomic_strands']))
  if len(tv)!=1 or len(st)!=1:amb.add(sid);continue
  coords.add(sid)
  if group!='overall' and group!=applicable:
   nonsigma.add(sid);continue
  if any(a['source']==src and a['sequence_id']==sid and a['element_type'] in ('minus10','minus35') for a in projected):annseq.add(sid)
 mapping_rows.append(dict(case=key,source=src,sigma_group=group,positive_sequences=len(pos),coordinate_unique_sequences=len(coords),qualified_element_annotated_sequences=len(annseq),coordinate_sequences_without_qualified_annotation=len((coords-nonsigma)-annseq),excluded_no_genomic_provenance=len(noctx),excluded_ambiguous_context=len(amb),excluded_sigma_not_applicable=len(nonsigma),coverage='mapped' if annseq else 'NA(no eligible source-matched element annotation)'))
writetsv(OUT/'mapping_coverage.tsv',mapping_rows)

# Find reported FIMO sites and real-interval overlaps. Count distinct intervals/sites, never pairwise duplicates.
site_cache={}; overlap_rows=[]; tss_rows=[]
# Annotated element windows from the verified 6-bp site position, descriptive only.
def expected_window(tp): return (-12.5,-4.5) if tp=='minus10' else (-35.5,-27.5)
for c in cases:
 key=c['key']; src=source_for(c); pos=cohort_data[key]; bymot=cand_by[key]
 for run in c.get('FIMO',[]):
  if 'mode' not in run: continue
  split=run['split']; mode=run['mode']
  if split not in ('discovery','development','holdout') or mode not in ('strict','diagnostic_raw_p'): continue
  fp=FOLLOW/key/f'{split}_{"strict" if mode=="strict" else "raw_p"}'/'fimo.tsv'
  if not fp.exists(): continue
  hits=[]
  split_positive=fasta(ROOT/c['cohort']/f'{split}.positive.fasta')
  for r in csv.DictReader(open(fp,encoding='utf-8'),delimiter='\t'):
   try:p=float(r['p-value']);q=float(r['q-value']);st=int(r['start']);en=int(r['stop'])
   except:continue
   if mode=='diagnostic_raw_p' and not p<1e-4:continue
   mid=r['motif_id']; candidate=bymot.get(mid,{})
   seqid=r['sequence_name'];
   if seqid not in split_positive: continue
   # FIMO may prepend an identifier namespace in some builds; current published data are exact IDs.
   hit={'case':key,'split':split,'mode':mode,'motif':mid,'native_sig':candidate.get('native_significant','False')=='True','seqid':seqid,'start':st,'stop':en,'strand':r['strand'],'center':(st+en)/2-61,'p':p,'q':q,'score':float(r['score'])}
   hits.append(hit)
  site_cache[(key,split,mode)]=hits
  # TSS position summaries: one row per motif including all and native-significant-only scopes.
  if split in ('discovery','development'):
   for scope,sel in [('all_reported',hits),('native_significant_only',[h for h in hits if h['native_sig']])]:
    motif_universe=sorted(bymot) if scope=='all_reported' else sorted(mid for mid,v in bymot.items() if v.get('native_significant')=='True')
    for motif in motif_universe:
     mh=[h for h in sel if h['motif']==motif]
     # Find corresponding candidate sequences denominator; all positives in the cohort split.
     splitfa=split_positive
     if not splitfa and split=='development': splitfa={}
     N=len(splitfa)
     win_app=(src=='RegulonDB' and c['group']=='Sigma70') or (src=='DBTBS' and c['group']=='SigA')
     for strand in ['+','-','both']:
      hs=mh if strand=='both' else [h for h in mh if h['strand']==strand]
      primary={}
      for h in hs:
       old=primary.get(h['seqid'])
       if old is None or (h['q'], -h['score'],h['start'],h['stop'],h['strand'])<(old['q'],-old['score'],old['start'],old['stop'],old['strand']):primary[h['seqid']]=h
      for tp in ['minus10','minus35']:
       lo,hi=expected_window(tp)
       nwin=sum(lo<=h['center']<=hi for h in primary.values())
       tss_rows.append(dict(case=key,source=src,sigma_group=c['group'],method=c['method'],control=c['control'],split=split,site_mode=mode,analysis_set=scope,motif_id=motif,hit_strand=strand,element_window=tp,center_window_relative_to_TSS=f'[{lo},{hi}]',window_applicable=win_app,positive_sequences=N,sites=len(hs),hit_sequences=len(set(h['seqid'] for h in hs)),primary_hit_sequences=len(primary),primary_hits_in_window=(nwin if win_app else 'NA(group has no uniform known sigma architecture)'),primary_hit_window_rate=((nwin/len(primary)) if primary and win_app else ('NA(no primary hits)' if win_app else 'NA(group has no uniform known sigma architecture)')),all_sites_center_median=(sorted(h['center'] for h in hs)[len(hs)//2] if hs else 'NA(no sites)'),denominator_note='one best q, then score, then leftmost site per sequence and motif; both-strand handled separately; q strict or raw p<1e-4 as labelled'))
  if src=='TJU':continue
  # Per-motif projection/overlap; zero-hit motifs are retained with explicit zero numerators.
  applicable='Sigma70' if src=='RegulonDB' else 'SigA'
  for scope in ('all_reported','native_significant_only'):
   motif_universe=sorted(bymot) if scope=='all_reported' else sorted(mid for mid,v in bymot.items() if v.get('native_significant')=='True')
   for mid in motif_universe:
    native=bymot[mid].get('native_significant')=='True'
    for tp in ['minus10','minus35']:
     anns=[a for a in projected if a['source']==src and a['element_type']==tp and (c['group']=='overall' or c['group']==applicable) and (c['group']=='overall' or a['sigma']==c['group'])]
     splitids=set(split_positive)
     coordids={sid for sid,m in cohort_master[key].items() if sid in splitids and len(set(json_arr(m['genomic_tss_values'])))==1 and len(set(json_arr(m['genomic_strands'])))==1 and (c['group']=='overall' or c['group'] in sigma_map[sid])}
     anns=[a for a in anns if a['sequence_id'] in splitids]
     annintervals={(a['sequence_id'],a['left'],a['right'],a['strand']) for a in anns}
     annseq={a['sequence_id'] for a in anns}
     testhits=[h for h in hits if h['motif']==mid and h['seqid'] in coordids and (scope=='all_reported' or h['native_sig'])]
     pred=set(); overlap=set()
     for h in testhits:
      m=cohort_master[key].get(h['seqid'])
      if not m:continue
      tss=int(next(iter(json_arr(m['genomic_tss_values'])))); gs=next(iter(json_arr(m['genomic_strands'])))
      if gs=='+': a=tss+h['start']-61;b=tss+h['stop']-61
      else: a=tss-(h['stop']-61);b=tss-(h['start']-61)
      ghit=gs if h['strand']=='+' else ('-' if gs=='+' else '+')
      ik=(h['seqid'],a,b,ghit);pred.add(ik)
      for ai in anns:
       if ai['sequence_id']==h['seqid'] and a<=ai['right'] and ai['left']<=b:
        overlap.add((ik,(ai['sequence_id'],ai['left'],ai['right'],ai['strand'])))
     unique_overlap_hits={x[0] for x in overlap}; unique_overlap_ann={x[1] for x in overlap}
     nseq=len(annseq);nann=len(annintervals);nhit=len(pred);nov=len(unique_overlap_hits)
     exact=sum(1 for ph,ah in overlap if ph[1]==ah[1] and ph[2]==ah[2])
     recip=sum(1 for ph,ah in overlap if (min(ph[2],ah[2])-max(ph[1],ah[1])+1)/(ph[2]-ph[1]+1)>=.5 and (min(ph[2],ah[2])-max(ph[1],ah[1])+1)/(ah[2]-ah[1]+1)>=.5)
     overlap_rows.append(dict(case=key,source=src,sigma_group=c['group'],method=c['method'],control=c['control'],split=split,site_mode=mode,analysis_set=scope,motif_id=mid,native_significant=native,element_type=tp,eligible_sequences=len(splitids) if (c['group']=='overall' or c['group']==applicable) else 0,coordinate_mapped_sequences=len(coordids),mapped_annotated_sequences=nseq,annotated_unique_sites=nann,predicted_unique_site_intervals=nhit,overlap_unique_predicted_intervals_ge1bp=nov,overlap_predicted_hit_rate=(nov/nhit if nhit else 'NA(no predicted hits)'),overlap_unique_annotation_sites_ge1bp=len(unique_overlap_ann),annotation_site_recovery_rate=(len(unique_overlap_ann)/nann if nann else 'NA(no eligible annotation sites)'),exact_boundary_predicted_intervals=exact,reciprocal_overlap_ge50_predicted_intervals=recip,literal_genomic_strand_matching_pairs=sum(1 for ph,ah in overlap if ph[3]==ah[3]),strand_note='input hit strand mapped through source genomic strand; PWM strand is not inherently biological; >=1bp overlap; dedup by genomic hit and annotated intervals; reference subset overlaps discovery source so descriptive'))
writetsv(OUT/'element_overlap.tsv',overlap_rows)
# Original frozen TJU 17 STREME candidates: reuse existing full raw-p diagnostic and strict FIMO tables.
# They are not refit and retain their original run IDs and site threshold provenance.
origrows=readtsv(PUBLIC/'original17_diagnostic_distribution.tsv')
orig_species=sorted(set(x['species'] for x in origrows))
for species in orig_species:
 cohort=V2/species/'main'; pos_by_split={s:fasta(cohort/f'{s}.positive.fasta') for s in ('discovery','development')}
 for split in ('discovery','development'):
  for mode,fp in [('diagnostic_raw_p',ROOT/'tmp/fimo_diagnostic_20261009'/species/split/'fimo.tsv'),('strict',ROOT/'tmp/m3_followup_20261008'/species/f'{split}_fimo'/'fimo.tsv')]:
   if not fp.exists(): continue
   bymot=defaultdict(list)
   for r in csv.DictReader(open(fp,encoding='utf-8'),delimiter='\t'):
    try:p=float(r['p-value']);q=float(r['q-value']);st=int(r['start']);en=int(r['stop']);score=float(r['score'])
    except:continue
    if mode=='diagnostic_raw_p' and not p<1e-4:continue
    sid=r['sequence_name']
    if sid not in pos_by_split[split]:continue
    bymot[r['motif_id']].append({'sid':sid,'start':st,'stop':en,'strand':r['strand'],'center':(st+en)/2-61,'q':q,'score':score})
   # Keep zero-hit old candidates by reading the precomputed per-motif inventory for this species.
   motif_ids=sorted(set(x['motif_id'] for x in origrows if x['species']==species))
   for motif in motif_ids:
    hits=bymot.get(motif,[])
    for strand in ['+','-','both']:
     hs=hits if strand=='both' else [h for h in hits if h['strand']==strand]
     primary={}
     for h in hs:
      old=primary.get(h['sid'])
      if old is None or (h['q'],-h['score'],h['start'],h['stop'],h['strand'])<(old['q'],-old['score'],old['start'],old['stop'],old['strand']):primary[h['sid']]=h
     for tp,(lo,hi) in [('minus10',(-12.5,-4.5)),('minus35',(-35.5,-27.5))]:
      nwin=sum(lo<=h['center']<=hi for h in primary.values())
      tss_rows.append(dict(case=f'original17_{species}',source='TJU',sigma_group='unknown',method='STREME',control='natural_control',split=split,site_mode=mode,analysis_set='original17_native_significant',motif_id=motif,hit_strand=strand,element_window=tp,center_window_relative_to_TSS=f'[{lo},{hi}]',window_applicable=False,positive_sequences=len(pos_by_split[split]),sites=len(hs),hit_sequences=len(set(h['sid'] for h in hs)),primary_hit_sequences=len(primary),primary_hits_in_window=(nwin if win_app else 'NA(group has no uniform known sigma architecture)'),primary_hit_window_rate=((nwin/len(primary)) if primary and win_app else ('NA(no primary hits)' if win_app else 'NA(group has no uniform known sigma architecture)')),all_sites_center_median=(sorted(h['center'] for h in hs)[len(hs)//2] if hs else 'NA(no sites)'),denominator_note='frozen original17; existing q<=0.05 strict table or raw p<1e-4 diagnostic table; one best q then score then leftmost per sequence'))
# Preserve absent development FIMO outputs (e.g., source-recorded insufficient negative count) as NA rather than zero hits.
for c in cases:
 for run in c.get('FIMO',[]):
  if run.get('split')!='development' or run.get('mode') is not None:continue
  for mid,x in cand_by[c['key']].items():
   scopes=['all_reported']+(['native_significant_only'] if x.get('native_significant')=='True' else [])
   for scope in scopes:
    for strand in ['+','-','both']:
     for tp in ['minus10','minus35']:
      tss_rows.append(dict(case=c['key'],source=source_for(c),sigma_group=c['group'],method=c['method'],control=c['control'],split='development',site_mode='strict',analysis_set=scope,motif_id=mid,hit_strand=strand,element_window=tp,center_window_relative_to_TSS='NA(FIMO not run)',window_applicable=False,positive_sequences=run.get('positive_count','NA'),sites='NA(FIMO not run: empty split)',hit_sequences='NA',primary_hit_sequences='NA',primary_hits_in_window='NA',primary_hit_window_rate='NA',all_sites_center_median='NA(FIMO not run)',FIMO_status=run.get('status','not_run'),denominator_note='Development partition empty; FIMO was not run; this is not a zero-hit result'))
for row in tss_rows:row.setdefault('FIMO_status','completed_existing_output')
writetsv(OUT/'tss_window_distribution.tsv',tss_rows)

# Holdout three exploratory findings: reproduce rows from frozen inferential table and add sequence composition / position descriptors.
three=[]
known_three=[('tjupan_bacillus_subtilis__overall__natural_control__meme','tjupan_bacillus_subtilis__overall__natural_control__meme__M01'),('tjupan_bacillus_subtilis__overall__natural_control__meme','tjupan_bacillus_subtilis__overall__natural_control__meme__M02'),('tjupan_bradyrhizobium__overall__dinucleotide_null__meme','tjupan_bradyrhizobium__overall__dinucleotide_null__meme__M02')]
for key,mid in known_three:
 stat=next(x for x in clusters if x['case']==key and x['motif_id']==mid and x['split']=='holdout' and x['family'] in ('natural_presence','paired_null_presence'))
 c=next(x for x in cases if x['key']==key); cohort=ROOT/c['cohort']; pos=fasta(cohort/'holdout.positive.fasta'); neg=fasta(cohort/'holdout.natural_control.fasta' if c['control']=='natural_control' else cohort/'holdout.dinucleotide_null.fasta')
 # compositions descriptive on full holdout positives and controls. For dinucleotide null, paired sequence statistics are also reported.
 def comp(d):
  seq=''.join(d.values()); n=len(seq);return {b:(seq.count(b)/n if n else None) for b in 'ACGT'}|{'GC':((seq.count('G')+seq.count('C'))/n if n else None),'AT':((seq.count('A')+seq.count('T'))/n if n else None)}
 pc,nc=comp(pos),comp(neg)
 hs=[]
 for mode in ('strict','diagnostic_raw_p'):
  hs.extend(site_cache.get((key,'holdout',mode),[]))
 motif_hits=[h for h in hs if h['motif']==mid]
 centers=[h['center'] for h in motif_hits]
 three.append(dict(case=key,motif_id=mid,consensus=next(x['consensus'] for x in cand if x['case']==key and x['motif_id']==mid),source=c['source'],method=c['method'],control=c['control'],family=stat['family'],positive_clusters=stat['positive_clusters'],positive_hit_clusters=stat['positive_hits'],positive_sequences_beyond_one_per_cluster=len(pos)-int(stat['positive_clusters']),control_clusters=stat['control_clusters'],control_hit_clusters=stat['control_hits'],control_sequences_beyond_one_per_cluster=len(neg)-int(stat['control_clusters']),mixed_clusters_excluded=stat['mixed_clusters_excluded'],rate_difference=stat['rate_difference'],CI95_lower=stat['rate_difference_CI95_lower'],CI95_upper=stat['rate_difference_CI95_upper'],p=stat['p'],BH_q=stat['BH_q'],positive_sequences=len(pos),control_sequences=len(neg),positive_GC=pc['GC'],control_GC=nc['GC'],positive_AT=pc['AT'],control_AT=nc['AT'],positive_A=pc['A'],positive_C=pc['C'],positive_G=pc['G'],positive_T=pc['T'],control_A=nc['A'],control_C=nc['C'],control_G=nc['G'],control_T=nc['T'],reported_strict_FIMO_sites=len(motif_hits),motif_site_center_median_relative_TSS=(sorted(centers)[len(centers)//2] if centers else 'NA(no sites)'),strict_site_centers_in_minus10_window=sum(-12.5<=v<=-4.5 for v in centers),strict_site_centers_in_minus35_window=sum(-35.5<=v<=-27.5 for v in centers),strict_site_centers_downstream_of_TSS=sum(v>0 for v in centers),homology_independent_test='frozen cluster-level test; see original cluster table; no posthoc change',annotation_relationship='TJU coordinates not transferable; Bradyrhizobium has no eligible coordinate-resolved reference annotation in this audit',interpretation='exploratory holdout association; sequence composition and positional summaries are descriptive'))
writetsv(OUT/'holdout_three_exploration.tsv',three)

# Reference coverage: distinguish qualified coordinate-resolved sites from text-only candidates.
coverage=[]
operons=json.load(open(REF/'regulondb_operons_all.json'))['data']['getAllOperon']['data']
up_candidates=[]
for op in operons:
 for tu in op.get('transcriptionUnits',[]):
  p=tu.get('promoter') or {}; note=p.get('note') or ''
  m=re.search(r'UP element, located between\s*(-?\d+)\s*and\s*(-?\d+) relative to the .*? transcriptional start site',note,re.I)
  if m:
   up_candidates.append({'name':p.get('name',''),'interval':(int(m.group(1)),int(m.group(2))),'tss':p.get('transcriptionStartSite',{}),'promoter':p,'record_ids':re.findall(r'RDBECOLIPRC\d+',note)})
for sp in ['escherichia_coli','bacillus_subtilis','acinetobacter_baumannii','bradyrhizobium','corynebacterium_diphtheriae','staphylococcus_aureus']:
 entries=[r for r in refsites if ('escherichia' in sp and r['source']=='RegulonDB') or ('bacillus' in sp and r['source']=='DBTBS')]
 for typ in ['minus10','minus35','UP','TFBS']:
  subset=[r for r in entries if r['element_type'].lower()==typ.lower()]
  pm=set()
  for r in subset:pm.update(json_arr(r['PMIDs']))
  candidates=up_candidates if (typ=='UP' and sp=='escherichia_coli') else []
  status='local qualified reference coverage' if subset else '尚无合格结构化位点（不表示公开数据库全局不存在）'
  detail=''
  candidate_pm=''
  if candidates:
   c0=candidates[0]
   detail=f"{c0['name']} relative interval {c0['interval']}; source note gives promoter sequence/TSS, but underlying element citation record has no resolved PMID in this snapshot; not included in qualified reference PWM"
   candidate_pm='RDBECOLIPRC22477 (PMID unresolved locally)'
   status='1 narrative coordinate candidate; not qualified for formal motif recovery until citation/evidence mapping is resolved'
  elif typ=='UP' and sp=='escherichia_coli':
   status='Narrative UP references exist, but no additional fully auditable coordinate record was resolved'
  coverage.append(dict(species=sp,element_type=typ,qualified_annotated_sites=len(subset),unique_PMIDs=len(pm),sources=','.join(sorted({r['source'] for r in subset})) if subset else ('RegulonDB narrative snapshot' if candidates else ''),candidate_count=len(candidates),candidate_reference=candidate_pm,candidate_context=detail,coverage_status=status,evidence_scope='curated short element in experimentally supported promoter; PMID/source record retained locally' if subset else 'Unqualified candidate is context only; do not calculate recovery rate'))
writetsv(OUT/'reference_coverage.tsv',coverage)

# Evidence-scope reconciliations, using already-frozen published aggregates.
# Pull full reported counts vs native-significant candidate counts; all 540 PWM receive FIMO, including nonsignificant ones.
evidence=[]
for method in ['meme','streme']:
 n=sum(1 for c in cases if c['method']==method)
 total=sum(int(c['motif_count']) for c in cases if c['method']==method)
 sig=sum(1 for r in cand if r['method']==method and r['native_significant']=='True')
 evidence.append(dict(scope='expanded_discovery',method=method,cases=n,all_reported_PWMs=total,native_significant_PWMs=sig,FIMO_policy='all reported PWMs scanned; downstream tables carry both scopes'))
# Existing summary facts explicitly preserve denominators and estimator distinction.
for d in [dict(scope='method_Tomtom',method='MEME_vs_STREME',pairs=61,q_le_0_05=61,both_sides_native_significant=9),dict(scope='known_reference_Tomtom',method='all_methods',pairs=320,q_le_0_05=15,both_native_significant_and_q=2),dict(scope='first_FIMO_diagnostic',method='old17',SEA_motif_partition_rows=34,SEA_native_significant=33,raw_p_lt_1e_4_sites=1494,strict_q_sites_original_aggregate=45,full_site_tables_strict=40,Staphylococcus_development_original_count=19,Staphylococcus_development_full_table_count=14,explanation='separate q-value estimator path/reservoir sizes; do not pool counts')]:evidence.append(d)
writetsv(OUT/'evidence_scopes.tsv',evidence)

# Four reproducible figures.
# A: TJU six species, natural control native significant counts.
spnames={'tjupan_bacillus_subtilis':'B. subtilis','tjupan_baumannii':'A. baumannii','tjupan_bradyrhizobium':'Bradyrhizobium','tjupan_diphtheria':'C. diphtheriae','tjupan_escherichia_coli':'E. coli','tjupan_staphylococcus':'Staphylococcus'}
order=list(spnames)
vals={m:[] for m in ['meme','streme']}
for s in order:
 for m in vals:
  c=next(x for x in cases if x['source']==s and x['control']=='natural_control' and x['method']==m)
  vals[m].append(sum(1 for x in cand if x['case']==c['key'] and x['native_significant']=='True'))
fig,ax=plt.subplots(figsize=(10,5));x=list(range(len(order)));w=.36
ax.bar([i-w/2 for i in x],vals['meme'],w,label='MEME (E≤0.05)');ax.bar([i+w/2 for i in x],vals['streme'],w,label='STREME (p≤0.05)')
ax.set_xticks(x,[spnames[s] for s in order],rotation=20,ha='right');ax.set_ylabel('Native-significant motifs');ax.set_title('TJU natural background: reported native thresholds differ');ax.legend();ax.text(.01,-.30,'Counts / species; MEME E-value and STREME internal p-value are method-specific, not directly comparable.',transform=ax.transAxes,fontsize=9)
plot_save(fig,'fig_a_tju_method_counts')
# B: similarity q matches and subset both native significant
pairs=[r for r in similarity if r.get('purpose')=='meme_vs_streme_exploratory' and r.get('q-value') and float(r['q-value'])<=.05]
cand_global={x['motif_id']:x for x in cand}
nboth=sum(1 for r in pairs if cand_global.get(r['Query_ID'],{}).get('native_significant')=='True' and cand_global.get(r['Target_ID'],{}).get('native_significant')=='True')
fig,ax=plt.subplots(figsize=(7,5));ax.bar(['Tomtom q≤0.05\nPWM pairs','Both motifs native\nsignificant'],[len(pairs),nboth],color=['#4c78a8','#f58518']);ax.set_ylabel('Pair count');ax.set_title('MEME–STREME PWM similarity');
for i,v in enumerate([len(pairs),nboth]):ax.text(i,v+.5,str(v),ha='center')
ax.text(.02,-.21,'Pairs are PWM similarities; not biological confirmation. Both-native subset uses each method’s own threshold.',transform=ax.transAxes,fontsize=9)
plot_save(fig,'fig_b_method_pwm_matches')
# C: verified annotations relative to TSS and element-overlap coverage
fig,axs=plt.subplots(1,2,figsize=(13,5.5))
for src,color in [('RegulonDB','#4c78a8'),('DBTBS','#f58518')]:
 dd=[a for a in projected if a['source']==src]
 centers=[(a['start']+a['stop'])/2-61 for a in dd]
 axs[0].scatter(centers,[0 if a['element_type']=='minus35' else 1 for a in dd],s=15,alpha=.55,label=f'{src} n={len(dd)}',color=color)
axs[0].set_yticks([0,1],['−35','−10']);axs[0].set_xlabel('Annotated element center relative to TSS (bp)');axs[0].set_title('Verified annotation positions');axs[0].legend()
# Show source-level annotation denominators and observed strict overlap results without adding motif-pair duplicates.
axs[1].axis('off');axs[1].set_title('Native-significant strict FIMO overlap (discovery)')
summary_lines=['0 predicted intervals overlapped a reference site by ≥1 bp.','', 'Eligible reference sites in discovery positives:']
for src,label in [('RegulonDB','E. coli'),('DBTBS','B. subtilis')]:
 parts=[]
 for tp,short in [('minus10','−10'),('minus35','−35')]:
  base=[a for a in projected if a['source']==src and a['element_type']==tp]
  # use the positive discovery IDs for the overall cohort
  c0=next(c for c in cases if c['source']==('regulondb_ecoli' if src=='RegulonDB' else 'dbtbs_bsub') and c['group']=='overall')
  ids=set(fasta(ROOT/c0['cohort']/'discovery.positive.fasta'))
  n=len({(a['sequence_id'],a['left'],a['right'],a['strand']) for a in base if a['sequence_id'] in ids and (src!='DBTBS' or a['sigma']=='SigA') and (src!='RegulonDB' or a['sigma']=='Sigma70')})
  parts.append(f'{short}: {n}')
 summary_lines.append(f'{label}:  '+ '   |   '.join(parts))
summary_lines += ['', 'FIMO strand was mapped through each promoter orientation;','overlap ignored strand and required ≥1 genomic base.','All totals use source-qualified, coordinate-projected sites.']
axs[1].text(.04,.88,'\n'.join(summary_lines),va='top',ha='left',fontsize=12,linespacing=1.6,transform=axs[1].transAxes)
fig.suptitle('Promoter-coordinate coverage; TJU has no genomic-coordinate overlap estimate');plot_save(fig,'fig_c_tss_and_annotation_coverage')
# D: frozen three cluster statistics with CI
fig,ax=plt.subplots(figsize=(9,5));yy=list(range(len(three)));effects=[float(r['rate_difference']) for r in three];lo=[float(r['CI95_lower']) for r in three];hi=[float(r['CI95_upper']) for r in three]
ax.errorbar(effects,yy,xerr=[[e-l for e,l in zip(effects,lo)],[h-e for e,h in zip(effects,hi)]],fmt='o',capsize=4,color='#4c78a8');ax.axvline(0,color='gray',lw=1);ax.set_yticks(yy,[r['consensus']+' · '+r['source'].replace('tjupan_','') for r in three]);ax.set_xlabel('Positive minus control hit-cluster rate (95% CI)');ax.set_title('Three exploratory holdout associations (frozen cluster tests)');plot_save(fig,'fig_d_holdout_cluster_effects')

# Hash inventory: all consumed inputs plus outputs. Store sensitive paths only in local manifest; public hash file is path-scrubbed.
inputs=[FOLLOW/'manifest.json',PUBLIC/'candidates.tsv',PUBLIC/'PWM_similarity.tsv',PUBLIC/'cluster_tests.tsv',REF/'reference_sites.tsv',REF/'regulondb_promoters.json',V2/'manifest.json']
inputs += [Path(__file__),PUBLIC/'original17_diagnostic_distribution.tsv']
for sp in ['tjupan_bacillus_subtilis','tjupan_baumannii','tjupan_bradyrhizobium','tjupan_diphtheria','tjupan_escherichia_coli','tjupan_staphylococcus']:
 for split in ['discovery','development']:
  for p in [ROOT/'tmp/fimo_diagnostic_20261009'/sp/split/'fimo.tsv',ROOT/'tmp/m3_followup_20261008'/sp/f'{split}_fimo'/'fimo.tsv']:
   if p.exists():inputs.append(p)
for c in cases:
 for s in ['discovery','development','holdout']:
  for mode in ['strict','raw_p']:
   p=FOLLOW/c['key']/f'{s}_{mode}'/'fimo.tsv'
   if p.exists():inputs.append(p)
for source in ['regulondb_ecoli','dbtbs_bsub']:
 inputs += [V2/source/'core'/'sequence_master.tsv',V2/source/'core'/'source_record_mapping.tsv',V2/source/'core'/'sigma_associations.tsv']
for c in cases:
 for p in [ROOT/c['cohort']/f'{split}.positive.fasta' for split in ['discovery','development','holdout']] + [ROOT/c['cohort']/f'{split}.{"natural_control" if c["control"]=="natural_control" else "dinucleotide_null"}.fasta' for split in ['discovery','development','holdout']]:
  if p.exists():inputs.append(p)
inputmap={str(p.relative_to(ROOT)):sha(p) for p in sorted(set(inputs))}
# Full local report retains exact local paths. Public hash manifest only hashes the immutable public/reference inputs and new public files.
outputs=[p for p in OUT.iterdir() if p.is_file() and p.name!='input_output_sha256.json']
local_receipt={'run_id':'m3_rq1_annotation_audit_20261010','read_only_inputs':inputmap,'output_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(outputs)},'counts':{'reference_sites':len(refsites),'projected_reference_sites':len(projected),'projection_failures':{f'{k[0]}:{k[1]}':v for k,v in projection_reasons.items()},'FIMO_cases':len(cases),'mapping_rows':len(mapping_rows),'overlap_rows':len(overlap_rows),'TSS_rows':len(tss_rows)}}
json.dump(local_receipt,open(LOCAL/'audit_receipt.json','w'),indent=2)
public_receipt={'run_id':'m3_rq1_annotation_audit_20261010','public_inputs_sha256':{str(p.relative_to(ROOT)):inputmap[str(p.relative_to(ROOT))] for p in [Path(__file__),FOLLOW/'manifest.json',PUBLIC/'candidates.tsv',PUBLIC/'PWM_similarity.tsv',PUBLIC/'cluster_tests.tsv',REF/'reference_sites.tsv',V2/'manifest.json']},'outputs_sha256':{str(p.relative_to(ROOT)):sha(p) for p in sorted(outputs) if p.suffix in ('.tsv','.png','.svg')},'aggregate_counts':local_receipt['counts']}
json.dump(public_receipt,open(OUT/'input_output_sha256.json','w'),indent=2)

# Summary with explicit evidence boundaries and reportable denominator facts.
lines=['# M3 RQ1 本地只读补证（2026-10-10）','',f"参考库位点：{len(refsites)}；按来源 ID、assembly、TSS、strand 和精确片段核验后可投影 {len(projected)} 条。投影失败：{dict((f'{k[0]}:{k[1]}',v) for k,v in projection_reasons.items()) or '0'}。",'','RegulonDB/DBTBS 的重合表按来源记录 ID 关联，并要求唯一 TSS/strand、参考版本一致、6bp 元件片段与冻结 81bp 窗口一致。一个预测位点与多个真实区间重叠时按唯一预测区间和唯一注释区间分别去重；≥1bp 为主定义。','',f"TJU 坐标覆盖逐条为 0/{sum(r['positive_sequences'] for r in mapping_rows if r['source']=='TJU')}；原因是 PromLoop CSV 未提供逐条基因组坐标，本报告只统计 FIMO 输入中的相对 TSS 位点，不借用同序列外部坐标。",'','## RQ1 初步结论（已观察 / 未证实 / 无法评估）','','**已观察：** 两种发现方法在六物种发现部分均有原生显著 motif；同源 PWM 比对产生 61 对 q≤0.05，其中两侧方法各自原生显著的 9 对。已知元件库 Tomtom 的 320 对中 15 对 q≤0.05，且同时满足 discovery 原生显著条件的 2 对均为 TJU B. subtilis SigA −10 相似候选。RegulonDB 与 DBTBS 注释元件坐标按来源 ID 成功投影，在发现序列内，RegulonDB −10/−35 的 18/13 个注释位点及 DBTBS SigA −10/−35 的 141/124 个注释位点中，native-significant strict FIMO 均未出现 ≥1bp 重叠（各运行/候选分母见逐 motif 表；注释位点重复展示不跨 motif 相加）。','','**未证实：** PWM 相似或 FIMO 区间相交不单独证明功能；三条预先报告的探索性 holdout 阳性不构成新的确认结论，原有确认性统计未被改写。B. subtilis 两条 natural-control holdout motif 的阳性集合 GC=32.6%、AT=67.4%，负样本 GC=46.4%、AT=53.6%；候选命中位点中心中位数分别为 −35.5bp 与 −24bp，说明组成差异可能参与自然负对照关联，需按探索性解释。位置窗口内命中率是描述统计。','','**无法评估：** TJU 缺逐条 genomic coordinate / strand，不能计算真实已知元件的坐标重合率；当前合格参考库没有合格 UP PWM；RegulonDB 本地 operon 注释含 1 条 guaBp UP 区间叙述（−59 至 −38），但底层 citation record `RDBECOLIPRC22477` 在快照中未解析出 PMID/证据对象，因此暂不纳入正式 recovery。结构化坐标化 TFBS 仍无合格覆盖，不计算恢复率。','',f"全输出 motif 分析范围：540 条报告 PWM；其中 MEME 原生显著 {sum(1 for r in cand if r['method']=='meme' and r['native_significant']=='True')}，STREME {sum(1 for r in cand if r['method']=='streme' and r['native_significant']=='True')}。位置与重合结果均在 `tss_window_distribution.tsv` / `element_overlap.tsv` 中按 `all_reported` 与 `native_significant_only`、split 和 FIMO 门槛区分。",'', '位置表的 −10/−35 中心窗口仅用于来源/σ分组确有统一参考架构的子集；TJU 缺逐条 σ、RegulonDB/DBTBS overall 是混合组，因此只给相对 TSS 中心分布并将架构窗口标 NA。DBTBS SigW development discovery 输入为空，相关 120 行标记为 FIMO 未运行，不记作零命中。','','## 证据边界与目录','','旧首次 FIMO 诊断的 SEA 34 项（33 项 native E≤0.05）、raw p<1e−4 的 1,494 个位点、原 strict 汇总 45 与完整 FIMO 表 40 分开保留。Staphylococcus development 的 19 与 14 是不同 q-value 估计路径/零分布规模下的两个证据口径，不合并。','','逐条映射及注释记录仅保存在本机 `tmp/m3_rq1_annotation_audit_20261010/`，本次公开文件只有不可逆汇总、短 motif 共识和统计图。输入/输出校验和见 `input_output_sha256.json`。']
(OUT/'summary.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')

# Include final summary hash in both receipts; receipt self-hash is intentionally omitted.
summary_path=OUT/'summary.md'
for receipt_path in [LOCAL/'audit_receipt.json',OUT/'input_output_sha256.json']:
 rec=json.load(open(receipt_path))
 target='output_sha256' if 'output_sha256' in rec else 'outputs_sha256'
 rec[target][str(summary_path.relative_to(ROOT))]=sha(summary_path)
 json.dump(rec,open(receipt_path,'w'),indent=2)
