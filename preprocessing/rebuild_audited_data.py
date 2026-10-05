"""Rebuild PR01-02 sequence tables and isolated FASTA inputs; never run motif tools."""
from __future__ import annotations
import csv, hashlib, json, random, re
from collections import Counter, defaultdict
from html.parser import HTMLParser
from pathlib import Path
import edlib
from preprocessing.provenance import sha256_file
from preprocessing.source_adapters import parse_gff3_attributes
from preprocessing.shuffle import shuffle_sequence, dinucleotide_counts
from experiments.prepare_tjupan_motif_inputs import SPECIES

ROOT=Path(__file__).resolve().parents[1]
SEED=20261005

def rc(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]
def canonical(s): return min(s,rc(s))
def genome(path): return ''.join(x.strip() for x in path.read_text().splitlines() if not x.startswith('>')).upper()
def tss_window(g,t,strand):
    start,end=(t-60,t+20) if strand=='+' else (t-20,t+60)
    if start<1 or end>len(g): raise ValueError('TSS window outside reference')
    s=g[start-1:end]
    return s if strand=='+' else rc(s)
def related(a,b):
    # Full 81nt query fitting either orientation, <=8 edit operations; target
    # aligned span necessarily >=73nt. Includes substitutions, shifts and indels.
    return any(edlib.align(x,y,mode='HW',task='distance',k=8)['editDistance']>=0
               for x,y in ((a,b),(rc(a),b),(b,a),(rc(b),a)))
def kmers(s): return {canonical(s[i:i+8]) for i in range(len(s)-7)}
def table(path,rows,columns=None):
    path.parent.mkdir(parents=True,exist_ok=True)
    if columns is None: columns=sorted(set().union(*(r.keys() for r in rows))) if rows else ['sequence_id']
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=columns,delimiter='\t');w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in r.items()})
def fasta(path,rows):
    with path.open('w',encoding='ascii') as f:
        for r in rows:f.write(f">{r['sequence_id']}\n{r['sequence']}\n")

class PageParser(HTMLParser):
    def __init__(self): super().__init__(convert_charrefs=True);self.rows=[];self.row=None;self.cell=None;self.bold=0
    def handle_starttag(self,tag,attrs):
        if tag=='tr':self.row=[]
        elif tag in ('td','th') and self.row is not None:self.cell={'text':'','bold':[]};self.bold=0
        elif tag in ('b','strong'):self.bold+=1
    def handle_endtag(self,tag):
        if tag in ('b','strong'):self.bold=max(0,self.bold-1)
        elif tag in ('td','th') and self.cell is not None:
            self.row.append(self.cell);self.cell=None
        elif tag=='tr' and self.row is not None:
            self.rows.append(self.row);self.row=None
    def handle_data(self,data):
        if self.cell is not None:
            before=len(re.sub(r'\s+','',self.cell['text']));clean=re.sub(r'\s+','',data)
            if self.bold:self.cell['bold']+=list(range(before+1,before+len(clean)+1))
            self.cell['text']+=data

def dbtbs_records(path,g):
    p=PageParser();p.feed(path.read_text(encoding='latin1'));direction=None;active=False;out=[]
    for rownum,cells in enumerate(p.rows):
        texts=[' '.join(c['text'].split()) for c in cells];norm=[re.sub('[^a-z]','',s.lower()) for s in texts]
        if 'direction' in norm and 'genomeposition' in norm:active=False;continue
        if direction is None and len(texts)>=4 and texts[2] in ('+','-') and re.fullmatch(r'\d+\.\.\d+',texts[3]):direction=texts[2]
        if any('bindingfactor' in x for x in norm) and any('bindingseq' in x for x in norm):active=True;continue
        if not active or len(texts)<6 or texts[1].lower()!='promoter' or not re.fullmatch(r'sig(?:ma)?[a-z0-9]+',texts[0],re.I):continue
        factor,_,loc,absolute,_,evidence=texts[:6];s=re.sub(r'\s+','',cells[4]['text']).upper()
        r={'source':'dbtbs_bsub','species':'bacillus_subtilis','record_id':f'{path.name}:{rownum}','source_file':(str(path.relative_to(ROOT)) if path.is_relative_to(ROOT) else str(path)),'source_file_sha256':sha256_file(path),'sigma':factor,'confidence':'experimental_page','evidence':evidence,'strand':direction,'location':loc,'absolute_position':absolute,'bold_positions':cells[4]['bold'],'assembly':'NC_000964.2','reason':''}
        lm=re.fullmatch(r'([+-]?\d+):([+-]?\d+)',loc);am=re.fullmatch(r'(\d+)\.\.(\d+)',absolute)
        if not lm or not am or direction is None:r['reason']='Location_ND_or_unresolved';out.append(r);continue
        a,b=map(int,lm.groups());lo,hi=map(int,am.groups());z=lambda x:x if x<0 else x-1
        t1=lo-z(a) if direction=='+' else hi+z(a);t2=hi-z(b) if direction=='+' else lo+z(b)
        expected=g[lo-1:hi];expected=expected if direction=='+' else rc(expected)
        if t1!=t2:r['reason']='inconsistent_location_endpoints'
        elif s!=expected:r['reason']='reference_binding_sequence_mismatch'
        elif cells[4]['bold'] and cells[4]['bold']!=[1-z(a)]:r['reason']='multiple_or_conflicting_TSS_marks'
        else:r.update(sequence=tss_window(g,t1,direction),tss=t1,input_orientation='transcription_forward',orientation_evidence='binding_sequence_and_location_verified',tier='core',label=1)
        out.append(r)
    return out

def collect():
    records=[];excluded=[];assets={}
    def asset(p):assets[str(p.relative_to(ROOT))]=sha256_file(p)
    for species,folder in SPECIES.items():
        base=ROOT/'data/raw/tju_pan_promoter/reg_and_gen/Datasets'/folder
        rows=list(csv.DictReader((base/'Dataset.csv').open()));byid={r['seq_id']:r for r in rows};conflicts=set()
        for name in ('train','dev','test'):
            p=base/f'{name}.csv';asset(p)
            for r in csv.DictReader(p.open()):
                if any(r[k]!=byid[r['seq_id']][k] for k in ('seq','label','seq_type')):conflicts.add(r['seq_id'])
        pp=base/'positive_samples.csv';asset(pp);positive={r['seq_id']:r['seq'] for r in csv.DictReader(pp.open())}
        if positive!={r['seq_id']:r['seq'] for r in rows if r['label']=='1'}:raise ValueError('positive_samples differs')
        asset(base/'Dataset.csv')
        for r in rows:
            s=r['seq'].upper();rec={'source':'tjupan_'+species,'species':species,'record_id':r['seq_id'],'sequence':s,'label':int(r['label']),'sigma':'unknown','confidence':'upstream_label','tier':'main','source_file':str((base/'Dataset.csv').relative_to(ROOT)),'source_file_sha256':assets[str((base/'Dataset.csv').relative_to(ROOT))],'tss':None,'strand':None,'assembly':None,'input_orientation':'upstream_transcription_forward_accepted','orientation_evidence':'PromLoop_author_statement_accepted_by_project','reason':''}
            if r['seq_id'] in conflicts:rec['reason']='split_record_conflict'
            elif len(s)!=81 or set(s)-set('ACGT'):rec['reason']='invalid_81nt_ACGT'
            (excluded if rec['reason'] else records).append(rec)
    gp=ROOT/'data/raw/regulondb_e_coli_k12_20260916.fna';rp=ROOT/'data/raw/regulondb_promoters_20260916.gff3';g=genome(gp);asset(gp);asset(rp)
    occurrences=Counter()
    for line in rp.read_text().splitlines():
        if not line or line.startswith('#'):continue
        f=line.split('\t');a=parse_gff3_attributes(f[8]);s=a.get('Sequence','').upper();name=a['name'];occurrences[name]+=1;t=int(f[3]);strand=f[6]
        r={'source':'regulondb_ecoli','species':'escherichia_coli','record_id':f'{name}#{occurrences[name]}','sequence':s,'label':1,'sigma':a.get('SigmaFactor') or 'unknown','confidence':a['Confidence'],'evidence':a.get('Evidence',''),'tier':'core' if a['Confidence'] in ('Strong','Confirmed') else 'extended','tss':t,'strand':strand,'assembly':f[0],'input_orientation':'transcription_forward','orientation_evidence':'exact_reference_window_match','source_file':str(rp.relative_to(ROOT)),'source_file_sha256':assets[str(rp.relative_to(ROOT))],'reason':''}
        if len(s)!=81 or set(s)-set('ACGT') or t<1:r['reason']='missing_or_invalid_sequence_TSS'
        elif s!=tss_window(g,t,strand):r['reason']='reference_window_mismatch'
        (excluded if r['reason'] else records).append(r)
    gp=ROOT/'data/raw/dbtbs_reference_20261005/NC_000964.2.fasta';asset(gp);g=genome(gp)
    for p in sorted((ROOT/'data/raw/dbtbs_v4.1_20261005/pages').glob('*.html')):
        asset(p)
        for r in dbtbs_records(p,g):(excluded if r['reason'] else records).append(r)
    return records,excluded,assets

def run(output=ROOT/'data/processed/pr01_02_data_v2'):
    if output.exists():raise FileExistsError(output)
    output.mkdir(parents=True)
    records,excluded,assets=collect();groups=defaultdict(list)
    for r in records:
        key=canonical(r['sequence']) if r['source'].startswith('tjupan') else r['sequence']
        groups[(r['source'],r['tier'],r['label'],key)].append(r)
    masters=[];relations=[];mapping=[]
    for (source,tier,label,key),rr in sorted(groups.items()):
        s=rr[0]['sequence'];sid=source+'_'+tier+'_'+hashlib.sha256((str(label)+key).encode()).hexdigest()[:16]
        row={'sequence_id':sid,'source':source,'tier':tier,'species':rr[0]['species'],'label':label,'sequence':s,'length':81,'local_tss_index':61,'window_offsets':'[-60,20]','orientation':rr[0]['input_orientation'],'coordinate_evidence':rr[0]['orientation_evidence'],'assembly':rr[0]['assembly'],'source_records':len(rr),'genomic_tss_values':sorted({r['tss'] for r in rr if r['tss'] is not None}),'genomic_strands':sorted({r['strand'] for r in rr if r['strand'] is not None})}
        masters.append(row)
        for r in rr:
            mapping.append({k:v for k,v in r.items() if k not in ('sequence','reason')}|{'sequence_id':sid})
            relations.append({'sequence_id':sid,'source_record_id':r['record_id'],'sigma_raw':r['sigma'],'sigma_id':r['species']+':'+r['sigma'],'promoter_confidence':r['confidence'],'sigma_association_evidence':'source_label_not_independently_verified','evidence':r.get('evidence','')})
    # Build one global leakage graph per species, across BOTH promoter/control
    # and all sources/tiers. An 8mer index is exhaustive for <=8 edits on 81nt;
    # nine disjoint 9nt blocks guarantee at least one preserved >=8nt substring.
    parent=list(range(len(masters)))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    def union(a,b):
        a,b=find(a),find(b)
        if a!=b:parent[max(a,b)]=min(a,b)
    links=[];by_species=defaultdict(list)
    for i,r in enumerate(masters):by_species[r['species']].append(i)
    for species,indices in sorted(by_species.items()):
        posting=defaultdict(list)
        for i in indices:
            candidates=set()
            for k in kmers(masters[i]['sequence']):candidates.update(posting[k])
            for j in sorted(candidates):
                if related(masters[i]['sequence'],masters[j]['sequence']):union(i,j);links.append({'a':masters[i]['sequence_id'],'b':masters[j]['sequence_id'],'method':'full_query_fitting_edit_distance_le8_both_orientations'})
            for k in kmers(masters[i]['sequence']):posting[k].append(i)
        # Native records at the same promoter locus or strongly overlapping
        # genomic windows remain on the same side even if their bases differ.
        loci=defaultdict(list)
        for i in indices:
            for t in masters[i]['genomic_tss_values']:
                if masters[i]['assembly']:loci[masters[i]['assembly']].append((t,i))
        for points in loci.values():
            points.sort()
            for n,(t,i) in enumerate(points):
                for u,j in points[n+1:]:
                    if u-t>8:break
                    if i!=j:union(i,j)
        print('leakage_graph',species,'sequences',len(indices),'links',len(links),flush=True)
    members=defaultdict(list)
    for i in range(len(masters)):members[find(i)].append(i)
    for indices in members.values():
        group=min(masters[i]['sequence_id'] for i in indices);bucket=int(hashlib.sha256(f'{SEED}:{group}'.encode()).hexdigest()[:8],16)%10
        split='discovery' if bucket<7 else 'development' if bucket==7 else 'holdout'
        for i in indices:masters[i].update(leakage_group=group,split=split)
    # Tests are a separate set of components. Source-disjoint external checks
    # are the auxiliary HOLDOUT components with no TJU positive/control member.
    has_tju={r['leakage_group'] for r in masters if r['source'].startswith('tjupan')}
    test_rows=[r for r in masters if r['split']!='discovery'];test_seq={canonical(r['sequence']) for r in test_rows};test_index=defaultdict(set)
    for r in test_rows:
        for k in kmers(r['sequence']):test_index[(r['species'],k)].add(r['sequence'])
    summary={};output_files={};rng=random.Random(SEED)
    for source,tier in sorted({(r['source'],r['tier']) for r in masters}):
        cohort=[r for r in masters if r['source']==source and r['tier']==tier];d=output/source/tier;d.mkdir(parents=True)
        table(d/'sequence_master.tsv',cohort);ids={r['sequence_id'] for r in cohort};table(d/'sigma_associations.tsv',[r for r in relations if r['sequence_id'] in ids]);table(d/'source_record_mapping.tsv',[r for r in mapping if r['sequence_id'] in ids])
        stats={}
        for split in ('discovery','development','holdout'):
            for label in (1,0):
                rows=[r for r in cohort if r['split']==split and r['label']==label]
                if not rows:continue
                name=f'{split}.{"positive" if label else "natural_control"}.fasta';fasta(d/name,rows);stats[name]=len(rows)
            positive=[r for r in cohort if r['split']==split and r['label']==1]
            null=[]
            for r in positive:
                for attempt in range(128):
                    s=shuffle_sequence(r['sequence'],rng,method='dinucleotide')
                    if s==r['sequence']:continue
                    if split=='discovery':
                        candidates=set()
                        for k in kmers(s):candidates.update(test_index[(r['species'],k)])
                        if any(related(s,t) for t in candidates):continue
                    if dinucleotide_counts(s)!=dinucleotide_counts(r['sequence']):raise AssertionError('invalid shuffle')
                    null.append({'sequence_id':r['sequence_id']+'__dinucleotide_null','sequence':s});break
                else:raise ValueError('cannot generate isolated dinucleotide control '+r['sequence_id'])
            if null:fasta(d/f'{split}.dinucleotide_null.fasta',null)
        external=[r for r in cohort if r['split']=='holdout' and r['label']==1 and r['leakage_group'] not in has_tju]
        if not source.startswith('tjupan') and external:fasta(d/'external_validation.positive.fasta',external);stats['external_validation.positive.fasta']=len(external)
        summary[source+'/'+tier]=stats
    table(output/'sequence_master.tsv',masters);table(output/'sigma_associations.tsv',relations);table(output/'source_record_mapping.tsv',mapping);table(output/'exclusions.tsv',excluded);table(output/'leakage_links.tsv',links)
    membership=defaultdict(set)
    for r in masters:membership[r['leakage_group']].add(r['split'])
    assert all(len(s)==1 for s in membership.values())
    assert not ({canonical(r['sequence']) for r in masters if r['split']=='discovery'}&test_seq)
    for p in sorted(output.rglob('*')):
        if p.is_file():output_files[str(p.relative_to(output))]=sha256_file(p)
    manifest={'schema':'pr01-02-data-v2','status':'data_rebuilt_no_motif_experiment','seed':SEED,'primary':'TJU six species; upstream TSS alignment accepted','auxiliary':'RegulonDB and DBTBS sigma annotation/analysis and source-disjoint holdout checks','tss_index_1based':61,'raw_sha256':assets,'output_sha256':output_files,'excluded_reasons':dict(Counter(r['reason'] for r in excluded)),'counts':summary,'statistical_units':len(masters),'source_records_retained':len(records),'sigma_associations':len(relations),'leakage_graph_links':len(links),'cross_split_leakage_groups':0,'known_limits':['TJU alignment relies on upstream statement, not per-record genomic verification','DBTBS .2 matches accepted binding intervals; exact 2005 genome flank identity unverified','promoter confidence does not prove sigma association confidence','small sigma groups require post-split eligibility check','fitting edit-distance threshold is an operational isolation rule, not proof of biological independence']}
    (output/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True);return manifest

if __name__=='__main__':run()
