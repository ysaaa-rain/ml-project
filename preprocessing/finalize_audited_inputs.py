"""Add sigma-specific views and verify all experimental input boundaries."""
import csv,json
from collections import defaultdict,Counter
from pathlib import Path
from preprocessing.rebuild_audited_data import ROOT,rc,canonical,kmers,related,table,fasta
from preprocessing.provenance import sha256_file

def read(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def readfa(p):
    out=[];name=None;s=''
    for line in p.read_text().splitlines():
        if line.startswith('>'):
            if name:out.append({'sequence_id':name,'sequence':s})
            name=line[1:];s=''
        else:s+=line.strip()
    if name:out.append({'sequence_id':name,'sequence':s})
    return out

def run(base=ROOT/'data/processed/pr01_02_data_v2'):
    masters=read(base/'sequence_master.tsv');rels=read(base/'sigma_associations.tsv');sigmas=defaultdict(set)
    for r in rels:
        if r['sigma_raw']!='unknown':sigmas[r['sequence_id']].add(r['sigma_raw'])
    summaries=[]
    for source in ('regulondb_ecoli','dbtbs_bsub'):
        d=base/source/'core';rows=[r for r in masters if r['source']==source and r['tier']=='core']
        for z in sorted(set().union(*(sigmas[r['sequence_id']] for r in rows))):
            chosen=[r for r in rows if sigmas[r['sequence_id']]=={z}];counts=Counter(r['split'] for r in chosen)
            eligibility='confirmation_candidate' if len(chosen)>=50 and counts['discovery']>=30 and counts['holdout']>=20 else 'exploratory_only' if len(chosen)>=20 else 'descriptive_only'
            summary={'source':source,'sigma':z,'single_sigma_sequences':len(chosen),'discovery':counts['discovery'],'development':counts['development'],'holdout':counts['holdout'],'eligibility':eligibility};summaries.append(summary)
            dest=d/'sigma'/z;dest.mkdir(parents=True,exist_ok=True)
            for split in ('discovery','development','holdout'):
                subset=[r for r in chosen if r['split']==split]
                if not subset:continue
                fasta(dest/f'{split}.positive.fasta',subset)
                lookup={r['sequence_id'].removesuffix('__dinucleotide_null'):r for r in readfa(d/f'{split}.dinucleotide_null.fasta')};fasta(dest/f'{split}.dinucleotide_null.fasta',[lookup[r['sequence_id']] for r in subset])
    for p in base.glob("*/*/external_validation.positive.fasta"):
        lookup = {r["sequence_id"].removesuffix("__dinucleotide_null"): r
                  for r in readfa(p.parent / "holdout.dinucleotide_null.fasta")}
        fasta(p.parent / "external_validation.dinucleotide_null.fasta",
              [lookup[r["sequence_id"]] for r in readfa(p)])
    table(base/'sigma_group_summary.tsv',summaries)
    auxiliary=defaultdict(list)
    for r in masters:
        if not r['source'].startswith('tjupan') and r['tier']=='core':auxiliary[(r['species'],canonical(r['sequence']))].append(r)
    annotations=[]
    for r in masters:
        if r['source'].startswith('tjupan') and r['label']=='1':
            for a in auxiliary[(r['species'],canonical(r['sequence']))]:
                annotations.append({'tju_sequence_id':r['sequence_id'],'auxiliary_sequence_id':a['sequence_id'],'source':a['source'],'same_orientation':r['sequence']==a['sequence'],'sigma_candidates':sorted(sigmas[a['sequence_id']]),'status':'exact_source_mapping_candidate_not_independent_validation'})
    table(base/'tju_sigma_annotation_candidates.tsv',annotations)
    # Independent verification, including generated controls: candidate index
    # never drops frequent kmers. No cross-side related sequence is allowed.
    train=defaultdict(set);test=defaultdict(set)
    for p in base.glob('*/*/*.fasta'):
        target=train if p.name.startswith('discovery.') else test
        source=p.parent.parent.name;species=next(r['species'] for r in masters if r['source']==source)
        for r in readfa(p):target[species].add(r['sequence'])
    leaks=[]
    for species,ss in train.items():
        posting=defaultdict(set)
        for t in test[species]:
            for k in kmers(t):posting[k].add(t)
        for s in ss:
            candidates=set()
            for k in kmers(s):candidates.update(posting[k])
            for t in candidates:
                if related(s,t):leaks.append((species,canonical(s),canonical(t)))
    if leaks:raise ValueError(f'cross-boundary leakage including generated controls: {len(leaks)}')
    manifest=json.loads((base/'manifest.json').read_text());manifest['sigma_groups']=summaries;manifest['annotation_candidate_links']=len(annotations);manifest['independent_fasta_boundary_verification']={'includes_generated_controls':True,'related_pairs':0,'training_sequences_including_controls':sum(map(len,train.values())),'test_sequences_including_controls':sum(map(len,test.values()))};manifest['output_sha256']={str(p.relative_to(base)):sha256_file(p) for p in sorted(base.rglob('*')) if p.is_file() and p.name!='manifest.json'};(base/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2)+'\n');print(json.dumps({'sigma_groups':summaries,'annotation_links':len(annotations),'verification':manifest['independent_fasta_boundary_verification']},indent=2))
if __name__=='__main__':run()
