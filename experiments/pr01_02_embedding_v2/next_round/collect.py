import json
import itertools
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from ..audit import BRANCH,ROOT,DATA,load,sha,save
from ..discovery import read_meme,bh,similarity
from .common import NEXT,SOURCES,PRESENTATION,start,finish


def main():
    inputs=list((NEXT/'real_data_comparison').glob('*/development_enrichment.tsv'))+[BRANCH/'summary/all_development_enrichment.tsv']
    target,manifest=start('second_round_collection',{'primary_seed':20261005,'family':'new exploratory family, original q preserved'},inputs)
    original=pd.read_csv(inputs[-1],sep='\t')
    original=original[(original.layer==6)&(original.width==10)].drop_duplicates(['source','seed','method','motif','control'])
    original=original.rename(columns={'q_six_species':'q_original'})
    original['method']=original.method.replace({'prokbert':'A0','onehot':'OH','tokenwindow':'A1'})
    records=[original]
    for path in inputs[:-1]:
        source,method,seed=path.parent.name.split('__');seed=int(seed.removeprefix('s'))
        d=pd.read_csv(path,sep='\t');d['source']=source;d['seed']=seed;records.append(d)
    d=pd.concat(records,ignore_index=True)
    for seed in [20261005,20261006,20261007]:
        for control in ['natural','shuffle']:
            mask=(d.seed==seed)&(d.control==control)&(d.method!='A0_recalc')
            d.loc[mask,'q_new_exploratory']=bh(d.loc[mask,'p'].values)
    out=NEXT/'real_data_comparison'
    d.to_csv(out/'all_methods_enrichment.tsv',sep='\t',index=False)
    base=d[(d.seed==20261005)&(d.method!='A0_recalc')]
    counts=base.groupby(['source','method','control']).agg(candidates=('motif','size'),
            significant_new_family=('q_new_exploratory',lambda x:int((x<=.05).sum()))).reset_index()
    counts.to_csv(out/'primary_seed_counts.tsv',sep='\t',index=False)
    # A2 repeated-seed full-length PWM stability, using the same one-to-one policy.
    stability=[]
    for source in SOURCES:
        paths=[out/f'{source}__A2__s{s}' for s in [20261005,20261006,20261007]]
        for a,b in itertools.combinations(paths,2):
            ma,mb=read_meme(a/'motifs.meme'),read_meme(b/'motifs.meme');na,nb=list(ma),list(mb)
            matrix=np.array([[similarity(ma[x],mb[y],10)[0] for y in nb] for x in na])
            ia,ib=linear_sum_assignment(-matrix)
            da=pd.read_csv(a/'instances.tsv',sep='\t');db=pd.read_csv(b/'instances.tsv',sep='\t')
            for i,j in zip(ia,ib):
                x=set(zip(da[da.motif==na[i]].sequence_id,da[da.motif==na[i]].start0));y=set(zip(db[db.motif==nb[j]].sequence_id,db[db.motif==nb[j]].start0))
                stability.append({'source':source,'run_a':a.name,'run_b':b.name,'motif_a':na[i],'motif_b':nb[j],
                        'full_pwm_similarity':matrix[i,j],'instance_exact_jaccard':len(x&y)/len(x|y)})
    pd.DataFrame(stability).to_csv(out/'A2_seed_stability.tsv',sep='\t',index=False)
    # Algebraically equivalent A0 recalculation is a numerical implementation control.
    numerical=[]
    for path in out.glob('*__A0_recalc__*'):
        source=path.name.split('__')[0]
        a=read_meme(BRANCH/'runs'/f'{source}__prokbert__L6W10__s20261005'/'motifs.meme');b=read_meme(path/'motifs.meme')
        matrix=np.array([[similarity(x,y,10)[0] for y in b.values()] for x in a.values()])
        ia,ib=linear_sum_assignment(-matrix)
        numerical.append({'source':source,'old_candidates':len(a),'recalculated_candidates':len(b),
                    'median_matched_pwm_similarity':float(np.median(matrix[ia,ib])),
                    'interpretation':'small floating point changes can alter clustering; no biological improvement claim'})
    pd.DataFrame(numerical).to_csv(out/'A0_numerical_control.tsv',sep='\t',index=False)
    # Traditional fixed PWM grammar, Development-only custom max-score scan.
    positions,pairs=[] ,[]
    for source in SOURCES:
        scan=pd.read_csv(BRANCH/'runs'/f'{source}__onehot__L6W10__s20261005'/'development_scans.tsv',sep='\t')
        h=scan[(scan.control=='positive')&scan.present&scan.motif.str.startswith('STREME_')].copy();h['source']=source
        positions.append(h)
        for seqid,g in h.groupby('sequence_id'):
            for (_,a),(_,b) in itertools.combinations(g.sort_values(['start1','motif']).iterrows(),2):
                pairs.append({'source':source,'sequence_id':seqid,'upstream':a.motif,'downstream':b.motif,
                             'gap':int(b.start1-a.stop1-1),'center_distance':b.center_tss-a.center_tss,
                             'upstream_strand':a.strand,'downstream_strand':b.strand})
    grammar=NEXT/'figures'
    pd.concat(positions).to_csv(grammar/'traditional_development_positions.tsv',sep='\t',index=False)
    pd.DataFrame(pairs).to_csv(grammar/'traditional_development_pairs.tsv',sep='\t',index=False)
    inventory=[]
    for path in (NEXT/'run_manifests').glob('*.json'):
        m=json.loads(path.read_text());inventory.append({'run_id':path.stem,'status':m['status'],'seconds':m.get('elapsed_seconds'),
              'process_peak_memory_bytes':m.get('process_peak_memory_bytes'),'manifest_sha256':sha(path)})
    pd.DataFrame(inventory).to_csv(NEXT/'run_manifests/run_inventory.tsv',sep='\t',index=False)
    finish(target,manifest,list(out.glob('*.tsv'))+[grammar/'traditional_development_positions.tsv',grammar/'traditional_development_pairs.tsv'])
    print(counts[counts.method=='A2'].to_string(index=False))


if __name__=='__main__':main()
