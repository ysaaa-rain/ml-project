"""Reuse first-round artifacts; do not refit or alter their statistical families."""
import json
import itertools
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from statsmodels.stats.contingency_tables import StratifiedTable
from ..audit import BRANCH,DATA,load,sha
from ..discovery import read_meme,similarity,bh
from .common import NEXT,SOURCES,start,finish


def stability():
    files=list((BRANCH/'runs').glob('*/manifest.json'))
    target,manifest=start('historical_stability',{'null_repetitions':200,'seed':20261009,'tolerance_bp':[0,1,2,3]},files)
    rng=np.random.default_rng(20261009)
    rows,families,offsets=[] ,[],[]
    for source in SOURCES:
        for method in ['onehot','prokbert','tokenwindow']:
            paths=[BRANCH/'runs'/f'{source}__{method}__L6W10__s{s}' for s in [20261005,20261006,20261007]]
            if not all(p.exists() for p in paths):continue
            matrices=[read_meme(p/'motifs.meme') for p in paths]
            instances=[pd.read_csv(p/'instances.tsv',sep='\t') for p in paths]
            scanned=[pd.read_csv(p/'development_scans.tsv',sep='\t') for p in paths]
            for first,second in itertools.combinations(range(3),2):
                a,b=matrices[first],matrices[second];an,bn=list(a),list(b)
                for strategy in ['full','local']:
                    comparisons=[[similarity(a[x],b[y],minimum_overlap=10 if strategy=='full' else 6) for y in bn] for x in an]
                    matrix=np.array([[s[0] for s in r] for r in comparisons])
                    ai,bi=linear_sum_assignment(-matrix)
                    for i,j in zip(ai,bi):
                        score,shift,orientation=comparisons[i][j]
                        namea,nameb=an[i],bn[j]
                        families.append({'source':source,'method':method,'seed_a':20261005+first,'seed_b':20261005+second,
                          'strategy':strategy,'motif_a':namea,'motif_b':nameb,'similarity':score,'offset':shift,'orientation':orientation,
                          'length_a':len(a[namea]),'length_b':len(b[nameb]),
                          **{f'passes_{t}':score>=t for t in [.6,.7,.8,.9]}})
                        for stage in ['discovery_instances','development_scans']:
                            if stage=='discovery_instances':
                                da=instances[first][instances[first].motif==namea].copy()
                                db=instances[second][instances[second].motif==nameb].copy()
                                da['pos']=da.start0;db['pos']=db.start0
                            else:
                                da=scanned[first][(scanned[first].motif==namea)&(scanned[first].control=='positive')&scanned[first].present].copy()
                                db=scanned[second][(scanned[second].motif==nameb)&(scanned[second].control=='positive')&scanned[second].present].copy()
                                da['pos']=da.start1-1;db['pos']=db.start1-1
                            assert not da.sequence_id.duplicated().any() and not db.sequence_id.duplicated().any()
                            joined=da[['sequence_id','pos','strand']].merge(db[['sequence_id','pos','strand']],on='sequence_id',suffixes=('_a','_b'))
                            if not len(joined):continue
                            canonical_b=joined.strand_b if orientation=='+' else joined.strand_b.map({'+':'-','-':'+'})
                            equal_strand=(joined.strand_a==canonical_b).to_numpy()
                            xa,xb=joined.pos_a.to_numpy(),joined.pos_b.to_numpy()
                            aligned_b=xb+np.where(canonical_b.to_numpy()=='+',-shift,shift)
                            delta=aligned_b-xa
                            raw_delta=xb-xa
                            total=len(da)+len(db)
                            near_tss=abs(xa+4.5-60)<=15
                            for tol in [0,1,2,3]:
                                match=int(((abs(delta)<=tol)&equal_strand).sum())
                                physical=int((abs(raw_delta)<=tol).sum())
                                random_starts=rng.integers(0,72,(200,len(joined)))
                                random_matches=((abs(random_starts-xa)<=tol)&equal_strand).sum(axis=1)
                                null=random_matches/(total-random_matches)
                                empirical_matches=np.array([int(((abs(rng.permutation(xb)-xa)<=tol)&equal_strand).sum()) for _ in range(200)])
                                rows.append({'source':source,'method':method,'strategy':strategy,'stage':stage,
                                  'seed_a':20261005+first,'seed_b':20261005+second,'motif_a':namea,'motif_b':nameb,
                                  'pwm_similarity':score,'alignment_offset':shift,'alignment_orientation':orientation,
                                  'tolerance_bp':tol,'n_a':len(da),'n_b':len(db),'common_sequence_ids':len(joined),
                                  'direction_aware_aligned_jaccard':match/(total-match),
                                  'physical_position_jaccard':physical/(total-physical),
                                  'uniform_null_jaccard_mean':float(null.mean()),'uniform_null_jaccard_p95':float(np.quantile(null,.95)),
                                  'empirical_position_null_mean':float(np.mean(empirical_matches/(total-empirical_matches))),
                                  'near_tss_match_fraction':float(((abs(delta[near_tss])<=tol)&equal_strand[near_tss]).mean()) if near_tss.any() else np.nan,
                                  'other_match_fraction':float(((abs(delta[~near_tss])<=tol)&equal_strand[~near_tss]).mean()) if (~near_tss).any() else np.nan,
                                  'mean_interval_iou':float(np.mean(np.maximum(0,10-abs(delta))/(20-np.maximum(0,10-abs(delta))))),
                                  'threshold_a':json.loads((paths[first]/'thresholds.json').read_text())['thresholds'][namea],
                                  'threshold_b':json.loads((paths[second]/'thresholds.json').read_text())['thresholds'][nameb]})
                            if strategy=='full' and score>=.8:
                                offsets.extend({'source':source,'method':method,'stage':stage,'motif_a':namea,'motif_b':nameb,
                                   'seed_a':20261005+first,'seed_b':20261005+second,'sequence_id':r.sequence_id,
                                   'delta_bp':int(d),'canonical_strand_match':bool(same),'start_a0':int(r.pos_a),'start_b0':int(r.pos_b)}
                                   for r,d,same in zip(joined.itertuples(),delta,equal_strand))
    out=NEXT/'stability_diagnostics'
    pd.DataFrame(rows).to_csv(out/'tolerance_stability.tsv',sep='\t',index=False)
    pd.DataFrame(families).to_csv(out/'pwm_family_matches.tsv',sep='\t',index=False)
    pd.DataFrame(offsets).to_csv(out/'site_offsets.tsv',sep='\t',index=False)
    finish(target,manifest,list(out.glob('*.tsv')))
    print('Historical stability diagnostics completed',flush=True)


def gc():
    inputs=[DATA/s/'main/sequence_master.tsv' for s in SOURCES]
    target,manifest=start('historical_gc_and_shuffle',{'gc_bins':[0,.2,.4,.6,.8,1],'cmh_min_per_class':5},inputs)
    gc_rows,overlap_rows,motif_rows,strata_rows,shuffle_rows=[],[],[],[],[]
    baseline=pd.read_csv(BRANCH/'summary/all_development_enrichment.tsv',sep='\t')
    base=baseline[(baseline.layer==6)&(baseline.width==10)&(baseline.seed==20261005)].drop_duplicates(['source','method','motif','control'])
    m=int((base.control=='shuffle').sum())
    for source in SOURCES:
        for split in ['discovery','development']:
            pos,neg=load(source,split,1),load(source,split,0)
            a=pos.sequence.str.count('[GC]')/81;b=neg.sequence.str.count('[GC]')/81
            ha=np.histogram(a,bins=np.linspace(0,1,21))[0]/len(a);hb=np.histogram(b,bins=np.linspace(0,1,21))[0]/len(b)
            overlap_rows.append({'source':source,'split':split,'positive_n':len(a),'natural_n':len(b),
                       'positive_gc_mean':a.mean(),'natural_gc_mean':b.mean(),'difference_percentage_points':100*(a.mean()-b.mean()),
                       'histogram_overlap_20_bins':float(np.minimum(ha,hb).sum())})
            gc_rows += [{'source':source,'split':split,'label':label,'sequence_id':r.sequence_id,'leakage_group':r.leakage_group,
                         'gc':r.sequence.count('G')/81+r.sequence.count('C')/81} for label,frame in [(1,pos),(0,neg)] for r in frame.itertuples()]
        frames={control:load(source,'development',label) for control,label in [('positive',1),('natural',0)]}
        seen_streme=False
        for method in ['onehot','prokbert','tokenwindow']:
            path=BRANCH/'runs'/f'{source}__{method}__L6W10__s20261005'
            scan=pd.read_csv(path/'development_scans.tsv',sep='\t');matrices=read_meme(path/'motifs.meme')
            from ..discovery import read_meme as read
            matrices.update({f'STREME_{k}':v for k,v in read(BRANCH.parents[1]/'results/motif/tju_streme_v2_20261008'/source/'selected_main.meme').items()})
            for name,pwm in matrices.items():
                is_streme=name.startswith('STREME_')
                if is_streme and seen_streme:continue
                candidate_method='STREME' if is_streme else method
                groups={}
                for control,frame in frames.items():
                    hits=scan[(scan.motif==name)&(scan.control==control)]
                    d=frame[['sequence_id','leakage_group','sequence']].merge(hits[['sequence_id','present']],on='sequence_id')
                    d['gc']=d.sequence.str.count('[GC]')/81
                    groups[control]=d.groupby('leakage_group').agg(gc=('gc','mean'),present=('present','max'))
                shared=set(groups['positive'].index)&set(groups['natural'].index)
                for control in groups:groups[control]=groups[control].drop(list(shared))
                tables=[];covered={c:0 for c in groups}
                for left,right in zip([0,.2,.4,.6,.8],[.2,.4,.6,.8,1]):
                    a=groups['positive'];b=groups['natural']
                    a=a[(a.gc>=left)&(a.gc<right if right<1 else a.gc<=right)]
                    b=b[(b.gc>=left)&(b.gc<right if right<1 else b.gc<=right)]
                    table=np.array([[a.present.sum(),len(a)-a.present.sum()],[b.present.sum(),len(b)-b.present.sum()]],float)
                    eligible=len(a)>=5 and len(b)>=5
                    strata_rows.append({'source':source,'method':candidate_method,'motif':name,'gc_left':left,'gc_right':right,
                       'positive_groups':len(a),'natural_groups':len(b),'positive_hit_rate':a.present.mean(),'natural_hit_rate':b.present.mean(),
                       'cmh_eligible':eligible,'risk_difference':a.present.mean()-b.present.mean()})
                    if eligible:tables.append(table);covered['positive']+=len(a);covered['natural']+=len(b)
                cmh_or,cmh_p=np.nan,np.nan
                if tables:
                    try:
                        import warnings
                        with warnings.catch_warnings():
                            warnings.simplefilter('ignore')
                            st=StratifiedTable(np.stack(tables,axis=2),shift_zeros=False)
                            cmh_or=float(st.oddsratio_pooled);cmh_p=float(st.test_null_odds().pvalue)
                    except ValueError:pass
                selected=base[(base.source==source)&(base.method==candidate_method)&(base.motif==name)]
                natural=selected[selected.control=='natural'].iloc[0];shuf=selected[selected.control=='shuffle'].iloc[0]
                motif_rows.append({'source':source,'method':candidate_method,'motif':name,'pwm_gc':float(pwm[:,1:3].sum(axis=1).mean()),
                  'natural_risk_difference':natural.positive_hits/natural.positive_groups-natural.control_hits/natural.control_groups,
                  'original_natural_q':natural.q_six_species,'cmh_or_exploratory':cmh_or,'cmh_p_exploratory':cmh_p,
                  'cmh_positive_coverage':covered['positive']/len(groups['positive']),
                  'cmh_natural_coverage':covered['natural']/len(groups['natural'])})
                d=int(shuf.discordant_pos_only+shuf.discordant_null_only)
                shuffle_rows.append({'source':source,'method':candidate_method,'motif':name,'paired_groups':int(shuf.positive_groups),
                     'pos_only':int(shuf.discordant_pos_only),'null_only':int(shuf.discordant_null_only),'discordant_groups':d,
                     'paired_risk_difference':(shuf.discordant_pos_only-shuf.discordant_null_only)/shuf.positive_groups,
                     'raw_p':shuf.p,'original_q':shuf.q_six_species,'family_m':m,
                     'minimum_attainable_one_sided_p':2.**(-d) if d else 1.,'first_rank_bh_cutoff':.05/m,
                     'discordance_too_small_for_first_rank_bh':2.**(-d)>.05/m})
            seen_streme=True
    out=NEXT/'gc_diagnostics'
    pd.DataFrame(gc_rows).to_csv(out/'gc_per_sequence.tsv',sep='\t',index=False)
    pd.DataFrame(overlap_rows).to_csv(out/'gc_distribution_overlap.tsv',sep='\t',index=False)
    d=pd.DataFrame(motif_rows);valid=d.cmh_p_exploratory.notna();d.loc[valid,'cmh_q_exploratory']=bh(d.loc[valid,'cmh_p_exploratory'])
    d.to_csv(out/'motif_gc_adjustment.tsv',sep='\t',index=False)
    pd.DataFrame(strata_rows).to_csv(out/'gc_strata_effects.tsv',sep='\t',index=False)
    pd.DataFrame(shuffle_rows).to_csv(out/'shuffle_discordance_diagnostic.tsv',sep='\t',index=False)
    finish(target,manifest,list(out.glob('*.tsv')),{'interpretation':'CMH is additional exploration, original q family preserved; attainable-p bounds are not post-hoc power estimates.'})
    print('GC and paired-shuffle diagnostics completed',flush=True)


if __name__=='__main__':stability();gc()
