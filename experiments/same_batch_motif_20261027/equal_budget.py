"""Re-evaluate every method at the same available budget within each paired run."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .pipeline import ROOT, HERE, BASE, ART, config, read, sha, candidate_cache, read_meme, similarity, sequence_windows, ALPHABET, evaluate_partial, save

OUT=ROOT/'reports/same_batch_motif_equal_20261027'

def main():
    p,c=config(); lock=read(HERE/'fit_lock.json')
    assert sha(HERE/'protocol.json')==lock['protocol_sha256'] and sha(HERE/'pipeline.py')==lock['pipeline_sha256']
    for rel,digest in lock['source_sha256'].items():assert sha(ROOT/rel)==digest
    for rel,digest in lock['output_sha256'].items():assert sha(HERE/rel)==digest
    rows=[]
    for scenario in p['scenario_ids']:
        for seed in p['data_seeds']:
            name=f'{scenario}__s{seed}'; old=BASE/'artifacts'/name;folder=ART/name;data=read(old/'data.json')
            prepared={k:sequence_windows(data[k],c['width']) for k in ['calibration_control','development','development_control','test','test_control']}
            reference=np.full((c['width'],4),data['scenario']['mutation']/3)
            for i,char in enumerate(c['synthetic_pattern']):reference[i,ALPHABET.index(char)]=1-data['scenario']['mutation']
            for alg in p['algorithm_seeds']:
                dirs={};rankings={}
                for method in p['methods']:
                    if method.startswith('TIED'):
                        directory=folder/f'a{alg}'/method; ranking=read(directory/'cap5.json')
                    else:
                        directory=old/f'{method}__a{alg}'/'BG_ONLY';ranking=[r['candidate'] for r in read(directory/'ranking.json')[:p['candidate_budget']]]
                    dirs[method]=directory;rankings[method]=ranking
                budget=min(p['candidate_budget'],*(len(v) for v in rankings.values()))
                for method in p['methods']:
                    motifs=read_meme(dirs[method]/'motifs.meme');cache,_=candidate_cache(motifs,data,c,prepared)
                    matches={n:similarity(pwm,reference,minimum_overlap=p['partial_min_overlap_bp']) for n,pwm in motifs.items()}
                    for split in ['development','test']:
                        rows.append({'scenario':scenario,'data_seed':seed,'algorithm_seed':alg,'method':method,'common_budget':budget,
                            **evaluate_partial(cache,matches,rankings[method][:budget],data,split,p['position_tolerance_bp'])})
            print('equal budget scored',name,flush=True)
    OUT.mkdir(parents=True,exist_ok=False);frame=pd.DataFrame(rows);assert len(frame)==450
    assert (frame.candidate_count==frame.common_budget).all()
    frame.to_csv(OUT/'runs.tsv',sep='\t',index=False,lineterminator='\n')
    per=frame.groupby(['scenario','data_seed','method','split'],as_index=False).agg(f1=('f1','mean'),blank_fpr=('family_blank_fpr','mean'),candidate_count=('candidate_count','mean'),overlap_bp=('overlap_bp','mean'))
    per.to_csv(OUT/'per_data_seed.tsv',sep='\t',index=False,lineterminator='\n')
    summary=per.groupby(['scenario','method','split'],as_index=False).agg(f1_mean=('f1','mean'),f1_min=('f1','min'),f1_max=('f1','max'),blank_fpr=('blank_fpr','mean'),n_data=('data_seed','nunique'))
    summary.to_csv(OUT/'summary.tsv',sep='\t',index=False,lineterminator='\n')
    save(OUT/'manifest.json',{'protocol':p,'fit_lock_sha256':sha(HERE/'fit_lock.json'),'equal_budget_script_sha256':sha(Path(__file__)),
        'budget_rule':'min(5, candidate counts of all five methods in each same data/algorithm seed)',
        'budget_counts':frame[frame.method==p['methods'][0]].groupby(['split','common_budget']).size().to_dict().__str__(),
        'run_rows':len(frame),'outputs_sha256':{x.name:sha(x) for x in OUT.iterdir() if x.is_file()}})
    print('EQUAL_SCORE_COMPLETED',flush=True)

if __name__=='__main__':main()
