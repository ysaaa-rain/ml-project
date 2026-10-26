"""Verify frozen experiment, reconstructing ranking PWMs at their original precision.

The preregistered verifier compares floating-point ranking records recomputed from
8-decimal MEME text with records computed before that serialization. This audit
keeps all scientific sources, models, rankings and scoring unchanged. It rebuilds
the original full-precision PWMs from saved instances and verifies exact ranks.
"""
import joblib
import numpy as np
import pandas as pd
import torch

from .pipeline import (ART,HERE,METHODS,PROTOCOL,ROOT,SUMMARY,LocalCNN,assert_locked,make_pwm,
    maxima,rank_candidates,read,read_meme,reverse_complement,save,setup,sha,views)


def main():
    cfg=read(PROTOCOL);setup(cfg);assert_locked();manifest=read(SUMMARY/'manifest.json')
    assert manifest['status']=='COMPLETED'
    for field in ['cnn_trains','random_models','cluster_fits','pwm_pools','evaluation_sets']:
        assert manifest[field]==cfg['expected_'+field]
    protected=read(SUMMARY/'protected_before.json');assert all(sha(ROOT/n)==d for n,d in protected.items())
    instances=0;pwms=0;rankings_checked=0;views_checked=0;truths=0;cnn_maps=0
    precision_deltas={};candidate_order_changes=0;max_pwm_error=0.
    for folder in sorted(ART.iterdir()):
        data=read(folder/'data.json');seen=set()
        for split in ['discovery','discovery_control','development','development_control','test','test_control','calibration_control']:
            canonical={min(s,reverse_complement(s)) for s in data[split]}
            assert len(canonical)==len(data[split]) and not canonical&seen;seen|=canonical
        for split in ['discovery','development','test']:
            for r in data[split+'_truth']:
                assert data[split][r['sequence_index']][r['start0']:r['end0']]==r['fragment'];truths+=1
        for parent in sorted(folder.glob('a*')):
            rankings={}
            for method in METHODS:
                out=parent/method;motifs=read_meme(out/'motifs.meme');frame=pd.read_csv(out/'instances.tsv',sep='\t')
                assert set(frame.motif)==set(motifs) and not frame.duplicated(['motif','sequence_index']).any()
                for r in frame.itertuples():
                    dna=data['discovery'][r.sequence_index][r.start0:r.end0]
                    assert r.end0-r.start0==cfg['width'] and 0<=r.start0<=81-cfg['width']
                    assert r.fragment==(dna if r.strand=='+' else reverse_complement(dna));instances+=1
                original={}
                for name,pwm in motifs.items():
                    rows=frame[frame.motif==name];assert len(rows)>=cfg['minimum_promoters']
                    _,original[name]=make_pwm(rows.fragment.tolist())
                    error=float(np.abs(original[name]-pwm).max());max_pwm_error=max(max_pwm_error,error)
                    assert error<=5.000001e-9;pwms+=1
                support=frame.groupby('motif').size().to_dict()
                rankings[method]=rank_candidates(original,support,data['discovery'],data['discovery_control'],cfg)
                saved=read(out/'ranking.json');assert rankings[method]==saved;rankings_checked+=1
                rounded=rank_candidates(motifs,support,data['discovery'],data['discovery_control'],cfg)
                candidate_order_changes+=int([r['candidate'] for r in rounded]!=[r['candidate'] for r in saved])
                saved_by_id={r['candidate']:r for r in saved}
                for row in rounded:
                    before=saved_by_id[row['candidate']]
                    for key,value in row.items():
                        if isinstance(value,(int,float)):
                            precision_deltas[key]=max(precision_deltas.get(key,0.),abs(value-before[key]))
                if method.startswith('CNN'):
                    model=LocalCNN(cfg['cnn']['filters'],cfg['width']);model.load_state_dict(torch.load(out/'model.pt',weights_only=True));model.eval()
                    score,start,forward=maxima(model,data['discovery']);null,_,_=maxima(model,data['discovery_control'])
                    saved_maps=np.load(out/'activation_maxima.npz');threshold=np.quantile(null,.95,axis=0,method='higher')
                    assert np.array_equal(score,saved_maps['score']) and np.array_equal(start,saved_maps['start']) and np.array_equal(forward,saved_maps['forward'])
                    assert np.array_equal(threshold,saved_maps['threshold'])
                    expected={f'F{j:02d}' for j in range(score.shape[1]) if (score[:,j]>threshold[j]).sum()>=cfg['minimum_promoters']}
                    assert expected==set(motifs)
                    for name in expected:
                        j=int(name[1:]);selected=frame[frame.motif==name]
                        assert set(selected.sequence_index)==set(np.flatnonzero(score[:,j]>threshold[j]))
                        for r in selected.itertuples():
                            assert r.start0==start[r.sequence_index,j] and (r.strand=='+')==forward[r.sequence_index,j]
                    cnn_maps+=1
                else:assert joblib.load(out/'kmeans.joblib').n_features_in_==48
            actual=read(parent/'views.json');assert actual==views(rankings,cfg['budget_cap']);views_checked+=len(actual)
    assert views_checked==cfg['expected_evaluation_sets'] and candidate_order_changes==0
    precision={'status':'EXPLAINED','original_verifier':'pipeline.verify exact equality from rounded PWM fails',
        'cause':'ranking originally uses full-precision PWM; MEME text stores8 decimal places, with strict threshold ties',
        'original_rankings_exactly_rebuilt_from_saved_instances':rankings_checked,
        'rounded_PWM_candidate_order_changes':candidate_order_changes,'max_PWM_absolute_error':max_pwm_error,
        'max_ranking_field_changes_after_PWM_rounding':precision_deltas,
        'scientific_sources_models_rankings_scores_modified':False,'audit_source_sha256':sha(HERE/'verify_artifacts.py')}
    save(SUMMARY/'precision_audit.json',precision)
    inventory={p.relative_to(HERE).as_posix():sha(p) for p in ART.rglob('*') if p.is_file()}
    inventory.update({p.relative_to(HERE).as_posix():sha(p) for p in SUMMARY.iterdir() if p.is_file() and p.name not in ['artifact_hashes.json','verification.json']})
    save(SUMMARY/'artifact_hashes.json',inventory)
    result={'status':'PASSED','artifact_hashes':len(inventory),'historical_files_unchanged':len(protected),
        'instances_rebuilt':instances,'pwms_rebuilt':pwms,'rankings_recomputed':rankings_checked,
        'cnn_activation_maps_recomputed':cnn_maps,'truth_sites_checked':truths,'evaluation_sets':views_checked,
        'exact_RC_duplicates_per_condition':0,'real_holdout_used':False,
        'verification_source_sha256':sha(HERE/'verify_artifacts.py'),
        'precision_policy':'ranking recomputed from original-precision PWM reconstructed from instances; rounded PWM order unchanged',
        'limits':'saved-model activation recomputation and artifact consistency, not independent retraining or near-homology audit'}
    save(SUMMARY/'verification.json',result);print(result,flush=True)


if __name__=='__main__':main()
