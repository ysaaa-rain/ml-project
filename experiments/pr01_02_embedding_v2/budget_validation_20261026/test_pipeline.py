"""Controls for candidate selection, equal budgets and replicate-level inference."""
import copy

import numpy as np

from . import pipeline as p


def ranked(names):
    return [{'candidate':name} for name in names]


def test_all_views_and_counts_are_predeclared():
    cfg=p.read(p.PROTOCOL)
    cases=len(cfg['scenarios'])*len(cfg['data_seeds'])*len(cfg['algorithm_seeds'])
    assert cases*len(cfg['representations'])==cfg['expected_cluster_fits']==135
    assert cases*len(cfg['representations'])*len(cfg['variants'])==cfg['expected_pwm_pools']==270
    assert cases*10*len(cfg['variants'])==cfg['expected_evaluation_sets']==900
    old=p.read(p.PREVIOUS/'protocol.json')
    assert not set(cfg['data_seeds'])&set(old['data_seeds'])


def test_matching_limits_actual_counts_and_keeps_rank_order():
    rankings={'EMB128':ranked(['E3','E1','E2','E0','E4','E5']),
              'OH24':ranked(['A1','A0','A2']),'OH48':ranked(['B2','B1','B0','B3','B4','B5'])}
    before=copy.deepcopy(rankings);views=p.memberships(rankings,5)
    assert len(views)==10 and rankings==before
    pair=[v for v in views if v['view']=='MATCH_OH24']
    assert all(v['budget']==3 for v in pair)
    assert pair[0]['candidates']==['E3','E1','E2']
    pair=[v for v in views if v['view']=='MATCH_OH48']
    assert all(v['budget']==5 for v in pair)
    assert pair[1]['candidates']==['B2','B1','B0','B3','B4']


def test_zero_match_does_not_erase_uncapped_or_fixed_cap_baseline():
    views=p.memberships({'EMB128':[],'OH24':ranked(['A']),'OH48':ranked(['B','C'])},5)
    pairs=[v for v in views if v['view'].startswith('MATCH_')]
    assert len(pairs)==4 and all(v['zero_budget'] and not v['candidates'] for v in pairs)
    fixed=[v for v in views if v['view']=='CAP5' and v['representation']=='OH48'][0]
    full=[v for v in views if v['view']=='FULL' and v['representation']=='OH48'][0]
    assert fixed['candidates']==full['candidates']==['B','C']
    assert fixed['underfilled']


def test_discovery_rank_uses_contrast_and_deterministic_ties():
    cfg=p.read(p.PROTOCOL);pattern='TTGACATATAAT'
    _,good=p.make_pwm([pattern]*30);_,bad=p.make_pwm(['C'*12]*30)
    positive=['C'*20+pattern+'C'*49]*16;blank=['C'*81]*16
    rows=p.rank_candidates({'b':good,'a':good,'bad':bad},{'a':30,'b':30,'bad':30},positive,blank,cfg)
    assert [r['candidate'] for r in rows]==['a','b','bad']
    assert rows[0]['score']==1.
    assert all(r['score']==r['positive_presence']-r['blank_presence'] for r in rows)


def test_reference_oracle_cannot_restore_candidate_excluded_by_budget():
    cfg=p.read(p.PROTOCOL)
    scans={'present':np.array([True]),'start':np.array([20]),'forward':np.array([True]),'blank':np.array([False])}
    cache={'selected':{'cosine':.4,'orientation':'+','test':scans},
           'excluded_reference_match':{'cosine':1.,'orientation':'+','test':scans}}
    data={'test':['C'*81],'test_control':['C'*81],'test_truth':[{'start0':20,'strand':'+'}]}
    rows=p.evaluate_view(cache,['selected'],data,cfg,'test')
    assert all(r['f1']==0 and r['candidate']=='selected' for r in rows)
    full=p.evaluate_view(cache,list(cache),data,cfg,'test')
    assert all(r['f1']==1 for r in full)


def test_family_false_positives_survive_reference_gate_failure():
    cfg=p.read(p.PROTOCOL)
    scans={'present':np.array([True]),'start':np.array([20]),'forward':np.array([True]),'blank':np.array([True])}
    cache={'bad':{'cosine':.2,'orientation':'+','test':scans}}
    data={'test':['C'*81],'test_control':['C'*81],'test_truth':[{'start0':20,'strand':'+'}]}
    rows=p.evaluate_view(cache,['bad'],data,cfg,'test')
    assert all(r['f1']==0 and r['family_blank_fpr']==1 for r in rows)
    empty=p.evaluate_view(cache,[],data,cfg,'test')
    assert all(r['f1']==0 and r['family_blank_fpr']==0 for r in empty)


def test_cache_evaluation_matches_previous_full_pool_evaluation():
    from ..refinement_20261026.pipeline import evaluate
    cfg=p.read(p.PROTOCOL);pattern=cfg['synthetic_pattern'];_,pwm=p.make_pwm([pattern]*30)
    data={'scenario':{'mutation':0.},'discovery_control':['ACGT'*20+'A'],
        'calibration_control':['C'*81]*8,'development':['C'*20+pattern+'C'*49],
        'development_control':['C'*81],'development_truth':[{'start0':20,'strand':'+'}],
        'test':['C'*20+pattern+'C'*49],'test_control':['C'*81],
        'test_truth':[{'start0':20,'strand':'+'}]}
    prepared={s:p.sequence_windows(data[s],12) for s in ['calibration_control','development','development_control','test','test_control']}
    cache,_=p.candidate_cache({'M':pwm},data,cfg,prepared)
    actual=p.evaluate_view(cache,['M'],data,cfg,'test')
    expected,_,_=evaluate({'M':pwm},data,cfg,'test',prepared)
    assert [r['f1'] for r in actual]==[r['f1'] for r in expected]
    assert [r['tp'] for r in actual]==[r['tp'] for r in expected]


def test_five_seeds_sign_flip_cannot_claim_p_below_point05():
    assert p.sign_flip(np.ones(5))==.0625
    assert p.sign_flip(np.zeros(5))==1.
    assert p.sign_flip(np.array([1.,0,0,0,0]))==1.


def test_complete_budget_scoring_retains_views_and_primary_pair(tmp_path,monkeypatch):
    import time
    import pandas as pd
    cfg=p.read(p.PROTOCOL);cfg['data_seeds']=[7];cfg['algorithm_seeds']=[9]
    protocol=tmp_path/'protocol.json';p.save(protocol,cfg)
    artifact=tmp_path/'artifacts';summary=tmp_path/'summary';summary.mkdir()
    monkeypatch.setattr(p,'PROTOCOL',protocol);monkeypatch.setattr(p,'ART',artifact)
    monkeypatch.setattr(p,'SUMMARY',summary);monkeypatch.setattr(p,'assert_locked',lambda:None)
    tag='mut30_random__s7';folder=artifact/tag;folder.mkdir(parents=True)
    pattern=cfg['synthetic_pattern'];positive=['C'*20+pattern+'C'*49]*4
    data={'scenario':{'id':'mut30_random','mutation':.3},'seed':7,
          'discovery_control':['ACGT'*20+'A']*4,'calibration_control':['C'*81]*4,
          'development':positive,'development_control':['C'*81]*4,
          'test':positive,'test_control':['C'*81]*4}
    for split in ['development','test']:
        data[split+'_truth']=[{'start0':20,'strand':'+'} for _ in range(4)]
    p.save(folder/'data.json',data);_,pwm=p.make_pwm([pattern]*30)
    rankings={r['id']:ranked(['C00']) for r in cfg['representations']}
    for representation in cfg['representations']:
        for variant in cfg['variants']:
            out=folder/f'{representation["id"]}__a9'/variant;out.mkdir(parents=True)
            p.write_meme(out/'motifs.meme',{'C00':pwm},np.full(4,.25))
    for variant in cfg['variants']:
        p.save(folder/f'views__a9__{variant}.json',p.memberships(rankings,5))
    p.save(summary/'manifest.json',{'status':'FITS_COMPLETED','started_unix':time.time(),'conditions':[tag]})
    p.score()
    rows=pd.read_csv(summary/'metrics.tsv',sep='\t')
    assert len(rows)==160 and set(rows['view'])=={'FULL','CAP5','MATCH_OH24','MATCH_OH48'}
    decision=p.read(summary/'primary_decision.json')
    assert decision['data_seeds']==1 and decision['mean_delta']==0. and not decision['advance_gate_passed']
