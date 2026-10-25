"""Scientific safeguards for nested reductions, paired filters and held-out ranking."""
from types import SimpleNamespace

import numpy as np

from . import pipeline as p


def test_protocol_counts_and_seed_separation():
    cfg=p.read(p.PROTOCOL)
    main=len(cfg['scenarios'])*len(cfg['data_seeds'])*len(cfg['algorithm_seeds'])*len(cfg['methods'])
    stress=len(cfg['stress_scenarios'])*len(cfg['stress_seeds'])*len(cfg['algorithm_seeds'])*len(cfg['stress_methods'])
    assert main+stress==cfg['expected_discovery_outputs']==324
    assert not set(cfg['data_seeds']+cfg['stress_seeds']) & set(p.read(p.OLD/'protocol.json')['data_seeds'])


def test_median_filter_reuses_assignments_and_cannot_add_instances():
    labels=np.array([0,0,0,0,1,1,1,1]); distance=np.array([.8,.1,.2,.9,.4,.3,.7,.6])
    full,stats=p.select_instances(labels,distance,2,1.0)
    filtered,filtered_stats=p.select_instances(labels,distance,2,.5)
    for label in full:
        assert set(filtered[label].items()) <= set(full[label].items())
    assert filtered[0]=={0:1,1:0}
    assert filtered_stats[0]['windows_after']==2
    assert all(s['windows_after']<=t['windows_after'] for s,t in zip(filtered_stats,stats))


def test_deduplicated_crop_encoding_keeps_every_core(tmp_path,monkeypatch):
    cfg=p.read(p.PROTOCOL); cfg['batch_size']=7
    sequence=''.join(np.random.default_rng(21).choice(list('ACGT'),81))
    calls=[]

    def fake_encode(sequences,model,tokenizer,kmer,batch_size,layer):
        calls.extend(sequences)
        return np.array([[[sum(map(ord,s[i:i+6])),i] for i in range(19)] for s in sequences],dtype=np.float32)
    monkeypatch.setattr(p,'encode',fake_encode)
    features=p.local_order([sequence],SimpleNamespace(config=SimpleNamespace(hidden_size=2)),None,cfg,tmp_path/'test.npy')
    assert len(calls)==58
    for start in range(70):
        left,offset=p.crop_spec(start,12,24)
        expected=np.array([[sum(map(ord,sequence[i:i+6])),i-left] for i in range(start,start+7)],dtype=np.float32).reshape(-1)
        assert np.array_equal(features[0,start],expected)


def test_reverse_truth_direction_is_evaluated_correctly():
    cfg=p.read(p.PROTOCOL); pattern=cfg['synthetic_pattern']; _,pwm=p.make_pwm([pattern]*100)
    data={'scenario':{'mutation':0.},'discovery_control':['ACGT'*20+'A'],
          'calibration_control':['C'*81]*8,'test':['C'*20+p.reverse_complement(pattern)+'C'*49],
          'test_control':['C'*81],'test_truth':[{'start0':20,'strand':'-'}]}
    metrics,_=p.evaluate({'exact':pwm},data,cfg,'test')
    assert all(row['f1']==1 for row in metrics)


def test_wrong_strand_does_not_get_localization_credit():
    cfg=p.read(p.PROTOCOL); pattern=cfg['synthetic_pattern']; _,pwm=p.make_pwm([pattern]*100)
    data={'scenario':{'mutation':0.},'discovery_control':['ACGT'*20+'A'],
          'calibration_control':['C'*81]*8,'test':['C'*20+pattern+'C'*49],
          'test_control':['C'*81],'test_truth':[{'start0':20,'strand':'-'}]}
    metrics,_=p.evaluate({'exact':pwm},data,cfg,'test')
    assert all(row['tp']==0 for row in metrics)


def test_probe_negative_examples_never_overlap_truth():
    cfg=p.read(p.PROTOCOL)
    # Use start index as the synthetic feature so sampled window coordinates remain observable.
    features=np.arange(70,dtype=float).reshape(1,70,1)
    x,y=p.probe_examples(features,[{'sequence_index':0,'start0':25}],cfg)
    assert x[y==1].ravel().tolist()==[25]
    assert all(start+12<=25 or start>=37 for start in x[y==0].ravel())
    assert len(set(x[y==0].ravel().tolist()))==4


def test_new_streams_reproduce_truth_and_do_not_depend_on_training_size():
    cfg=p.read(p.PROTOCOL); scenario=cfg['stress_scenarios'][2]
    data=p.create_data(cfg,scenario,cfg['stress_seeds'][0])
    changed=p.create_data({**cfg,'discovery_n':16},scenario,cfg['stress_seeds'][0])
    assert data['test']==changed['test']
    assert {item['strand'] for item in data['test_truth']}=={'+','-'}
    for split in ['discovery','development','test']:
        for truth in data[split+'_truth']:
            assert data[split][truth['sequence_index']][truth['start0']:truth['end0']]==truth['fragment']
