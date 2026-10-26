"""Coordinate, selection and evaluation safeguards for fragment refinements."""
import numpy as np

from . import pipeline as p


def test_protocol_counts_and_new_seeds():
    cfg=p.read(p.PROTOCOL)
    count=len(cfg['scenarios'])*len(cfg['data_seeds'])*len(cfg['algorithm_seeds'])*len(cfg['representations'])
    assert count==cfg['expected_cluster_fits']==48
    assert count*len(cfg['variants'])==cfg['expected_discovery_outputs']==240
    old=p.read(p.PREVIOUS/'protocol.json')
    assert not set(cfg['data_seeds'])&set(old['data_seeds']+old['stress_seeds'])


def test_alignment_recovers_reverse_complement_without_parent_mutation():
    cfg=p.read(p.PROTOCOL);pattern='ACGTTGACATAT';sequences=['C'*20+pattern+'C'*49,'C'*20+p.reverse_complement(pattern)+'C'*49]
    before=list(sequences)
    raw=p.raw_instances(sequences,{0:20,1:20},np.arange(140,dtype=float),12)
    aligned=p.align_instances(sequences,raw,cfg,np.full(4,.25),False)
    assert [r['fragment'] for r in aligned]==[pattern,pattern]
    assert aligned[1]['strand']=='-' and all(r['shift']==0 for r in aligned)
    assert sequences==before


def test_alignment_offsets_stay_in_parent_coordinates_and_are_not_cumulative():
    cfg=p.read(p.PROTOCOL);pattern='ACGTTGACATAT'
    sequences=['C'*20+pattern+'C'*49,'C'*20+pattern+'C'*49]
    distance=np.ones(140);distance[20]=0
    raw=p.raw_instances(sequences,{0:20,1:18},distance,12)
    aligned=p.align_instances(sequences,raw,cfg,np.full(4,.25),True)
    assert aligned[1]['start0']==20 and aligned[1]['shift']==2
    assert all(abs(r['start0']-r['original_start0'])<=2 for r in aligned)
    for row in aligned:
        actual=sequences[row['sequence_index']][row['start0']:row['end0']]
        if row['strand']=='-':actual=p.reverse_complement(actual)
        assert actual==row['fragment']


def test_edge_options_never_wrap_or_shorten_fragment():
    cfg=p.read(p.PROTOCOL);sequences=['ACGTTGACATAT'+'C'*69,'C'*69+'ACGTTGACATAT']
    raw=p.raw_instances(sequences,{0:0,1:69},np.ones(140),12)
    aligned=p.align_instances(sequences,raw,cfg,np.full(4,.25),True)
    assert all(0<=r['start0']<=69 and len(r['fragment'])==12 for r in aligned)


def test_rejection_is_strict_and_does_not_promote_background():
    cfg=p.read(p.PROTOCOL);_,pwm=p.make_pwm(['A'*12]*20)
    rows=[{'fragment':'A'*12},{'fragment':'C'*12}]
    kept,threshold=p.reject_background(rows,pwm,p.sequence_windows(['C'*81]*8,12),cfg,np.full(4,.25))
    assert kept==[rows[0]]
    kept,_=p.reject_background([rows[0]],pwm,p.sequence_windows(['A'*81]*8,12),cfg,np.full(4,.25))
    assert kept==[]  # Equal to the threshold is rejected.


def test_fast_scanner_agrees_with_previous_scanner():
    from ..bottleneck_20261025.pilot import scan
    sequences=['C'*20+'TTGACATATAAT'+'C'*49,'C'*20+p.reverse_complement('TTGACATATAAT')+'C'*49]
    _,pwm=p.make_pwm(['TTGACATATAAT']*50)
    expected=scan(sequences,pwm,np.full(4,.25))
    actual=p.scan_windows(p.sequence_windows(sequences,12),pwm,np.full(4,.25))
    assert all(np.allclose(a,b) for a,b in zip(expected,actual))


def test_family_fpr_remains_visible_when_reference_gate_fails():
    cfg=p.read(p.PROTOCOL);cfg['synthetic_pattern']='C'*12
    _,pwm=p.make_pwm(['A'*12]*20)
    data={'scenario':{'mutation':0.},'discovery_control':['ACGT'*20+'A'],
          'calibration_control':['C'*81]*8,'test':['A'*81], 'test_control':['A'*81],
          'test_truth':[{'start0':0,'strand':'+'}]}
    prepared={key:p.sequence_windows(data[key],12) for key in ['calibration_control','test','test_control']}
    metrics,nulls,family=p.evaluate({'unrelated':pwm},data,cfg,'test',prepared)
    assert all(r['f1']==0 and not r['reference_eligible'] for r in metrics)
    assert nulls[0]['raw_blank_fpr']==family['raw_family_blank_fpr']==1.


def test_pairing_raw_and_background_never_adds_instances(tmp_path):
    cfg=p.read(p.PROTOCOL);cfg['minimum_promoters']=2
    sequences=['C'*20+'TTGACATATAAT'+'C'*49]*3
    labels=np.zeros(210,dtype=int);distance=np.ones(210)
    for row in range(3):distance[row*70+20]=0
    p.refine_clusters(sequences,['C'*81]*8,labels,distance,cfg,tmp_path)
    import pandas as pd
    for a,b in [('RAW','BG_ONLY'),('RC_SHIFT','RC_SHIFT_BG')]:
        raw=pd.read_csv(tmp_path/a/'instances.tsv',sep='\t')
        kept=pd.read_csv(tmp_path/b/'instances.tsv',sep='\t')
        assert set(kept.sequence_index)<=set(raw.sequence_index)
        assert set(kept.start0)<=set(raw.start0)
