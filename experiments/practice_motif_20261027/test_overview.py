"""Coordinate projection and candidate-set evaluation sanity checks."""
import numpy as np

from .build import evaluate_partial,project_reference_start


def test_reference_projection_tracks_scan_direction():
    starts=np.array([10,10,10,10]);forward=np.array([True,False,True,False])
    assert project_reference_start(starts,forward,2).tolist()==[12,8,12,8]
    assert project_reference_start(starts,forward,-2).tolist()==[8,12,8,12]


def test_partial_evaluation_uses_frozen_candidates_and_counts_blank_calls():
    data={'test':['A'*20,'C'*20],'test_control':['A'*20,'C'*20],
          'test_truth':[{'start0':7,'strand':'+'},{'start0':12,'strand':'-'}]}
    good={'test':{'present':np.array([True,True]),'blank':np.array([True,False]),
                  'start':np.array([5,14]),'forward':np.array([True,False])}}
    other={'test':{'present':np.array([False,False]),'blank':np.array([False,True]),
                   'start':np.array([0,0]),'forward':np.array([True,True])}}
    matches={'good':(.2,2,'+'),'other':(.9,0,'+')}
    row=evaluate_partial({'good':good,'other':other},matches,['good'],data,'test',0)
    assert row['candidate']=='good' and row['tp']==2 and row['prediction_count']==3
    assert row['f1']==.8 and row['family_blank_fpr']==.5
    assert row['partial_cosine']==.2 and row['overlap_bp']==10


def test_empty_candidate_view_is_retained():
    data={'test':['A'*20],'test_control':['C'*20],'test_truth':[{'start0':3,'strand':'+'}]}
    row=evaluate_partial({}, {}, [], data,'test',3)
    assert row['candidate']=='' and row['f1']==0 and row['prediction_count']==0
    assert row['family_blank_fpr']==0 and row['truth_n']==1
