"""Coordinate, supervision-control, budget and frozen-subset acceptance tests."""
import copy

import numpy as np
import pytest
import torch

from .pipeline import LocalCNN,cnn_pool,evaluate_view,maxima,read_meme,tensor,train,views
from ..bottleneck_20261025.pilot import reverse_complement


def detector(pattern):
    model=LocalCNN(1,len(pattern))
    with torch.no_grad():
        model.conv.weight.zero_();model.conv.bias.zero_()
        for i,c in enumerate(pattern):model.conv.weight[0,'ACGT'.index(c),i]=1
    return model


@pytest.mark.parametrize('sequence,start,strand',[
    ('CCCAACGCCCC',3,True),('CCCGTTCCCCCC',2,False),('AACGCCCCCCC',0,True),('CCCCCCCCCGTT',8,False)])
def test_receptive_field_and_reverse_mapping(sequence,start,strand):
    score,starts,forward=maxima(detector('AACG'),[sequence])
    assert score[0,0]==4 and starts[0,0]==start and forward[0,0]==strand
    fragment=sequence[start:start+4]
    assert (fragment if strand else reverse_complement(fragment))=='AACG'


def test_classification_reverse_complement_invariance():
    torch.manual_seed(1);model=LocalCNN(5,4)
    dna=['AACGCCCCCCC','TTTTGCAAAA']
    # Same-length sequences are required by tensor; batch each separately here.
    for sequence in dna:
        assert torch.allclose(model(tensor([sequence])),model(tensor([reverse_complement(sequence)])))


def test_training_updates_only_trained_copy_and_permutation_control():
    torch.set_num_threads(1);torch.manual_seed(10);initial=LocalCNN(3,2)
    actual=copy.deepcopy(initial);shuffled=copy.deepcopy(initial)
    cfg={'epochs':3,'batch_size':4,'learning_rate':.01,'weight_decay':.0001}
    seq=['AAAAAA']*8+['CCCCCC']*8;labels=[1]*8+[0]*8
    _,truth_hash=train(actual,seq,labels,cfg,7)
    _,shuffle_hash=train(shuffled,seq,labels,cfg,7,True)
    assert truth_hash!=shuffle_hash
    assert any(not torch.equal(a,b) for a,b in zip(actual.parameters(),initial.parameters()))
    assert any(not torch.equal(a,b) for a,b in zip(actual.parameters(),shuffled.parameters()))


@pytest.mark.parametrize('n_cnn,n_oh',[(0,7),(7,0),(2,8),(7,9)])
def test_matching_budget_preserves_empty_and_underfilled(n_cnn,n_oh):
    rankings={m:[{'candidate':str(i)} for i in range(n)] for m,n in
              [('CNN_LABEL',n_cnn),('OH48_BG',n_oh),('CNN_RANDOM',4),('CNN_SHUFFLED',3)]}
    all_views=views(rankings);matched=[v for v in all_views if v['view']=='MATCH_OH48']
    assert len(all_views)==10 and [v['budget'] for v in matched]==[min(5,n_cnn,n_oh)]*2
    assert all(v['underfilled']==(v['budget']<5) and v['zero_budget']==(v['budget']==0) for v in matched)


def test_empty_cnn_pool_is_saved_without_fabricating_candidates(tmp_path):
    cfg={'width':4,'minimum_promoters':2};model=detector('AACG')
    # Identical positives/blanks give equality to threshold, which must be rejected.
    cnn_pool(model,['CCAACGCC']*3,['CCAACGCC']*3,cfg,tmp_path)
    assert read_meme(tmp_path/'motifs.meme')=={}


def test_reference_cannot_recover_candidate_excluded_from_budget():
    cfg={'evaluation':{'reference_cosine_gate':.8,'position_tolerances':[2]}}
    fake={'cosine':1.,'orientation':'+','test':{'blank':np.array([True,False]),
          'present':np.array([True]),'start':np.array([1]),'forward':np.array([True])}}
    data={'test':['AAAA'],'test_control':['CCCC','GGGG'],'test_truth':[{'start0':1,'strand':'+'}]}
    empty=evaluate_view({'excluded':fake},[],data,cfg,'test')[0]
    assert empty['f1']==0 and empty['candidate_count']==0


def test_family_false_positives_are_not_hidden_by_reference_gate():
    cfg={'evaluation':{'reference_cosine_gate':.8,'position_tolerances':[2]}}
    candidate={'cosine':.2,'orientation':'+','test':{'blank':np.array([True,False]),
          'present':np.array([True]),'start':np.array([1]),'forward':np.array([True])}}
    data={'test':['AAAA'],'test_control':['CCCC','GGGG'],'test_truth':[{'start0':1,'strand':'+'}]}
    row=evaluate_view({'bad':candidate},['bad'],data,cfg,'test')[0]
    assert row['f1']==0 and row['family_blank_fpr']==.5
