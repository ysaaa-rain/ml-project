"""Boundary, orientation and fair-budget checks for the new CNN pilot."""
import numpy as np
import pytest
import torch

from .pipeline import (ForwardCNN,LocalCNN,choose_shift,direction_maxima,shifted_fragment,
                       tensor,views)
from ..bottleneck_20261025.pilot import reverse_complement


@pytest.mark.parametrize('forward,shift,start,expected_start',[(True,2,5,7),(True,-2,5,3),(False,2,5,3),(False,-2,5,7)])
def test_shift_uses_filter_orientation(forward,shift,start,expected_start):
    seq='ACGTACGTACGTACGTACGT'
    actual,fragment=shifted_fragment(seq,start,forward,shift,4)
    assert actual==expected_start
    dna=seq[actual:actual+4]
    assert fragment==(dna if forward else reverse_complement(dna))


def test_boundary_intersection_drops_same_sequences_for_all_offsets():
    sequences=['A'*20,'C'*20,'G'*20,'T'*20]
    records=[{'sequence_index':i,'start0':s,'forward':True} for i,s in enumerate([0,2,14,16])]
    eligible,_,scores=choose_shift(sequences,records,np.full(4,.25),4)
    assert [r['sequence_index'] for r in eligible]==[1,2]
    assert set(scores)=={'-2','-1','0','1','2'}


def test_discovery_information_recovers_a_uniform_known_offset():
    rng=np.random.default_rng(42);pattern='ACGTTGCAACGT';sequences=[];records=[]
    for i in range(256):
        dna=rng.choice(list('ACGT'),81);dna[25:37]=list(pattern);sequences.append(''.join(dna))
        records.append({'sequence_index':i,'start0':26,'forward':True})
    eligible,shift,_=choose_shift(sequences,records,np.full(4,.25),12)
    assert len(eligible)==256 and shift==-1


def test_information_shift_requires_no_reference_or_truth_argument():
    # Only the supplied sequence windows, candidate locations and background are exposed.
    import inspect
    assert list(inspect.signature(choose_shift).parameters)==['sequences','instances','bg','width']


def test_tied_and_forward_paths_share_weights_but_forward_uses_one_strand():
    torch.manual_seed(11);tied=LocalCNN(4,4);forward=ForwardCNN(4,4);forward.load_state_dict(tied.state_dict())
    dna=['AACGCCCCCCCCC','TTTTGCAAAAAAA']
    for sequence in dna:
        x=tensor([sequence]);plus,_=forward.activations(x)
        assert torch.allclose(forward(x),forward.head(plus.amax(-1)).squeeze(-1))
        assert torch.allclose(tied(x),tied(tensor([reverse_complement(sequence)])))
        _,_,strand=direction_maxima(forward,[sequence],'FWD');assert strand.all()


@pytest.mark.parametrize('raw,shift,oh',[(0,4,5),(2,6,7),(8,9,1),(5,5,5)])
def test_matched_budgets_include_zero_and_underfill(raw,shift,oh):
    rankings={m:[{'candidate':str(i)} for i in range(n)] for m,n in
        [('TIED_RAW',raw),('TIED_SHIFT',shift),('FWD_RAW',3),('FWD_SHIFT',4),('OH48_BG',oh)]}
    out=views(rankings);assert len(out)==14
    shift_pair=[v for v in out if v['view']=='MATCH_SHIFT'];oh_pair=[v for v in out if v['view']=='MATCH_OH']
    assert [v['budget'] for v in shift_pair]==[min(5,raw,shift)]*2
    assert [v['budget'] for v in oh_pair]==[min(5,shift,oh)]*2
    assert all(v['zero_budget']==(v['budget']==0) for v in shift_pair+oh_pair)
