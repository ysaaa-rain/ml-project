import numpy as np
import pytest
from .encoder import spans, token_to_base, windows
from .discovery import make_pwm, max_scan, background, similarity, reverse_pwm, bh
from .audit import load


def test_coordinates_and_edges():
    hidden = np.arange(76,dtype=np.float32)[:,None]
    base = token_to_base(hidden)
    assert spans(81)[-1] == (75,81)
    assert base[0,0] == 0 and base[-1,0] == 75
    assert base[60,0] == np.mean(np.arange(55,61))
    assert windows(base,10).shape == (72,1)
    assert 61-61 == 60-60 == 0
    with pytest.raises(ValueError):
        token_to_base(np.zeros((78,4)))


def test_holdout_rejected():
    with pytest.raises(ValueError):
        load('tjupan_bacillus_subtilis','holdout')


def test_synthetic_recovery_both_strands_and_pwm():
    _,pwm = make_pwm(['TTGACA']*50)
    seq = 'C'*25+'TTGACA'+'G'*50
    scan = max_scan([seq],pwm,np.full(4,.25))[0]
    assert scan[1:] == (25,'+')
    rc = 'C'*25+'TGTCAA'+'G'*50
    assert max_scan([rc],pwm,np.full(4,.25))[0][1:] == (25,'-')
    assert np.allclose(pwm.sum(axis=1),1)
    assert similarity(pwm,reverse_pwm(pwm),6)[0] == pytest.approx(1)
    assert np.allclose(bh([.01,.04,.03]),[.03,.04,.04])


def test_shuffle_preserves_dinucleotides():
    from collections import Counter
    from .audit import fasta, DATA
    directory = DATA/'tjupan_bacillus_subtilis/main'
    original = dict(fasta(directory/'discovery.positive.fasta'))
    for key,sequence in fasta(directory/'discovery.dinucleotide_null.fasta'):
        parent = original[key.removesuffix('__dinucleotide_null')]
        assert Counter(a+b for a,b in zip(sequence,sequence[1:])) == Counter(a+b for a,b in zip(parent,parent[1:]))
