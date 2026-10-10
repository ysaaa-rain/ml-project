import numpy as np
from .common import pool_weights,pool
from ..encoder import token_to_base,windows


def test_pooling_definition_and_equivalence():
    token=np.random.default_rng(1).normal(size=(2,76,8)).astype('float32')
    for width in [6,10,16]:
        assert np.allclose(pool(token,width,'A0'),windows(token_to_base(token),width),atol=1e-6)
        a0,a2=pool_weights(width,'A0'),pool_weights(width,'A2')
        assert np.allclose(a0[5:81-width-4],a2[5:81-width-4],atol=1e-7)
        assert np.max(abs(a0[0]-a2[0]))>0.01
        assert np.allclose(a2.sum(axis=1),1)
        assert np.allclose(pool_weights(width,'A1').sum(axis=1),1)


def test_full_token_and_overlap_weights():
    a1=pool_weights(10,'A1')
    assert np.allclose(a1[20,20:25],.2)
    assert a1[20,:20].sum()==0 and a1[20,25:].sum()==0
    a2=pool_weights(10,'A2')[20]
    assert np.allclose(a2[15:30],np.r_[np.arange(1,6),np.repeat(6,5),np.arange(5,0,-1)]/60)


def test_fast_scan_and_projection():
    from .benchmark import fast_scan,projected_start
    from ..discovery import make_pwm,max_scan
    _,pwm=make_pwm(['TTGACA']*20)
    seqs=['C'*20+'TTGACA'+'G'*55,'C'*20+'TGTCAA'+'G'*55]
    assert fast_scan(seqs,pwm,np.ones(4)/4)==max_scan(seqs,pwm,np.ones(4)/4)
    assert projected_start(27,'+',10,12,-2,'+')==(25,'+')
    assert projected_start(25,'-',10,12,-2,'-')==(25,'+')
