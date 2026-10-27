import numpy as np
import pandas as pd
from .pipeline import groups, rates, load_split

def test_group_presence_and_mixed_groups():
    pos=pd.DataFrame({'leakage_group':['a','a','mixed']})
    neg=pd.DataFrame({'leakage_group':['b','mixed']})
    assert groups(pos,np.array([False,True,False])).to_dict()=={'a':True,'mixed':False}
    assert rates(pos,neg,np.array([False,True,False]),np.array([False,True]),{'mixed'})==(1.,0.,1,1)

def test_frozen_real_data_and_paired_shuffles():
    for split in ['discovery','development']:
        pos,neg,shuffle,overlap=load_split(split)
        assert len(shuffle)==len(pos) and len(neg)>0
        for seq,shuf in zip(pos.sequence,shuffle):
            assert len(seq)==len(shuf)==81
            assert sorted(seq)==sorted(shuf)
        assert not overlap
