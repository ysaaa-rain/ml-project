"""Independent-data and initialization aggregation checks."""
import pandas as pd
import pytest

from .pipeline import PROTOCOL,calculate_init_effect,evaluate_gate,read,views


@pytest.mark.parametrize('raw,shift,oh',[(0,7,6),(2,7,3),(7,7,7)])
def test_candidate_budgets_keep_zero_and_underfill(raw,shift,oh):
    rankings={m:[{'candidate':str(i)} for i in range(n)] for m,n in
              [('TIED_RAW',raw),('TIED_SHIFT',shift),('OH48_BG',oh)]}
    all_views=views(rankings);assert len(all_views)==10
    shift_pair=[r for r in all_views if r['view']=='MATCH_SHIFT']
    oh_pair=[r for r in all_views if r['view']=='MATCH_OH']
    assert [r['budget'] for r in shift_pair]==[min(5,raw,shift)]*2
    assert [r['budget'] for r in oh_pair]==[min(5,shift,oh)]*2
    assert all(r['zero_budget']==(r['budget']==0) for r in shift_pair+oh_pair)


def test_init_effect_uses_paired_data_and_does_not_count_models_as_data():
    rows=[]
    for data_seed in range(5):
        for init in range(6):
            for name,value in [('TIED_RAW',.3),('TIED_SHIFT',.3+init/100)]:
                rows.append({'scenario':'mut30_random','data_seed':data_seed,'algorithm_seed':init,
                             'view':'MATCH_SHIFT','representation':name,'f1':value})
    effects=calculate_init_effect(pd.DataFrame(rows)).set_index('algorithm_seed')
    assert len(effects)==6 and set(effects.data_seeds)=={5}
    assert effects.loc[3,'mean_delta']==pytest.approx(.03)
    assert effects.loc[0,'data_wins']==0 and effects.loc[3,'data_wins']==5


def test_gate_requires_three_new_initializations_and_checks_blank_fpr():
    paired=pd.DataFrame({'delta':[.1,.1,.1,.1,0.], 'zero_budget_count':[0]*5,
                         'family_fpr_shift':[.1]*5,'family_fpr_baseline':[.2]*5})
    init=pd.DataFrame({'scenario':['mut30_random']*6,'algorithm_seed':[12,13,14,15,16,17],
                       'mean_delta':[.4,0.,.2,.1,.1,-.1]})
    strong=pd.DataFrame({'delta':[0.]*5})
    gate=evaluate_gate(paired,paired,strong,init,[14,15,16,17]);assert all(gate.values())
    init.loc[init.algorithm_seed==16,'mean_delta']=0.
    assert not evaluate_gate(paired,paired,strong,init,[14,15,16,17])['at_least3_of4_new_initializations_positive']
    paired['family_fpr_shift']=.3
    assert not evaluate_gate(paired,paired,strong,init,[14,15,16,17])['family_fpr_not_above_raw']


def test_data_seeds_are_fresh_and_initialization_subgroups_disjoint():
    cfg=read(PROTOCOL);previous=read(PROTOCOL.parents[1]/'cnn_boundary_20261027/protocol.json')
    assert set(cfg['data_seeds']).isdisjoint(previous['data_seeds'])
    assert set(cfg['reused_initialization_seeds'])==set(previous['algorithm_seeds'])
    assert set(cfg['new_initialization_seeds']).isdisjoint(previous['algorithm_seeds'])
    assert set(cfg['algorithm_seeds'])==set(cfg['new_initialization_seeds']+cfg['reused_initialization_seeds'])
