from pathlib import Path
from unittest.mock import patch
import pytest
from experiments import follow_expanded_m3 as f
from experiments.summarize_expanded_m3 import bh,representatives,compare

def test_sigma_pool_uses_q_not_raw_p(tmp_path):
 p=tmp_path/'input';p.write_text('test');dest=tmp_path/'out'
 with patch.object(f,'execute',return_value={'exit_code':1}) as run:
  f.scan('case',p,p,p,dest,'strict_sigma_pool','discovery')
 cmd=list(map(str,run.call_args.args[0]));assert '--qv-thresh' in cmd;assert cmd[cmd.index('--thresh')+1]=='0.05'

def test_raw_p_never_q_threshold(tmp_path):
 p=tmp_path/'input';p.write_text('test')
 with patch.object(f,'execute',return_value={'exit_code':1}) as run:f.scan('case',p,p,p,tmp_path/'out','diagnostic_raw_p','development')
 cmd=list(map(str,run.call_args.args[0]));assert '--qv-thresh' not in cmd;assert cmd[cmd.index('--thresh')+1]=='0.0001'

def test_meme_probability_matrix_and_native_e(tmp_path):
 (tmp_path/'meme.xml').write_text('<MEME><motifs><motif id="motif_1" name="A" width="1" sites="10" p_value="0.001" e_value="0.02"><probabilities><alphabet_matrix><alphabet_array><value letter_id="A">0.7</value><value letter_id="C">0.1</value><value letter_id="G">0.1</value><value letter_id="T">0.1</value></alphabet_array></alphabet_matrix></probabilities></motif></motifs></MEME>')
 rows,mx=f.parse_candidates(tmp_path,'meme','new');assert rows[0]['native_significant'];assert rows[0]['motif_id']=='new__M01';assert mx[0][-1][0]['A']==.7

def test_small_streme_test_not_significant(tmp_path):
 (tmp_path/'streme.xml').write_text('<STREME><model><test_positives count="2"/><test_negatives count="2"/></model><motifs><motif id="1-AAA" width="1" test_pvalue="0.001" test_pos_count="2" test_neg_count="2"><pos A="0.7" C="0.1" G="0.1" T="0.1"/></motif></motifs></STREME>')
 rows,_=f.parse_candidates(tmp_path,'streme','group');assert not rows[0]['native_significant'];assert not rows[0]['internal_test_valid']

def test_cluster_representatives_not_chosen_by_hits():
 master={sid:{'leakage_group':cl} for sid,cl in [('a','g1'),('b','g1'),('c','g2'),('d','g2'),('e','g3')]}
 reps,mixed=representatives({'a':1,'b':1,'c':0,'d':1,'e':0},master);assert reps=={'a':1,'e':0};assert mixed==1

def test_null_comparison_is_paired():
 master={str(i):{'leakage_group':str(i)} for i in range(10)};mapping={sid:1 for sid in master}
 r=compare(mapping,master,set(mapping),True);assert r['test']=='paired_exact_binomial';assert r['discordant_positive_only']==10;assert r['p']==pytest.approx(2/1024)

def test_bh_families_and_splits_are_separate():
 rows=[{'split':'discovery','family':'sigma','p':.02},{'split':'discovery','family':'sigma','p':.08},{'split':'holdout','family':'sigma','p':.01},{'split':'discovery','family':'null','p':'NA'}]
 result=bh(rows);assert result[0]['BH_q']==.04;assert result[2]['BH_q']==.01;assert result[2]['family_size']==1;assert 'BH_q' not in result[3]


def test_streme_validity_uses_test_partition_not_hit_count(tmp_path):
 (tmp_path/'streme.xml').write_text('<STREME><model><test_positives count="100"/><test_negatives count="100"/></model><motifs><motif id="1-AAA" width="1" test_pvalue="0.001" test_pos_count="20" test_neg_count="0"><pos A="0.7" C="0.1" G="0.1" T="0.1"/></motif></motifs></STREME>')
 rows,_=f.parse_candidates(tmp_path,'streme','group');assert rows[0]['native_significant'];assert rows[0]['internal_test_valid']


def test_missing_empty_partition_requires_skip_evidence(tmp_path):
 from experiments.summarize_expanded_m3 import partition
 assert partition(tmp_path/'missing.fasta',True)==[]
 with pytest.raises(FileNotFoundError):partition(tmp_path/'missing.fasta')


def test_no_discordant_pairs_still_has_nonzero_uncertainty():
 master={str(i):{'leakage_group':str(i)} for i in range(10)}
 r=compare({sid:1 for sid in master},master,set(),True)
 assert r['rate_difference']==0
 assert r['rate_difference_CI95_lower']<0<r['rate_difference_CI95_upper']
 assert 'discordant_marginals' in r['CI_method']
