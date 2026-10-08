from experiments.analyze_frozen_m3_grammar import bh,gap,choose_sites,fisher
from preprocessing.build_reference_library_v1 import RedParser
import pytest

def test_bh_preserves_original_order_and_monotonicity():
 assert bh([.04,.001,.03,.9])==pytest.approx([.0533333333,.004,.0533333333,.9])

def test_gap_adjacent_overlap_and_nested_sites():
 def s(a,b):return {'start':str(a),'stop':str(b),'strand':'+'}
 assert gap(s(4,9),s(10,15))==0
 assert gap(s(4,10),s(8,13))==-3
 assert gap(s(4,20),s(9,14))==-12

def test_primary_site_priority_q_score_start():
 def r(q,score,start):return {'sequence_name':'s','motif_id':'m','q-value':str(q),'score':str(score),'start':str(start),'stop':str(start+5),'strand':'+'}
 rs=[r(.01,8,20),r(.01,10,15),r(.01,10,5),r(.1,50,1)]
 assert choose_sites(rs)[('s','m')]['start']=='5'

def test_no_hits_is_not_positive_effect():
 r=fisher(0,100,0,80)
 assert r['raw_p']==1 and r['odds_ratio']=='NA'

def test_red_short_span_preserves_offset_and_publication_link():
 p=RedParser();p.feed('<table><tr><td>SigA</td><td>Promoter</td><td>-12:+1</td><td>1..13</td><td><tt>AA<font color=red>TA TAAT</font>GG</tt></td><td><a href="?list_uids=123">study</a>: PE</td></tr></table>')
 assert p.rows[0][4]['red']==list(range(2,8))
 assert p.rows[0][5]['links']==['?list_uids=123']
