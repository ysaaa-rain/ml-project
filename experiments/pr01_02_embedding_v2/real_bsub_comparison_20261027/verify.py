"""Verify source hashes, candidate ranking and group-level hit summaries."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .run import ROOT,OUT,HERE,BRANCH,CNN_ART,STREME,read,read_meme,sha,load_split,background,rank,sequence_windows,scan_windows,groups,rates

def main():
    m=read(OUT/'manifest.json');frozen=read(OUT/'input_lock.json')
    assert m['status']=='COMPLETED' and not m['holdout_used']
    assert frozen['source_sha256']==m['source_sha256']
    assert frozen['ranking_sha256']==sha(OUT/'discovery_rankings.tsv')
    assert sha(HERE/'report.py')==m['report_source_sha256']
    for p,h in m['source_sha256'].items():assert sha(ROOT/p)==h,p
    for p,h in m['outputs_sha256'].items():assert sha(OUT/p)==h,p
    ranking=pd.read_csv(OUT/'discovery_rankings.tsv',sep='\t',dtype={'seed':str})
    family=pd.read_csv(OUT/'family.tsv',sep='\t',dtype={'seed':str})
    counts=pd.read_csv(OUT/'cnn_discovery_positions.tsv',sep='\t',dtype={'seed':str})
    pos,neg,_,mixed=load_split('discovery');dp,dn,shuffle,devmixed=load_split('development');bg=background(neg.sequence.tolist())
    pools=[('STREME','existing',STREME)]
    for seed in [20261005,20261006]:
        pools.append(('CNN_SHIFT',str(seed),CNN_ART/f'seed_{seed}'/'SHIFT'/'motifs.meme'))
        for method in ['onehot','prokbert']:
            pools.append((method.upper(),str(seed),BRANCH/'runs'/f'tjupan_bacillus_subtilis__{method}__L6W10__s{seed}'/'motifs.meme'))
    assert len(pools)==7 and len(family)==14
    for method,seed,path in pools:
        motifs=read_meme(path);computed=rank(motifs,pos,neg,bg,.9,mixed)
        table=ranking[(ranking.method==method)&(ranking.seed==seed)].sort_values('rank')
        assert table.candidate.tolist()==[r['candidate'] for r in computed]
        assert np.allclose(table.threshold,[r['threshold'] for r in computed])
        assert np.allclose(table.discovery_score,[r['discovery_score'] for r in computed])
        for view,chosen in [('TOP5',computed[:5]),('ALL',computed)]:
            flags=[np.zeros(len(dp),bool),np.zeros(len(dn),bool),np.zeros(len(dp),bool)]
            for r in chosen:
                pwm=motifs[r['candidate']]
                for i,seq in enumerate([dp.sequence.tolist(),dn.sequence.tolist(),shuffle]):
                    flags[i]|=scan_windows(sequence_windows(seq,len(pwm)),pwm,bg)[0]>r['threshold']
            pr,nr,pn,nn=rates(dp,dn,flags[0],flags[1],devmixed);pg=groups(dp,flags[0]);sg=groups(dp,flags[2])
            row=family[(family.method==method)&(family.seed==seed)&(family.view==view)].iloc[0]
            assert row.candidates==len(chosen) and row.positive_groups==pn and row.natural_groups==nn
            assert np.allclose([row.positive_rate,row.natural_rate,row.shuffle_rate],[pr,nr,sg.mean()])
            assert row.paired_positive_only==int((pg&~sg).sum()) and row.paired_shuffle_only==int((sg&~pg).sum())
        if method=='CNN_SHIFT':
            instances=pd.read_csv(CNN_ART/f'seed_{seed}'/'SHIFT'/'instances.tsv',sep='\t')
            for item in computed[:3]:
                expected=instances[instances.candidate==item['candidate']];actual=counts[(counts.seed==seed)&(counts.candidate==item['candidate'])]
                assert int(actual.instance_n.sum())==len(expected)
                for r in actual.itertuples():
                    n=((expected.start0+(item['width']-1)/2-60==r.center_tss)&(expected.strand==r.strand)).sum()
                    assert n==r.instance_n
    summary={'status':'PASSED','source_files':len(m['source_sha256']),'outputs':len(m['outputs_sha256']),
       'candidate_pools':len(pools),'family_views_recomputed':len(family),'cnn_instance_position_groups_checked':6,'holdout_used':False}
    (OUT/'verification.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
