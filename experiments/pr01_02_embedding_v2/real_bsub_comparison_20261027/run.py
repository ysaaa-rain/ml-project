"""Rank and scan all existing real-data candidates with one common rule."""
import hashlib
import json
import math
from pathlib import Path
import numpy as np
import pandas as pd

from ..cnn_real_pilot_20261027.pipeline import (ROOT,DATA,ART as CNN_ART,SUMMARY as CNN_SUMMARY,
    background,groups,load_split,rank,rates,read,read_meme,scan_windows,sequence_windows,sha,validate_lock)
from ..discovery import similarity

HERE=Path(__file__).resolve().parent;OUT=HERE/'summary'
STREME=ROOT/'results/motif/tju_streme_v2_20261008/tjupan_bacillus_subtilis/natural/streme.txt'
STREME_MANIFEST=ROOT/'results/motif/tju_streme_v2_20261008/public_evidence/run_manifest.json'
BRANCH=HERE.parent;PROTOCOL=HERE/'protocol.json'

def main():
    cfg=read(PROTOCOL);lock=validate_lock();source={}
    def keep(p):source[p.relative_to(ROOT).as_posix()]=sha(p);return p
    for p in [PROTOCOL,Path(__file__),DATA/'sequence_master.tsv',DATA/'discovery.dinucleotide_null.fasta',
              DATA/'development.dinucleotide_null.fasta',STREME,STREME_MANIFEST,CNN_SUMMARY/'fit_lock.json']:
        keep(p)
    original=read(STREME_MANIFEST);assert sha(STREME)==original['output_sha256']['tjupan_bacillus_subtilis/natural/streme.txt']
    pwms=read_meme(STREME);assert len(pwms)==10
    pools=[('STREME','existing',pwms)]
    for seed in cfg['algorithm_seeds']:
        path=CNN_ART/f'seed_{seed}'/'SHIFT'/'motifs.meme';keep(path);keep(path.parent/'instances.tsv');pools.append(('CNN_SHIFT',str(seed),read_meme(path)))
        for method in ['onehot','prokbert']:
            folder=BRANCH/'runs'/f'tjupan_bacillus_subtilis__{method}__L6W10__s{seed}'
            manifest=read(keep(folder/'manifest.json'));path=keep(folder/'motifs.meme')
            assert sha(path)==manifest['output_sha256']['motifs.meme'] and not manifest['holdout_used']
            key=str((DATA/'sequence_master.tsv').relative_to(ROOT)).replace('/','\\')
            assert sha(DATA/'sequence_master.tsv')==manifest['input_hashes'][key]
            pools.append((method.upper(),str(seed),read_meme(path)))
    pos,neg,_,mixed=load_split('discovery');devpos,devneg,shuffle,devmixed=load_split('development')
    bg=background(neg.sequence.tolist());rankings=[]
    for method,seed,motifs in pools:
        ranking=rank(motifs,pos,neg,bg,.9,mixed)
        assert len(ranking)>=5
        rankings.append((method,seed,motifs,ranking))
    OUT.mkdir(exist_ok=False)
    rank_rows=[{'method':method,'seed':seed,'rank':i,**r} for method,seed,_,ranking in rankings for i,r in enumerate(ranking,1)]
    pd.DataFrame(rank_rows).to_csv(OUT/'discovery_rankings.tsv',sep='\t',index=False,lineterminator='\n')
    lock_path=OUT/'input_lock.json';lock_path.write_text(json.dumps({'source_sha256':source,'ranking_sha256':sha(OUT/'discovery_rankings.tsv'),
       'data_note':'Discovery and Development were previously observed; retrospective exploration'},indent=2)+'\n',encoding='utf8')
    rows=[];individual=[];similar=[];locations=[]
    for method,seed,motifs,ranking in rankings:
        hits={};wins={}
        for item in ranking:
            name=item['candidate'];pwm=motifs[name];width=len(pwm)
            if width not in wins:wins[width]=[sequence_windows(x,width) for x in [devpos.sequence.tolist(),devneg.sequence.tolist(),shuffle]]
            hits[name]=[scan_windows(x,pwm,bg)[0]>item['threshold'] for x in wins[width]]
            pr,nr,pn,nn=rates(devpos,devneg,hits[name][0],hits[name][1],devmixed)
            pg=groups(devpos,hits[name][0]);sg=groups(devpos,hits[name][2]);individual.append({'method':method,'seed':seed,'candidate':name,
              'positive_groups':pn,'natural_groups':nn,'positive_rate':pr,'natural_rate':nr,'shuffle_rate':float(sg.mean()),'paired_groups':len(pg)})
        for view,names in [('TOP5',[x['candidate'] for x in ranking[:5]]),('ALL',[x['candidate'] for x in ranking])]:
            flags=[np.logical_or.reduce([hits[n][i] for n in names]) for i in range(3)]
            pr,nr,pn,nn=rates(devpos,devneg,flags[0],flags[1],devmixed);pg=groups(devpos,flags[0]);sg=groups(devpos,flags[2])
            rows.append({'method':method,'seed':seed,'view':view,'candidates':len(names),'positive_groups':pn,'natural_groups':nn,
              'positive_rate':pr,'natural_rate':nr,'natural_difference':pr-nr,'shuffle_rate':float(sg.mean()),
              'paired_positive_only':int((pg&~sg).sum()),'paired_shuffle_only':int((sg&~pg).sum()),'paired_groups':len(pg)})
        if method=='CNN_SHIFT':
            instance_path=CNN_ART/f'seed_{seed}'/'SHIFT'/'instances.tsv';keep(instance_path)
            inst=pd.read_csv(instance_path,sep='\t')
            for r in ranking[:5]:
                name=r['candidate'];pwm=motifs[name];matches=[]
                for refname,ref in pools[0][2].items():
                    min_overlap=max(5,math.ceil(.6*min(len(pwm),len(ref))))
                    score,offset,strand=similarity(pwm,ref,minimum_overlap=min_overlap)
                    overlap=max(0,min(len(pwm),offset+len(ref))-max(0,offset))
                    matches.append((score,refname,overlap,offset,strand,len(ref)))
                best=max(matches,key=lambda x:(x[0],x[1]));similar.append({'seed':seed,'candidate':name,'streme':best[1],'cosine':best[0],
                    'overlap_bp':best[2],'offset':best[3],'orientation':best[4],'reference_width':best[5]})
            for r in ranking[:3]:
                rows_inst=inst[inst.candidate==r['candidate']]
                for (start,strand),g in rows_inst.groupby(['start0','strand']):
                    locations.append({'seed':seed,'candidate':r['candidate'],'center_tss':float(start+(r['width']-1)/2-60),
                      'strand':strand,'instance_n':len(g)})
    for name,frame in [('family.tsv',pd.DataFrame(rows)),('candidates.tsv',pd.DataFrame(individual)),('cnn_streme_similarity.tsv',pd.DataFrame(similar)),('cnn_discovery_positions.tsv',pd.DataFrame(locations))]:
        frame.to_csv(OUT/name,sep='\t',index=False,lineterminator='\n')
    outputs={x.name:sha(x) for x in OUT.iterdir() if x.is_file()}
    manifest={'status':'COMPLETED','source_sha256':source,'outputs_sha256':outputs,'source_count':len(source),'candidate_pools':len(pools),
      'family_rows':len(rows),'holdout_used':False,'limitations':cfg['limitations']}
    (OUT/'manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False)+'\n',encoding='utf8')
    print(pd.DataFrame(rows).query("view == 'TOP5'").to_string(index=False),flush=True)

if __name__=='__main__':main()
