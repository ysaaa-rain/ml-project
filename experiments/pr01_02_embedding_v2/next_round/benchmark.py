"""Controlled difficulty matrix, with ground truth strictly evaluation-only."""
import json
import itertools
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from sklearn.cluster import MiniBatchKMeans
from sklearn.decomposition import PCA
from threadpoolctl import threadpool_limits
from ..audit import BRANCH,sha,save
from ..encoder import load_model
from ..discovery import make_pwm,background,similarity,write_meme
from .common import NEXT,start,finish,pool

SCENARIOS=[
 {'id':'S1_strong_fixed','gc':.5,'mutation':0,'position':'fixed','motifs':['TTGACATATAAT']},
 {'id':'S2_mut15_fixed','gc':.5,'mutation':.15,'position':'fixed','motifs':['TTGACATATAAT']},
 {'id':'S2_mut30_fixed','gc':.5,'mutation':.30,'position':'fixed','motifs':['TTGACATATAAT']},
 {'id':'S3_strong_random','gc':.5,'mutation':0,'position':'random','motifs':['TTGACATATAAT']},
 {'id':'S4_gc25_random','gc':.25,'mutation':0,'position':'random','motifs':['TTGACATATAAT']},
 {'id':'S4_gc75_random','gc':.75,'mutation':0,'position':'random','motifs':['TTGACATATAAT']},
 {'id':'S5_pair_fixed','gc':.5,'mutation':0,'position':'random','motifs':['TTGACA','TATAAT'],'gap':[17,17]},
 {'id':'S6_pair_variable','gc':.5,'mutation':0,'position':'random','motifs':['TTGACA','TATAAT'],'gap':[10,24]},
 {'id':'S7_markov_random','gc':.5,'mutation':0,'position':'random','motifs':['TTGACATATAAT'],'persistence':.7}
]


def fast_scan(sequences,pwm,bg):
    array=np.array([['ACGT'.index(c) for c in s] for s in sequences],dtype=np.int8)
    window=np.lib.stride_tricks.sliding_window_view(array,len(pwm),axis=1)
    forward=np.log2(pwm/bg)[np.arange(len(pwm)),window].sum(axis=-1)
    reverse=np.log2(pwm[::-1,::-1]/bg)[np.arange(len(pwm)),window].sum(axis=-1)
    both=np.concatenate([forward,reverse],axis=1)
    idx=both.argmax(axis=1); n=forward.shape[1]
    return [(float(both[r,i]),int(i%n),'+' if i<n else '-') for r,i in enumerate(idx)]


def generate(n,scenario,rng,plant):
    gc=scenario['gc'];prob=np.array([(1-gc)/2,gc/2,gc/2,(1-gc)/2])
    sequences=[]; truth=[]
    for i in range(n):
        bases=rng.choice(list('ACGT'),81,p=prob)
        persistence=scenario.get('persistence',0)
        for j in range(1,81):
            if rng.random()<persistence:bases[j]=bases[j-1]
        if plant:
            motifs=scenario['motifs']
            if len(motifs)==1:
                positions=[25 if scenario['position']=='fixed' else int(rng.integers(0,82-len(motifs[0])))]
            else:
                gap=int(rng.integers(scenario['gap'][0],scenario['gap'][1]+1))
                left=int(rng.integers(0,82-sum(map(len,motifs))-gap))
                positions=[left,left+len(motifs[0])+gap]
            for m,(motif,pos) in enumerate(zip(motifs,positions)):
                planted=list(motif)
                for j,c in enumerate(planted):
                    if rng.random()<scenario['mutation']:planted[j]=rng.choice([x for x in 'ACGT' if x!=c])
                bases[pos:pos+len(motif)]=planted
                truth.append({'sequence_index':i,'family':m,'start0':pos,'end0':pos+len(motif),'strand':'+','planted_sequence':''.join(planted)})
        sequences.append(''.join(bases))
    return sequences,truth


def encode(sequences,model,tokenizer,batch_size=32):
    import torch
    rows=[]
    with torch.inference_mode():
        for start in range(0,len(sequences),batch_size):
            seqs=sequences[start:start+batch_size]
            inputs=tokenizer(seqs,return_tensors='pt',padding=True)
            for ids,seq in zip(inputs['input_ids'].tolist(),seqs):
                assert ids[1:-1]==[tokenizer.convert_tokens_to_ids(seq[i:i+6]) for i in range(76)]
            rows.append(model(**inputs,output_hidden_states=True).hidden_states[6][:,1:-1].numpy())
    return np.concatenate(rows)


def candidates(features,sequences,width,seed):
    flat=features.reshape(-1,features.shape[-1]);nw=82-width
    with threadpool_limits(limits=4):
        pca=PCA(n_components=min(24,flat.shape[1]),svd_solver='randomized',random_state=seed)
        z=pca.fit_transform(flat)
        km=MiniBatchKMeans(n_clusters=24,n_init=5,batch_size=2048,random_state=seed)
        lab=km.fit_predict(z)
        dist=np.linalg.norm(z-km.cluster_centers_[lab],axis=1)
    motifs,instance_rows={},[]
    for cluster in range(24):
        ids=np.flatnonzero(lab==cluster);selected={}
        for idx in ids[np.argsort(dist[ids],kind='stable')]:selected.setdefault(int(idx//nw),int(idx%nw))
        if len(selected)<20:continue
        name=f'C{cluster:02d}'
        _,pwm=make_pwm([sequences[i][pos:pos+width] for i,pos in selected.items()])
        motifs[name]=pwm
        instance_rows += [{'motif':name,'sequence_index':i,'start0':pos,'end0':pos+width,
                           'fragment':sequences[i][pos:pos+width]} for i,pos in selected.items()]
    return motifs,instance_rows


def projected_start(start,scan_strand,width,target_length,offset,alignment_strand):
    if scan_strand=='+':position=start+offset
    else:position=start+width-offset-target_length
    direction='+' if scan_strand==alignment_strand else '-'
    return position,direction


def evaluate(motifs,scenario,data,width,bg,strict=False):
    names=list(motifs); references=[]
    for sequence in scenario['motifs']:
        _,p=make_pwm([sequence]*256)
        mutation=scenario['mutation']
        if mutation:
            p=np.full((len(sequence),4),mutation/3)
            for i,c in enumerate(sequence):p[i,'ACGT'.index(c)]=1-mutation
        references.append(p)
    comparisons=[[similarity(motifs[n],r,minimum_overlap=min(len(r),width) if strict else min(6,len(r),width)) for n in names] for r in references]
    scores=np.array([[x[0] for x in row] for row in comparisons])
    families,indices=linear_sum_assignment(-scores)
    mapping={int(f):(names[int(i)],comparisons[int(f)][int(i)]) for f,i in zip(families,indices)}
    metrics,predictions,matches=[],[],[]
    for family,(name,(score,offset,orientation)) in mapping.items():
        pwm=motifs[name]
        threshold=float(np.quantile([x[0] for x in fast_scan(data['calibration_control'],pwm,bg)],.95,method='higher'))
        eligible=score>=.8
        matches.append({'family':family,'candidate':name,'cosine':score,'offset':offset,'alignment_strand':orientation,
                        'eligible':eligible,'threshold':threshold,'target_length':len(references[family]),'candidate_width':width})
        for split in ['development','test']:
            truth={r['sequence_index']:r for r in data[split+'_truth'] if r['family']==family}
            pos=fast_scan(data[split],pwm,bg);neg=fast_scan(data[split+'_control'],pwm,bg)
            n_pred_pos,n_pred_control=0,0;counts={t:0 for t in [0,1,2,3]};overlaps=[]
            for i,(value,start,strand) in enumerate(pos):
                predicted,canonical=projected_start(start,strand,width,len(references[family]),offset,orientation)
                yes=bool(value>threshold and eligible);n_pred_pos+=yes
                actual=truth[i];delta=predicted-actual['start0']
                intersection=max(0,min(predicted+len(references[family]),actual['end0'])-max(predicted,actual['start0']))
                iou=intersection/(2*len(references[family])-intersection)
                for t in counts:counts[t]+=int(yes and canonical==actual['strand'] and abs(delta)<=t)
                if yes:overlaps.append(iou)
                predictions.append({'split':split,'label':1,'sequence_index':i,'family':family,'candidate':name,
                       'present':yes,'raw_window_start0':start,'projected_truth_start0':predicted,'canonical_strand':canonical,
                       'truth_start0':actual['start0'],'delta':delta,'score':value,'interval_iou':iou})
            for i,(value,start,strand) in enumerate(neg):
                yes=bool(value>threshold and eligible);n_pred_control+=yes
                predictions.append({'split':split,'label':0,'sequence_index':i,'family':family,'candidate':name,'present':yes,'score':value})
            for t,tp in counts.items():
                precision=tp/(n_pred_pos+n_pred_control) if n_pred_pos+n_pred_control else 0
                recall=tp/len(truth)
                metrics.append({'split':split,'family':family,'candidate':name,'tolerance_bp':t,'precision':precision,
                    'recall':recall,'f1':2*precision*recall/(precision+recall) if precision+recall else 0,
                    'true_positives':tp,'prediction_count_positive':n_pred_pos,'prediction_count_control':n_pred_control,
                    'truth_count':len(truth),'pwm_similarity':score,'reference_eligible':eligible,
                    'mean_interval_iou_of_predictions':float(np.mean(overlaps)) if overlaps else 0})
    grammar=[]
    if len(references)==2:
        preds=pd.DataFrame(predictions)
        for split in ['development','test']:
            sub=preds[(preds.split==split)&(preds.label==1)]
            for i,g in sub.groupby('sequence_index'):
                if len(g)!=2:continue
                g=g.set_index('family');a,b=g.loc[0],g.loc[1]
                truthgap=b.truth_start0-a.truth_start0-len(references[0])
                recovered=bool(a.present and b.present and a.canonical_strand=='+' and b.canonical_strand=='+')
                gap=b.projected_truth_start0-a.projected_truth_start0-len(references[0])
                grammar.append({'split':split,'sequence_index':i,'both_detected':recovered,'true_gap':truthgap,
                  'predicted_gap':gap if recovered else np.nan,'correct_order':bool(recovered and b.projected_truth_start0>a.projected_truth_start0),
                  'gap_error':gap-truthgap if recovered else np.nan})
    return metrics,predictions,matches,grammar


def main():
    config=json.loads((NEXT/'protocol.json').read_text())['synthetic']
    target,manifest=start('synthetic_difficulty',config,[NEXT/'protocol.json',BRANCH/'config.json'])
    cfg=json.loads((BRANCH/'config.json').read_text(encoding='utf-8'))
    model,tokenizer=load_model(cfg)
    metric_rows=[];output_files=[]
    for scenario in SCENARIOS:
        for seed in config['seeds']:
            path=NEXT/'synthetic_benchmark'/f'{scenario["id"]}__s{seed}'
            path.mkdir(exist_ok=False)
            rng=np.random.default_rng(seed)
            data={}
            for split,n in [('discovery',256),('development',128),('test',128)]:
                data[split],data[split+'_truth']=generate(n,scenario,rng,True)
                data[split+'_control'],_=generate(n,scenario,rng,False)
            data['calibration_control'],_=generate(256,scenario,rng,False)
            save(path/'data.json',{'scenario':scenario,'seed':seed,**data})
            output_files.append(path/'data.json')
            tokens=encode(data['discovery'],model,tokenizer)
            np.save(path/'token_L6.npy',tokens)
            output_files.append(path/'token_L6.npy')
            widths=[6] if len(scenario['motifs'])==2 else [10]
            if scenario['id']=='S3_strong_random':widths=[6,10,16]
            bg=background(data['calibration_control'])
            for width in widths:
                for method in config['methods']:
                    run=path/f'{method}__w{width}';run.mkdir()
                    if method=='OH':
                        features=np.array([np.eye(4,dtype=np.float32)[['ACGT'.index(c) for c in s[i:i+width]]].ravel()
                                           for s in data['discovery'] for i in range(82-width)]).reshape(256,82-width,-1)
                    else:features=pool(tokens,width,method)
                    motifs,instances=candidates(features,data['discovery'],width,seed)
                    write_meme(run/'motifs.meme',motifs,bg)
                    pd.DataFrame(instances).to_csv(run/'instances.tsv',sep='\t',index=False)
                    metrics,preds,matches,grammar=evaluate(motifs,scenario,data,width,bg)
                    for row in metrics:row.update(scenario=scenario['id'],seed=seed,method=method,width=width,candidate_count=len(motifs))
                    metric_rows+=metrics
                    pd.DataFrame(metrics).to_csv(run/'metrics.tsv',sep='\t',index=False)
                    pd.DataFrame(preds).to_csv(run/'predictions.tsv',sep='\t',index=False)
                    pd.DataFrame(matches).to_csv(run/'reference_matches.tsv',sep='\t',index=False)
                    if grammar:pd.DataFrame(grammar).to_csv(run/'grammar.tsv',sep='\t',index=False)
                    output_files += [p for p in run.iterdir() if p.is_file()]
            print(f'Completed synthetic {scenario["id"]} seed={seed}',flush=True)
    pd.DataFrame(metric_rows).to_csv(NEXT/'synthetic_benchmark/all_metrics.tsv',sep='\t',index=False)
    output_files.append(NEXT/'synthetic_benchmark/all_metrics.tsv')
    finish(target,manifest,output_files,{'scenarios':SCENARIOS,'conditions':27,'method_fits':132,
        'truth_usage':'evaluation-only PWM mapping and positional metrics; never PCA, clustering or instance selection',
        'synthetic_test_limit':'within this version, no selection from test; future changes require new independent synthetic seeds'})


if __name__=='__main__':main()
