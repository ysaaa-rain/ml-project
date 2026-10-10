"""CPU baseline. Reuses only verified cache; fails closed on loading errors."""
import argparse
import itertools
import json
from pathlib import Path
import subprocess
import time
import traceback
import joblib
import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.cluster import MiniBatchKMeans
from threadpoolctl import threadpool_limits
from .audit import ROOT, BRANCH, DATA, load, fasta, save, sha
from .encoder import load_model, extract
from .discovery import (background,make_pwm,max_scan,similarity,write_meme,read_meme,
                        group_presence,test_enrichment,bh)


def logo(path, pwm, name):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.textpath import TextPath
    from matplotlib.patches import PathPatch
    from matplotlib.transforms import Affine2D
    colors = dict(A='#228B22',C='#2469c8',G='#e1a500',T='#d13e44')
    heights = pwm*(2+(pwm*np.log2(pwm)).sum(axis=1))[:,None]
    fig,ax = plt.subplots(figsize=(max(5,len(pwm)*.45),3))
    for x,row in enumerate(heights):
        bottom = 0
        for idx in np.argsort(row):
            char,h = 'ACGT'[idx],row[idx]
            glyph = TextPath((0,0),char,size=1,prop={'weight':'bold'})
            bounds = glyph.get_extents()
            trans = Affine2D().translate(-bounds.x0,-bounds.y0).scale(.9/bounds.width,h/bounds.height).translate(x+.05,bottom)
            ax.add_patch(PathPatch(glyph, transform=trans+ax.transData, color=colors[char],lw=0))
            bottom += h
    ax.set(xlim=(0,len(pwm)),ylim=(0,2),ylabel='Information (bits)',xlabel='Aligned position',title=name)
    ax.set_xticks(np.arange(len(pwm))+.5,np.arange(1,len(pwm)+1))
    fig.tight_layout(); fig.savefig(path,dpi=160); plt.close(fig)


def run(source, seed, representation, model=None, tokenizer=None, layer=None, width=None):
    cfg = json.loads((BRANCH/'config.json').read_text(encoding='utf-8'))
    cfg['seed'] = seed
    if layer is not None:
        cfg['baseline_layer'] = layer
    if width is not None:
        cfg['baseline_window'] = width
    if representation == 'multiscale':
        cfg['baseline_window'] = max(cfg['windows'])
    tag = f'{source}__{representation}__L{cfg["baseline_layer"]}W{cfg["baseline_window"]}__s{seed}'
    out = BRANCH/'runs'/tag
    out.mkdir(parents=True,exist_ok=False)
    started = time.time()
    manifest = {'status':'RUNNING','config':cfg,'source':source,'representation':representation,
                'holdout_used':False,'started_unix':started,
                'device':'cpu','environment_sha256':sha(BRANCH/'audit/environment.json'),
                'dependencies':json.loads((BRANCH/'audit/environment.json').read_text(encoding='utf-8'))['packages'],
                'git_head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                'code_sha256':{p.name:sha(p) for p in BRANCH.glob('*.py')},
                'data_manifest_sha256':sha(DATA/'manifest.json')}
    save(out/'manifest.json',manifest)
    try:
        discovery = load(source,'discovery',1)
        natural_discovery = load(source,'discovery',0)
        dev = load(source,'development',1)
        natural_dev = load(source,'development',0)
        ids = {r.sequence_id:r.leakage_group for r in dev.itertuples()}
        null_records = fasta(DATA/source/'main/development.dinucleotide_null.fasta')
        null_dev = pd.DataFrame([{'sequence_id':i,'sequence':s,'leakage_group':ids[i.removesuffix('__dinucleotide_null')]}
                                 for i,s in null_records])
        width = cfg['baseline_window']
        if representation in ('prokbert','multiscale','tokenwindow'):
            cache = BRANCH/'cache'/source
            if (cache/'cache_manifest.json').exists():
                cached = json.loads((cache/'cache_manifest.json').read_text(encoding='utf-8'))
                assert cached['ordered_ids'] == discovery.sequence_id.tolist()
                for key in ['revision','code_revision','layers','windows']:
                    assert cached['config'][key] == cfg[key]
                assert cached['files']['index.tsv'] == sha(cache/'index.tsv')
                path = cache/f'window_L{cfg["baseline_layer"]}_W{width}.npy'
                assert cached['files'][path.name] == sha(path)
            else:
                extract(discovery,cache,cfg,model,tokenizer)
                path = cache/f'window_L{cfg["baseline_layer"]}_W{width}.npy'
            if representation == 'tokenwindow':
                tokenpath = cache/f'token_L{cfg["baseline_layer"]}.npy'
                assert cached['files'][tokenpath.name] == sha(tokenpath)
                token = np.load(tokenpath,mmap_mode='r')
                features = np.stack([token[:,i:i+width-5].mean(axis=1) for i in range(82-width)],axis=1).reshape(-1,384)
            elif representation == 'multiscale':
                aligned = []
                for w in cfg['windows']:
                    p = cache/f'window_L{cfg["baseline_layer"]}_W{w}.npy'
                    assert cached['files'][p.name] == sha(p)
                    shift = (width-w)//2
                    a = np.load(p,mmap_mode='r')[:,shift:shift+82-width]
                    aligned.append(a)
                features = np.concatenate(aligned,axis=-1).reshape(len(discovery)*(82-width),-1)
            else:
                features = np.load(path,mmap_mode='r').reshape(-1,model.config.hidden_size if model else 384)
        else:
            alphabet = 'ACGT'
            features = np.stack([np.eye(4,dtype=np.float32)[[alphabet.index(c) for c in seq[i:i+width]]].ravel()
                                 for seq in discovery.sequence for i in range(82-width)])
        with threadpool_limits(limits=cfg['threads']):
            pca = PCA(n_components=min(cfg['pca_components'],features.shape[1]),svd_solver='randomized',random_state=seed)
            reduced = pca.fit_transform(features)
            km = MiniBatchKMeans(n_clusters=cfg['clusters'],batch_size=2048,n_init=5,random_state=seed)
            labels = km.fit_predict(reduced)
            distances = np.linalg.norm(reduced-km.cluster_centers_[labels],axis=1)
        joblib.dump({'pca':pca,'kmeans':km},out/'fitted.joblib')
        bg = background(natural_discovery.sequence)
        motif_rows, instance_rows, motifs, merges = [], [], {}, []
        n_windows = 82-width
        for cluster in range(cfg['clusters']):
            members = np.flatnonzero(labels == cluster)
            # One closest window per promoter and per cluster; overlapping windows never inflate PWM counts.
            selected = {}
            for idx in members[np.argsort(distances[members],kind='stable')]:
                selected.setdefault(int(idx//n_windows),int(idx))
            if len(selected) < cfg['min_promoters']:
                continue
            name = f'EMB_{representation}_C{cluster:02d}'
            instances = []
            for row,idx in selected.items():
                start = idx % n_windows
                r = discovery.iloc[row]
                seq = r.sequence[start:start+width]
                instances.append(seq)
                instance_rows.append({'motif':name,'sequence_id':r.sequence_id,'leakage_group':r.leakage_group,
                      'split':'discovery','start0':start,'end0_exclusive':start+width,'start1':start+1,
                      'stop1':start+width,'center_tss':start+(width-1)/2-60,'strand':'+',
                      'sequence':seq,'distance':float(distances[idx])})
            counts,pwm = make_pwm(instances,cfg['pseudocount'])
            motifs[name] = pwm
            np.savetxt(out/f'{name}.pfm.tsv',counts,delimiter='\t')
            motif_rows.append({'motif':name,'promoters':len(instances),'windows':len(members),
                               'consensus': ''.join('ACGT'[i] for i in pwm.argmax(axis=1)),
                               'information_bits':float((2+(pwm*np.log2(pwm)).sum(axis=1)).sum())})
        # Conservative full-length similarity avoids collapsing arbitrary short overlaps.
        retained = {}
        for name,pwm in motifs.items():
            match = next((other for other,m in retained.items() if similarity(pwm,m,minimum_overlap=width)[0]>=.95),None)
            if match:
                merges.append({'motif':name,'representative':match})
            else:
                retained[name] = pwm
        write_meme(out/'motifs.meme',retained,bg)
        pd.DataFrame(instance_rows).to_csv(out/'instances.tsv',sep='\t',index=False)
        pd.DataFrame(motif_rows).to_csv(out/'candidates.tsv',sep='\t',index=False)
        save(out/'merge_map.json',merges)
        for name,pwm in retained.items():
            logo(out/f'{name}.png',pwm,name)
        evaluation, hits, thresholds, comparisons = [], [], {}, []
        traditional_path = ROOT/'results/motif/tju_streme_v2_20261008'/source/'selected_main.meme'
        traditional = read_meme(traditional_path)
        all_motifs = {**retained, **{f'STREME_{k}':v for k,v in traditional.items()}}
        for name,pwm in all_motifs.items():
            calibration = np.array([v[0] for v in max_scan(natural_discovery.sequence,pwm,bg)])
            threshold = float(np.quantile(calibration,1-cfg['development_control_fpr'],method='higher'))
            thresholds[name] = threshold
            presence = {}
            for control, frame in [('positive',dev),('natural',natural_dev),('shuffle',null_dev)]:
                scan = max_scan(frame.sequence,pwm,bg)
                present = [v[0]>threshold for v in scan]
                presence[control] = group_presence(frame,present)
                for row,(score,start,strand),yes in zip(frame.itertuples(),scan,present):
                    hits.append({'motif':name,'control':control,'sequence_id':row.sequence_id,
                      'leakage_group':row.leakage_group,'score':score,'present':yes,'start1':start+1,
                      'stop1':start+len(pwm),'strand':strand,'center_tss':start+(len(pwm)-1)/2-60})
            for control in ['natural','shuffle']:
                shared = set(presence['positive'].index)&set(presence[control].index)
                if control == 'natural':
                    # Mixed-label related groups cannot supply independent cells to Fisher.
                    pos = presence['positive'].drop(list(shared))
                    neg = presence[control].drop(list(shared))
                    stats = test_enrichment(pos,neg)
                    stats['shared_groups_excluded'] = len(shared)
                else:
                    # Shuffle is bound to each promoter group: paired exact McNemar test.
                    from scipy.stats import binomtest
                    pos = presence['positive'].loc[sorted(shared)].values
                    neg = presence[control].loc[sorted(shared)].values
                    b, c = int((pos&~neg).sum()), int((~pos&neg).sum())
                    stats = {'positive_groups':len(pos),'control_groups':len(neg),'positive_hits':int(pos.sum()),
                             'control_hits':int(neg.sum()), 'p':float(binomtest(b,b+c,.5,alternative='greater').pvalue) if b+c else 1.,
                             'discordant_pos_only':b,'discordant_null_only':c}
                evaluation.append({'motif':name,'method':'STREME' if name.startswith('STREME_') else representation,
                         'control':control,'threshold':threshold,'test':'paired_binomial' if control=='shuffle' else 'group_fisher',**stats})
        evaluation = pd.DataFrame(evaluation)
        for control in ['natural','shuffle']:
            mask = evaluation.control == control
            evaluation.loc[mask,'q'] = bh(evaluation.loc[mask,'p'].values)
        evaluation.to_csv(out/'development_enrichment.tsv',sep='\t',index=False)
        hits = pd.DataFrame(hits)
        hits.to_csv(out/'development_scans.tsv',sep='\t',index=False)
        save(out/'thresholds.json',{'calibration_split':'discovery_natural','comparison':'>',
                         'target_sequence_fpr':cfg['development_control_fpr'],'thresholds':thresholds,
                         'scan':'maximum zero-order log-odds over both strands; not FIMO site p/q'})
        for name,pwm in retained.items():
            for other,m in traditional.items():
                sim, offset, strand = similarity(pwm,m)
                comparisons.append({'embedding':name,'streme':other,'cosine_similarity':sim,'offset':offset,'strand':strand,
                                    'interpretation':'descriptive alignment only; no significance'})
        pd.DataFrame(comparisons).to_csv(out/'streme_pwm_similarity.tsv',sep='\t',index=False)
        grammar = []
        pos_hits = hits[(hits.control == 'positive') & hits.present & hits.motif.isin(retained)]
        for seqid, g in pos_hits.groupby('sequence_id'):
            for (_,a),(_,b) in itertools.combinations(g.sort_values(['start1','motif']).iterrows(),2):
                grammar.append({'sequence_id':seqid,'upstream':a.motif,'downstream':b.motif,
                                'gap':int(b.start1-a.stop1-1),'center_distance':float(b.center_tss-a.center_tss),
                                'up_strand':a.strand,'down_strand':b.strand})
        pd.DataFrame(grammar,columns=['sequence_id','upstream','downstream','gap','center_distance','up_strand','down_strand']).to_csv(out/'development_grammar_descriptive.tsv',sep='\t',index=False)
        manifest.update(status='COMPLETED',elapsed_seconds=time.time()-started,
                        discovery_n=len(discovery),development_n=len(dev),motifs=len(retained),
                        pca_explained_variance=float(pca.explained_variance_ratio_.sum()),
                        input_hashes={str(p.relative_to(ROOT)):sha(p) for p in
                           [DATA/source/'main/sequence_master.tsv',DATA/source/'main/development.dinucleotide_null.fasta',traditional_path]},
                        output_sha256={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!='manifest.json'})
    except Exception:
        manifest.update(status='FAILED',elapsed_seconds=time.time()-started,error=traceback.format_exc())
        save(out/'manifest.json',manifest)
        raise
    save(out/'manifest.json',manifest)
    print(json.dumps({'run':tag,'status':manifest['status'],'seconds':manifest['elapsed_seconds'],'motifs':manifest['motifs']}),flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--sources',nargs='+',default=['tjupan_bacillus_subtilis'])
    parser.add_argument('--seeds',nargs='+',type=int,default=[20261005])
    parser.add_argument('--representations',nargs='+',choices=['prokbert','onehot','multiscale','tokenwindow'],default=['prokbert','onehot'])
    parser.add_argument('--layer',type=int,choices=[3,6])
    parser.add_argument('--width',type=int,choices=[6,8,10,12,16])
    args = parser.parse_args()
    cfg = json.loads((BRANCH/'config.json').read_text(encoding='utf-8'))
    assert json.loads((BRANCH/'audit/audit.json').read_text(encoding='utf-8'))['status']=='COMPLETED'
    model,tokenizer = load_model(cfg) if set(args.representations)&{'prokbert','multiscale','tokenwindow'} else (None,None)
    for source in args.sources:
        for representation in args.representations:
            for seed in args.seeds:
                run(source,seed,representation,model,tokenizer,args.layer,args.width)


if __name__ == '__main__':
    main()
