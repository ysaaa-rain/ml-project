"""Recompute persisted ranking inputs and audit new synthetic partitions."""
import pandas as pd
import joblib

from .pipeline import (ART,HERE,PROTOCOL,SUMMARY,assert_locked,rank_candidates,read,
                       read_meme,save,sha)
from ..bottleneck_20261025.pilot import reverse_complement


def main():
    assert_locked();cfg=read(PROTOCOL);truths=0;sequences=0;rankings=0;dimensions={r['id']:set() for r in cfg['representations']}
    for folder in ART.iterdir():
        data=read(folder/'data.json');seen=set()
        for split in ['discovery','discovery_control','development','development_control','test','test_control','calibration_control']:
            canonical={min(s,reverse_complement(s)) for s in data[split]}
            assert len(canonical)==len(data[split]) and not seen & canonical
            seen.update(canonical);sequences+=len(data[split])
        for split in ['discovery','development','test']:
            for r in data[split+'_truth']:
                assert data[split][r['sequence_index']][r['start0']:r['end0']]==r['fragment'];truths+=1
        for representation in cfg['representations']:
            for algorithm_seed in cfg['algorithm_seeds']:
                parent=folder/f'{representation["id"]}__a{algorithm_seed}'
                fitted=joblib.load(parent/'kmeans.joblib')
                expected=48 if representation['pca'] is None else representation['pca']
                assert fitted.n_features_in_==expected;dimensions[representation['id']].add(int(fitted.n_features_in_))
                for variant in cfg['variants']:
                    out=folder/f'{representation["id"]}__a{algorithm_seed}'/variant
                    motifs=read_meme(out/'motifs.meme');frame=pd.read_csv(out/'instances.tsv',sep='\t')
                    support=frame.groupby('motif').size().to_dict();saved=read(out/'ranking.json')
                    rebuilt=rank_candidates(motifs,support,data['discovery'],data['discovery_control'],cfg)
                    assert saved==rebuilt and len(saved)==len({r['candidate'] for r in saved});rankings+=1
        for path in folder.glob('views__*.json'):
            for view in read(path):assert len(view['candidates'])==len(set(view['candidates']))
    result={'status':'PASSED','source_sha256':sha(HERE/'qa.py'),'truth_sites_checked':truths,
            'synthetic_sequences_checked':sequences,'discovery_rankings_recomputed':rankings,
            'exact_RC_partition_duplicates':0,'real_holdout_used':False,
            'cluster_input_dimensions':{name:sorted(values) for name,values in dimensions.items()},
            'duplicate_scope':'within each condition; blank streams can repeat across difficulty conditions sharing a seed',
            'limits':'persisted-artifact consistency; no independent encoder/clustering refit or near-homology audit'}
    save(SUMMARY/'independent_qa.json',result);print(result,flush=True)


if __name__=='__main__':main()
