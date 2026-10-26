"""Audit control weights/labels and independently reconstruct OH48 retained instances."""
import hashlib

import numpy as np
import pandas as pd
import torch

from .pipeline import (ART,HERE,PROTOCOL,SUMMARY,LocalCNN,assert_locked,background,
    make_pwm,raw_instances,read,reject_background,save,select_instances,sequence_windows,setup,sha)


def main():
    cfg=read(PROTOCOL);setup(cfg);assert_locked();randoms=0;labels=0;oh_pools=0
    for folder in sorted(ART.iterdir()):
        data=read(folder/'data.json')
        for parent in sorted(folder.glob('a*')):
            seed=int(parent.name[1:]);torch.manual_seed(seed)
            expected=LocalCNN(cfg['cnn']['filters'],cfg['width']).state_dict()
            random=torch.load(parent/'CNN_RANDOM/model.pt',weights_only=True)
            assert expected.keys()==random.keys() and all(torch.equal(expected[k],random[k]) for k in expected);randoms+=1
            original=np.array([1]*len(data['discovery'])+[0]*len(data['discovery_control']),dtype=np.float32)
            shuffled=original[np.random.default_rng(seed).permutation(len(original))]
            for method,y in [('CNN_LABEL',original),('CNN_SHUFFLED',shuffled)]:
                training=read(parent/method/'training.json')
                assert hashlib.sha256(y.tobytes()).hexdigest()==training['labels_sha256']
                weights=torch.load(parent/method/'model.pt',weights_only=True)
                assert any(not torch.equal(weights[k],random[k]) for k in weights);labels+=1
            out=parent/'OH48_BG';assignments=np.load(out/'assignments.npz')
            selections,_=select_instances(assignments['labels'],assignments['distance'],82-cfg['width'],1.0)
            frame=pd.read_csv(out/'instances.tsv',sep='\t');bg=background(data['discovery_control'])
            blank=sequence_windows(data['discovery_control'],cfg['width']);expected_rows=[]
            for label,selected in selections.items():
                rows=raw_instances(data['discovery'],selected,assignments['distance'],cfg['width'])
                if len(rows)<cfg['minimum_promoters']:continue
                _,pwm=make_pwm([r['fragment'] for r in rows]);kept,_=reject_background(rows,pwm,blank,cfg,bg)
                if len(kept)>=cfg['minimum_promoters']:
                    expected_rows.extend((f'C{label:02d}',r['sequence_index'],r['start0'],r['fragment']) for r in kept)
            actual=[(r.motif,r.sequence_index,r.start0,r.fragment) for r in frame.itertuples()]
            assert sorted(actual)==sorted(expected_rows);oh_pools+=1
    result={'status':'PASSED','source_sha256':sha(HERE/'qa.py'),'random_initial_states_matched':randoms,
        'label_streams_and_trained_weight_changes_checked':labels,'OH48_pools_reconstructed':oh_pools,
        'real_holdout_used':False,'limits':'saved-state consistency; no independent full retraining or near-homology audit'}
    save(SUMMARY/'independent_qa.json',result);print(result,flush=True)


if __name__=='__main__':main()
