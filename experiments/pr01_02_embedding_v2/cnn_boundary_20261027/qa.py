"""Additional saved-state controls, OH48 reconstruction and boundary accounting."""
import hashlib
import joblib
import numpy as np
import pandas as pd
import torch

from .pipeline import (ART,HERE,PROTOCOL,SUMMARY,LocalCNN,assert_locked,background,make_pwm,read,
    save,sequence_windows,setup,sha)
from ..cnn_features_20261026.pipeline import raw_instances,reject_background,select_instances


def main():
    cfg=read(PROTOCOL);setup(cfg);assert_locked();initials=labels=oh_pools=0;records=[]
    for folder in sorted(ART.iterdir()):
        data=read(folder/'data.json')
        for parent in sorted(folder.glob('a*')):
            seed=int(parent.name[1:]);torch.manual_seed(seed)
            initial=LocalCNN(cfg['cnn']['filters'],cfg['width'])
            initial_hash=hashlib.sha256(b''.join(p.detach().numpy().tobytes() for p in initial.parameters())).hexdigest()
            info=read(parent/'model_info.json');assert info['initial_state_sha256']==initial_hash;initials+=1
            y=np.array([1]*len(data['discovery'])+[0]*len(data['discovery_control']),dtype=np.float32)
            expected=hashlib.sha256(y.tobytes()).hexdigest()
            for direction in ['TIED','FWD']:
                assert info[direction+'_labels_sha256']==expected
                state=torch.load(parent/f'{direction}_model.pt',weights_only=True)
                assert any(not torch.equal(state[k],v) for k,v in initial.state_dict().items());labels+=1
                maxima=np.load(parent/f'{direction}_maxima.npz')
                for j,row in enumerate(read(parent/f'{direction}_decisions.json')):
                    assert row['threshold']==float(maxima['threshold'][j])
                    records.append({'scenario':data['scenario']['id'],'data_seed':data['seed'],'algorithm_seed':seed,
                        'direction':direction,**{k:v for k,v in row.items() if k!='scores'}})
            out=parent/'OH48_BG';assignments=np.load(out/'assignments.npz')
            assert joblib.load(out/'kmeans.joblib').n_features_in_==48
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
    pd.DataFrame(records).to_csv(SUMMARY/'filter_decisions.tsv',sep='\t',index=False)
    result={'status':'PASSED','source_sha256':sha(HERE/'qa.py'),'paired_initial_states_recomputed':initials,
        'label_streams_and_weight_changes_checked':labels,'OH48_pools_reconstructed':oh_pools,
        'real_holdout_used':False,'limits':'saved-state consistency, no independent full retraining or near-homology audit'}
    save(SUMMARY/'independent_qa.json',result);print(result,flush=True)


if __name__=='__main__':main()
