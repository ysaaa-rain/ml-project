"""Read-only evidence audit; writes only a new second-round verification file."""
import json
import numpy as np
import pandas as pd
from ..audit import BRANCH,ROOT,load,sha,save
from ..discovery import read_meme
from .common import NEXT

def main():
    result={'status':'COMPLETED','holdout_used':False,'hashes_checked':0,
            'historical_runs':0,'new_real_runs':0,'new_real_instances':0,
            'synthetic_conditions':0,'synthetic_fits':0,'synthetic_instances':0,
            'synthetic_truth_sites':0,'exact_cross_split_duplicates':0}
    for path in (BRANCH/'runs').glob('*/manifest.json'):
        m=json.loads(path.read_text(encoding='utf-8'))
        assert m['status']=='COMPLETED' and m['holdout_used'] is False
        for name,expected in m['output_sha256'].items():
            assert sha(path.parent/name)==expected,(path,name)
            result['hashes_checked']+=1
        result['historical_runs']+=1
    manifests=[]
    for path in (NEXT/'run_manifests').glob('*.json'):
        m=json.loads(path.read_text(encoding='utf-8'))
        assert m['status']=='COMPLETED' and m['holdout_used'] is False,path
        for name,expected in m.get('input_sha256',{}).items():
            assert sha(ROOT/name)==expected,(path,'input',name)
            result['hashes_checked']+=1
        for name,expected in m.get('output_sha256',{}).items():
            assert sha(ROOT/name)==expected,(path,'output',name)
            result['hashes_checked']+=1
        manifests.append({'run_id':path.stem,'status':m['status'],'seconds':m.get('elapsed_seconds'),
                          'process_peak_memory_bytes':m.get('process_peak_memory_bytes'),'manifest_sha256':sha(path)})
    for path in (NEXT/'real_data_comparison').glob('*__*/instances.tsv'):
        source=path.parent.name.split('__')[0]
        parent=load(source,'discovery').set_index('sequence_id')
        d=pd.read_csv(path,sep='\t')
        assert not d.duplicated(['motif','sequence_id']).any(),path
        for row in d.itertuples():
            p=parent.loc[row.sequence_id]
            assert row.sequence==p.sequence[row.start0:row.end0_exclusive]
            assert row.start1==row.start0+1 and row.stop1==row.end0_exclusive
            assert row.center_tss==(row.start1+row.stop1)/2-61
            assert row.leakage_group==p.leakage_group and row.split=='discovery'
        for pwm in read_meme(path.parent/'motifs.meme').values():
            assert np.all(pwm>0) and np.allclose(pwm.sum(axis=1),1)
        result['new_real_runs']+=1;result['new_real_instances']+=len(d)
    for path in (NEXT/'synthetic_benchmark').glob('*/data.json'):
        data=json.loads(path.read_text(encoding='utf-8'))
        groups=[data[x] for x in ['discovery','development','test','discovery_control','development_control','test_control','calibration_control']]
        for i,a in enumerate(groups):
            assert all(len(s)==81 and set(s)<=set('ACGT') for s in a)
            for b in groups[i+1:]:
                result['exact_cross_split_duplicates']+=len(set(a)&set(b))
        for split in ['discovery','development','test']:
            for site in data[split+'_truth']:
                assert data[split][site['sequence_index']][site['start0']:site['end0']]==site['planted_sequence']
                result['synthetic_truth_sites']+=1
        for inst in path.parent.glob('*__w*/instances.tsv'):
            d=pd.read_csv(inst,sep='\t')
            assert not d.duplicated(['motif','sequence_index']).any()
            for row in d.itertuples():
                assert row.fragment==data['discovery'][row.sequence_index][row.start0:row.end0]
            result['synthetic_instances']+=len(d);result['synthetic_fits']+=1
            matches=pd.read_csv(inst.parent/'strict/reference_matches.tsv',sep='\t')
            assert matches.candidate.is_unique and matches.family.is_unique
        result['synthetic_conditions']+=1
    metrics=pd.read_csv(NEXT/'synthetic_benchmark/strict_metrics.tsv',sep='\t')
    assert ((metrics[['precision','recall','f1']]>=0)&(metrics[['precision','recall','f1']]<=1)).all().all()
    assert (metrics.true_positives<=metrics.truth_count).all()
    assert not metrics.duplicated(['scenario','seed','method','width','split','family','tolerance_bp']).any()
    assert result['exact_cross_split_duplicates']==0
    assert result['historical_runs']==53 and result['new_real_runs']==24
    assert result['synthetic_conditions']==27 and result['synthetic_fits']==132
    result['limits']=['hash audit verifies frozen inputs and outputs, not historical execution independently',
                      'model fitting and truth separation reviewed in benchmark.py; no exhaustive homology audit rerun',
                      'metrics are reference-assisted localization; not unsupervised biological discovery precision']
    pd.DataFrame(manifests).to_csv(NEXT/'run_manifests/run_inventory_final.tsv',sep='\t',index=False)
    save(NEXT/'verification.json',result)
    print(json.dumps(result,ensure_ascii=False))

if __name__=='__main__':main()
