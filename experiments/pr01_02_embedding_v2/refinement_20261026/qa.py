"""Independent checks from persisted instances, not fitting-time variables."""
import numpy as np
import pandas as pd

from .pipeline import (ART,HERE,PROTOCOL,SUMMARY,assert_locked,make_pwm,read,
                       read_meme,reverse_complement,save,sha)


def instance_keys(frame):
    return {(r.motif,r.sequence_index,r.original_start0,r.start0,r.end0,r.strand,r.shift,r.fragment)
            for r in frame.itertuples()}


def main():
    assert_locked();cfg=read(PROTOCOL);manifest=read(SUMMARY/'manifest.json')
    lock=read(SUMMARY/'fit_lock.json')
    assert manifest['source_sha256']==lock['source_sha256']
    assert manifest['protocol_sha256']==lock['protocol_sha256']
    pwms=0;instances=0;pairs=0
    for folder in ART.iterdir():
        data=read(folder/'data.json')
        for parent in [p for p in folder.iterdir() if p.is_dir()]:
            frames={}
            for variant in cfg['variants']:
                frame=pd.read_csv(parent/variant/'instances.tsv',sep='\t');frames[variant]=frame
                motifs=read_meme(parent/variant/'motifs.meme')
                assert set(frame.motif)==set(motifs)
                assert not frame.duplicated(['motif','sequence_index']).any()
                for row in frame.itertuples():
                    sequence=data['discovery'][row.sequence_index]
                    assert 0<=row.original_start0<=len(sequence)-cfg['width']
                    assert 0<=row.start0<row.end0<=len(sequence)
                    assert row.strand in ['+','-'] and row.end0-row.start0==cfg['width']
                    assert row.start0-row.original_start0==row.shift
                    if variant in ['RAW','BG_ONLY']:assert row.strand=='+' and row.shift==0
                    if variant=='RC_ONLY':assert row.shift==0
                    fragment=sequence[row.start0:row.end0]
                    if row.strand=='-':fragment=reverse_complement(fragment)
                    assert row.fragment==fragment;instances+=1
                for name,pwm in motifs.items():
                    fragments=frame[frame.motif==name].fragment.tolist()
                    assert len(fragments)>=cfg['minimum_promoters']
                    _,rebuilt=make_pwm(fragments)
                    assert np.allclose(rebuilt,pwm,rtol=0,atol=1e-8);pwms+=1
            for source,filtered in [('RAW','BG_ONLY'),('RC_SHIFT','RC_SHIFT_BG')]:
                assert instance_keys(frames[filtered])<=instance_keys(frames[source]);pairs+=1
    result={'status':'PASSED','source_sha256':sha(HERE/'qa.py'),'candidate_pwms_rebuilt':pwms,
            'instances_checked':instances,'background_subset_pairs':pairs,
            'bounds_orientation_unique_per_sequence':True,'real_holdout_used':False,
            'limits':'persisted-artifact consistency; not independent model refitting'}
    save(SUMMARY/'independent_qa.json',result);print(result,flush=True)


if __name__=='__main__':main()
