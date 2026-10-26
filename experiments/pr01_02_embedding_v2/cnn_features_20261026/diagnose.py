"""Post-score descriptive boundary diagnostics; never changes scores or models."""
import numpy as np
import pandas as pd

from .pipeline import ALPHABET,ART,HERE,PROTOCOL,SUMMARY,assert_locked,read,read_meme,save,sha
from ..discovery import similarity


def main():
    assert_locked();cfg=read(PROTOCOL)
    metrics=pd.read_csv(SUMMARY/'metrics.tsv',sep='\t')
    metrics=metrics[(metrics.split=='test')&(metrics.tolerance_bp==2)&(metrics.view=='CAP5')&(metrics.representation=='CNN_LABEL')]
    rows=[]
    for folder in sorted(ART.iterdir()):
        data=read(folder/'data.json');mutation=data['scenario']['mutation'];reference=np.full((cfg['width'],4),mutation/3)
        for i,c in enumerate(cfg['synthetic_pattern']):reference[i,ALPHABET.index(c)]=1-mutation
        for parent in sorted(folder.glob('a*')):
            seed=int(parent.name[1:]);out=parent/'CNN_LABEL';motifs=read_meme(out/'motifs.meme')
            names=[v['candidates'] for v in read(parent/'views.json') if v['representation']=='CNN_LABEL' and v['view']=='CAP5'][0]
            frame=pd.read_csv(out/'instances.tsv',sep='\t')
            candidates=[(similarity(motifs[name],reference,minimum_overlap=8),name) for name in names]
            match,winner=max(candidates,key=lambda item:(item[0][0],item[1])) if candidates else ((0.,0,'+'),None)
            original=metrics[(metrics.scenario==data['scenario']['id'])&(metrics.data_seed==data['seed'])&(metrics.algorithm_seed==seed)].iloc[0]
            selected=frame[frame.motif==winner];errors=[];directions=[]
            for r in selected.itertuples():
                truth=data['discovery_truth'][r.sequence_index];errors.append(abs(r.start0-truth['start0']))
                directions.append((r.strand==truth['strand'])==(match[2]=='+'))
            rows.append({'scenario':data['scenario']['id'],'data_seed':data['seed'],'algorithm_seed':seed,
                'primary_full_width_candidate':original.candidate,'primary_full_width_cosine':original.pwm_cosine,
                'primary_reference_eligible':bool(original.reference_eligible),'primary_f1':original.f1,
                'exploratory_partial_candidate':winner,'partial_cosine':match[0],'partial_offset':match[1],
                'overlap_bp':cfg['width']-abs(match[1]),'partial_orientation':match[2],
                'discovery_median_abs_start_error':float(np.median(errors)) if errors else None,
                'discovery_direction_agreement':float(np.mean(directions)) if directions else None})
    pd.DataFrame(rows).to_csv(SUMMARY/'boundary_diagnostics.tsv',sep='\t',index=False)
    save(SUMMARY/'diagnostics_policy.json',{'phase':'post-Test descriptive analysis, not preregistered endpoint',
        'minimum_overlap':8,'reference_used':True,'source_sha256':sha(HERE/'diagnose.py'),
        'selection':'only already-locked CNN_LABEL CAP5 candidates; no model/PWM/score/budget changes',
        'limits':'partial matches and Discovery offsets are diagnostic, not successful full-width recovery or validation'})
    print('Post-score boundary diagnostics saved',flush=True)


if __name__=='__main__':main()
