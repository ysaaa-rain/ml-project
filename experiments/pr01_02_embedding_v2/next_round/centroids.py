import pandas as pd
import numpy as np
import joblib
from ..audit import BRANCH
from .common import NEXT,start,finish


def main():
    target,manifest=start('centroid_diagnostic',{'space':'PCA inverse transform into original feature coordinates'},[NEXT/'stability_diagnostics/pwm_family_matches.tsv'])
    d=pd.read_csv(NEXT/'stability_diagnostics/pwm_family_matches.tsv',sep='\t')
    d=d[d.strategy=='full'];cache={};rows=[]
    for r in d.itertuples():
        centers=[]
        for seed,name in [(r.seed_a,r.motif_a),(r.seed_b,r.motif_b)]:
            key=(r.source,r.method,seed)
            if key not in cache:
                obj=joblib.load(BRANCH/'runs'/f'{r.source}__{r.method}__L6W10__s{seed}'/'fitted.joblib')
                cache[key]=(obj['pca'].inverse_transform(obj['kmeans'].cluster_centers_)-obj['pca'].mean_)
            centers.append(cache[key][int(name.rsplit('C',1)[-1])])
        a,b=centers
        rows.append({'source':r.source,'method':r.method,'seed_a':r.seed_a,'seed_b':r.seed_b,
                     'motif_a':r.motif_a,'motif_b':r.motif_b,'pwm_similarity':r.similarity,'pwm_orientation':r.orientation,
                     'center_cosine_original_feature_space':float(a@b/(np.linalg.norm(a)*np.linalg.norm(b))),
                     'center_l2_original_feature_space':float(np.linalg.norm(a-b)),
                     'limit':'native feature centers, no assumed hidden-vector RC transform; descriptive only'})
    out=NEXT/'stability_diagnostics/cluster_centroid_comparison.tsv'
    pd.DataFrame(rows).to_csv(out,sep='\t',index=False)
    finish(target,manifest,[out])


if __name__=='__main__':main()
