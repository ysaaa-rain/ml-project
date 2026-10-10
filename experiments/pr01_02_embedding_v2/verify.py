"""Verify real outputs against frozen parent DNA and cached manifests."""
import json
import numpy as np
import pandas as pd
from .audit import BRANCH, load, sha, save
from .discovery import read_meme


def main():
    total = 0
    checked = []
    for path in sorted((BRANCH/'runs').iterdir()):
        manifest = json.loads((path/'manifest.json').read_text(encoding='utf-8'))
        assert manifest['status'] == 'COMPLETED',path
        assert manifest['holdout_used'] is False
        for name,expected in manifest['output_sha256'].items():
            assert sha(path/name) == expected, path/name
        parent = load(manifest['source'],'discovery').set_index('sequence_id')
        instances = pd.read_csv(path/'instances.tsv',sep='\t')
        assert not instances.duplicated(['motif','sequence_id']).any()
        for row in instances.itertuples():
            p = parent.loc[row.sequence_id]
            assert row.sequence == p.sequence[row.start0:row.end0_exclusive]
            assert row.start1 == row.start0+1 and row.stop1 == row.end0_exclusive
            assert row.center_tss == (row.start1+row.stop1)/2-61
            assert row.leakage_group == p.leakage_group
            assert row.split == 'discovery'
        for pwm in read_meme(path/'motifs.meme').values():
            assert np.all(pwm>0) and np.allclose(pwm.sum(axis=1),1)
        total += len(instances)
        checked.append(path.name)
    result = {'status':'COMPLETED','runs':len(checked),'instances_exactly_verified':total,
              'checks':['output hashes','frozen source DNA substrings','1based/0based/TSS coordinate equivalence',
                        'no duplicate promoter per cluster','group and split provenance','PWM normalization','no holdout evaluation'],
              'run_names':checked}
    save(BRANCH/'summary/verification.json',result)
    print(json.dumps({k:v for k,v in result.items() if k!='run_names'}))


if __name__=='__main__':
    main()
