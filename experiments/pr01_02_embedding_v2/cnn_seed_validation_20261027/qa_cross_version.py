"""Audit exact/RC sequences across the new and three preceding synthetic versions."""
from .pipeline import ART,HERE,PREVIOUS,SUMMARY,assert_locked,read,reverse_complement,save,sha

SPLITS=['discovery','discovery_control','development','development_control','test','test_control','calibration_control']


def canonical_sequences(folder):
    data=read(folder/'data.json');result=set()
    for split in SPLITS:
        result.update(min(s,reverse_complement(s)) for s in data[split])
    return result


def main():
    assert_locked();new=set()
    for folder in ART.iterdir():new.update(canonical_sequences(folder))
    old_versions=['cnn_boundary_20261027','cnn_features_20261026','budget_validation_20261026']
    rows=[]
    for version in old_versions:
        old=set();old_root=PREVIOUS.parent/version/'artifacts'
        folders=[p for p in old_root.iterdir() if p.is_dir() and (p/'data.json').is_file()]
        assert folders,version
        for folder in folders:old.update(canonical_sequences(folder))
        duplicates=len(new&old);rows.append({'version':version,'old_conditions':len(folders),
            'distinct_canonical_sequences':len(old),'exact_or_reverse_complement_overlaps':duplicates})
        assert duplicates==0,(version,duplicates)
    result={'status':'PASSED','source_sha256':sha(HERE/'qa_cross_version.py'),
        'new_distinct_canonical_sequences':len(new),'versions':rows,'real_holdout_used':False,
        'limits':'only exact/RC matches audited; not a near-homology audit or proof of biological independence'}
    save(SUMMARY/'cross_version_qa.json',result);print(result,flush=True)


if __name__=='__main__':main()
