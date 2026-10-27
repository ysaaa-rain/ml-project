"""Verify provenance, equal budgets and published summaries."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
from .pipeline import ROOT, HERE, ART, BASE, sha

SYN=ROOT/'reports/same_batch_motif_equal_20261027'
REAL=ROOT/'reports/same_batch_real_equal_20261027'

def read(p):return json.loads(p.read_text(encoding='utf8'))

def main():
    lock=read(HERE/'fit_lock.json')
    assert lock['cnn_fits']==45 and len(list(ART.glob('*/*/model.pt')))==45
    assert lock['pipeline_sha256']==sha(HERE/'pipeline.py') and lock['protocol_sha256']==sha(HERE/'protocol.json')
    assert lock['prior_lock_sha256']==sha(BASE/'summary/fit_lock.json')
    for name,digest in lock['source_sha256'].items():assert sha(ROOT/name)==digest,name
    for name,digest in lock['output_sha256'].items():assert sha(HERE/name)==digest,name
    sm=read(SYN/'manifest.json');rm=read(REAL/'manifest.json')
    assert sm['fit_lock_sha256']==sha(HERE/'fit_lock.json')
    assert sm['equal_budget_script_sha256']==sha(HERE/'equal_budget.py')
    assert sm['report_source_sha256']==sha(HERE/'report.py')
    assert rm['script_sha256']==sha(HERE/'real_compare.py') and rm['report_script_sha256']==sha(HERE/'real_report.py')
    for directory,manifest,key in [(SYN,sm,'outputs_sha256'),(REAL,rm,'output_sha256')]:
        for name,digest in manifest[key].items():assert sha(directory/name)==digest,name
    for name,digest in rm['source_sha256'].items():assert sha(ROOT/name)==digest,name
    syn=pd.read_csv(SYN/'runs.tsv',sep='\t');real=pd.read_csv(REAL/'runs.tsv',sep='\t')
    assert len(syn)==450 and len(real)==42
    keys=['scenario','data_seed','algorithm_seed','split']
    assert (syn.groupby(keys).method.nunique()==5).all()
    assert (syn.groupby(keys).common_budget.nunique()==1).all()
    assert (syn.candidate_count==syn.common_budget).all()
    assert (real.groupby('source').common_budget.nunique()==1).all()
    assert (real.candidate_count==real.common_budget).all()
    assert np.allclose(real.positive_rate,real.positive_hits/real.positive_n)
    assert np.allclose(real.control_rate,real.control_hits/real.control_n)
    assert np.allclose(real.rate_difference,real.positive_rate-real.control_rate)
    old=pd.read_csv(ROOT/'reports/practice_motif_20261027/synthetic_partial_metrics.tsv',sep='\t')
    old=old[(old.version=='budget_validation_20261026')&(old.variant=='BG_ONLY')&(old.view=='CAP5')].rename(columns={'representation':'method'})
    same=syn[(syn.method.isin(['EMB128','OH24','OH48']))&(syn.common_budget==5)].merge(old,on=keys+['method'],suffixes=('_new','_old'))
    assert len(same)>0 and np.allclose(same.f1_new,same.f1_old,atol=1e-12)
    summary={'status':'PASSED','cnn_models':45,'synthetic_rows':len(syn),'real_runs':len(real),
      'old_budget5_scores_reproduced':len(same),'synthetic_budget_counts':syn[syn.method=='EMB128'].groupby('common_budget').size().to_dict(),
      'real_species':int(real.source.nunique()),'source_hashes_checked':len(lock['source_sha256'])+len(rm['source_sha256']),
      'cnn_artifact_hashes_checked':len(lock['output_sha256'])}
    print(json.dumps(summary,ensure_ascii=False,indent=2))
    (SYN/'verification.json').write_text(json.dumps(summary,indent=2)+'\n',encoding='utf8')

if __name__=='__main__':main()
