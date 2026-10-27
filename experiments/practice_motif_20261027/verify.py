"""Check candidate coverage, frozen memberships, hashes and optional browser UI."""
import json
from pathlib import Path

import numpy as np
import pandas as pd

from .build import BRANCH,HERE,OUT,ROOT,sha


def browser_check():
    try:
        from playwright.sync_api import sync_playwright
    except ImportError:return {'status':'UNAVAILABLE','reason':'playwright not installed'}
    paths=[Path('C:/Program Files/Google/Chrome/Application/chrome.exe'),
           Path('C:/Program Files (x86)/Microsoft/Edge/Application/msedge.exe')]
    executable=next((p for p in paths if p.is_file()),None)
    if executable is None:return {'status':'UNAVAILABLE','reason':'Chrome/Edge executable not found'}
    qa=ROOT/'tmp/practice_motif_20261027_qa';qa.mkdir(parents=True,exist_ok=True)
    with sync_playwright() as p:
        browser=p.chromium.launch(headless=True,executable_path=str(executable))
        page=browser.new_page(viewport={'width':1440,'height':900},device_scale_factor=1)
        errors=[];page.on('pageerror',lambda err:errors.append(str(err)))
        page.goto((OUT/'index.html').resolve().as_uri(),wait_until='load')
        assert page.locator('#body tr').count()==50 and '1829' in page.locator('#count').inner_text()
        page.locator('#next').click();assert '51–100' in page.locator('#count').inner_text()
        page.locator('#kind').select_option(label='传统方法');assert '540' in page.locator('#count').inner_text()
        page.locator('#kind').select_option(label='Embedding历史');assert '1289' in page.locator('#count').inner_text()
        page.locator('#method').select_option(label='prokbert');assert '436' in page.locator('#count').inner_text()
        assert 'L3 / W10 / s20261005' in page.locator('#body').inner_text()
        logo=page.get_by_role('link',name='Logo').first;assert logo.count()==1
        assert (OUT/logo.get_attribute('href')).resolve().is_file()
        page.locator('#search').fill('this-candidate-does-not-exist');assert '没有匹配结果' in page.locator('#count').inner_text()
        page.locator('#search').fill('');page.screenshot(path=str(qa/'desktop_final.png'),full_page=True)
        page.set_viewport_size({'width':390,'height':844});page.screenshot(path=str(qa/'mobile_final.png'),full_page=True)
        assert not errors,errors;browser.close()
    return {'status':'PASSED','browser':str(executable),'desktop_and_mobile_screenshots':True,
        'rows_all':1829,'traditional':540,'embedding':1289,'prokbert':436,'pagination_search_logo_layer_labels':'PASSED'}


def main():
    manifest=json.loads((OUT/'manifest.json').read_text(encoding='utf8'));assert manifest['status']=='COMPLETED'
    for name,digest in manifest['source_hashes'].items():assert sha(ROOT/name)==digest,name
    for name,digest in manifest['output_sha256'].items():assert sha(OUT/name)==digest,name
    trad=pd.read_csv(OUT/'traditional_candidates.tsv',sep='\t');embedding=pd.read_csv(OUT/'embedding_candidates.tsv',sep='\t')
    original=pd.read_csv(ROOT/'results/motif/m3_expanded_20261009/candidates.tsv',sep='\t')
    assert len(trad)==len(original)==540 and set(trad.motif_id)==set(original.motif_id)
    assert trad.motif_id.is_unique and trad.original_main17.sum()==17
    source=pd.read_csv(BRANCH/'summary/all_development_enrichment.tsv',sep='\t');keys=['run','motif','method','source','layer','width','seed']
    assert len(embedding)==1289 and not embedding.duplicated(keys).any()
    for control in ['natural','shuffle']:
        source_control=source[source.control==control]
        joined=source_control.merge(embedding,on=keys,validate='one_to_one')
        assert len(joined)==len(source_control)==1289
        for column in ['positive_groups','control_groups','positive_hits','control_hits','p','q','q_six_species']:
            assert np.allclose(joined[column],joined[f'{control}_{column}'],equal_nan=True),f'{control} {column}'
    metrics=pd.read_csv(OUT/'synthetic_partial_metrics.tsv',sep='\t');assert len(metrics)==4800
    assert metrics.version.nunique()==4 and (metrics.groupby('split').size()==2400).all()
    assert metrics.f1.between(0,1).all() and metrics.family_blank_fpr.between(0,1).all()
    nonempty=metrics[metrics.candidate_count>0];assert nonempty.overlap_bp.between(8,12).all()
    empty=metrics[metrics.candidate_count==0];assert (empty.f1==0).all()
    cache={};membership_checked=0
    for r in metrics.itertuples():
        base=BRANCH/r.version;folder=base/'artifacts'/f'{r.scenario}__s{r.data_seed}'
        if r.version=='budget_validation_20261026':path=folder/f'views__a{r.algorithm_seed}__{r.variant}.json'
        else:path=folder/f'a{r.algorithm_seed}'/'views.json'
        if path not in cache:cache[path]=json.loads(path.read_text(encoding='utf8'))
        views=[v for v in cache[path] if v['representation']==r.representation and v['view']==r.view]
        assert len(views)==1 and len(views[0]['candidates'])==r.candidate_count
        assert r.candidate_count==0 or r.candidate in views[0]['candidates'];membership_checked+=1
    for file in ['six_species_profiles.png','synthetic_partial_comparison.png']:
        assert (OUT/file).is_file() and (OUT/file).stat().st_size>10000
    result={'status':'PASSED','verify_source_sha256':sha(HERE/'verify.py'),
        'source_hashes_checked':len(manifest['source_hashes']),'output_hashes_checked':len(manifest['output_sha256']),
        'traditional_candidates_checked':len(trad),'embedding_rows_checked':len(embedding),
        'synthetic_split_rows_checked':len(metrics),'frozen_memberships_checked':membership_checked,
        'synthetic_overlap_lengths':{str(k):int(v) for k,v in metrics[metrics.split=='test'].overlap_bp.value_counts().sort_index().items()},
        'browser':browser_check(),'limits':'existing post-hoc synthetic evaluation and prior aggregate outputs; no real Holdout access or new biological confirmation'}
    (OUT/'verification.json').write_text(json.dumps(result,ensure_ascii=False,indent=2)+'\n',encoding='utf8')
    print(json.dumps(result,ensure_ascii=False),flush=True)


if __name__=='__main__':main()
