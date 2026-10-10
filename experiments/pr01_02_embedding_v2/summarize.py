"""Aggregate only completed, hash-verified runs; no planned results."""
import json
from pathlib import Path
import itertools
import numpy as np
import pandas as pd
from scipy.optimize import linear_sum_assignment
from .audit import BRANCH, load, save, sha
from .discovery import read_meme, similarity, bh


def main():
    out = BRANCH/'summary'
    out.mkdir(exist_ok=True)
    runs, evaluations, stability, gc_rows = [], [], [], []
    completed = []
    for path in sorted((BRANCH/'runs').iterdir()):
        manifest = json.loads((path/'manifest.json').read_text(encoding='utf-8'))
        if manifest['status'] != 'COMPLETED':
            continue
        for name,expected in manifest['output_sha256'].items():
            assert sha(path/name) == expected, path/name
        cfg = manifest['config']
        row = {'run':path.name,'source':manifest['source'],'method':manifest['representation'],
               'layer':cfg['baseline_layer'],'width':cfg['baseline_window'],'seed':cfg['seed'],
               'seconds':manifest['elapsed_seconds'],'motifs':manifest['motifs'],
               'discovery_n':manifest['discovery_n'],'development_n':manifest['development_n']}
        runs.append(row)
        completed.append((path,row))
        eval = pd.read_csv(path/'development_enrichment.tsv',sep='\t')
        for key in ['run','source','layer','width','seed']:
            eval[key] = row[key]
        evaluations.append(eval)
        if row['layer']==6 and row['width']==10 and row['seed']==20261005:
            scans = pd.read_csv(path/'development_scans.tsv',sep='\t')
            for label,control in [(1,'positive'),(0,'natural')]:
                frame = load(row['source'],'development',label)
                frame['gc'] = frame.sequence.str.count('[GC]')/81
                subset = scans[(scans.control==control)&~scans.motif.str.startswith('STREME_')].merge(frame[['sequence_id','gc']],on='sequence_id')
                subset['gc_bin'] = pd.cut(subset.gc,bins=[0,.2,.4,.6,.8,1],include_lowest=True).astype(str)
                for (motif,bin),g in subset.groupby(['motif','gc_bin']):
                    gc_rows.append({'source':row['source'],'method':row['method'],'motif':motif,
                                    'control':control,'gc_bin':bin,'n_sequences':len(g),
                                    'hit_rate':float(g.present.mean()),'scope':'descriptive composition audit; no matching/resplitting'})
    all_eval = pd.concat(evaluations,ignore_index=True)
    # Main frozen configuration across six species, per seed; methods together with one copy of STREME.
    all_eval['q_six_species'] = np.nan
    for seed in sorted(set(r['seed'] for r in runs)):
        base = all_eval[(all_eval.layer==6)&(all_eval.width==10)&(all_eval.seed==seed)&all_eval.method.isin(['prokbert','onehot','tokenwindow','STREME'])]
        for control in ['natural','shuffle']:
            family = base[base.control==control].drop_duplicates(['source','method','motif'])
            q = bh(family.p.to_numpy())
            for (_,r),value in zip(family.iterrows(),q):
                mask = (all_eval.source==r.source)&(all_eval.method==r.method)&(all_eval.motif==r.motif)&(all_eval.seed==seed)&(all_eval.control==control)&(all_eval.layer==6)&(all_eval.width==10)
                all_eval.loc[mask,'q_six_species'] = value
    all_eval.to_csv(out/'all_development_enrichment.tsv',sep='\t',index=False)
    pd.DataFrame(runs).to_csv(out/'run_summary.tsv',sep='\t',index=False)
    pd.DataFrame(gc_rows).to_csv(out/'gc_composition_diagnostic.tsv',sep='\t',index=False)
    # Stability is a discovery-only one-to-one PWM match, never duplicate biological evidence.
    keys = {(r['source'],r['method'],r['layer'],r['width']) for r in runs}
    for source,method,layer,width in sorted(keys):
        selected = [(p,r) for p,r in completed if (r['source'],r['method'],r['layer'],r['width'])==(source,method,layer,width)]
        for (pa,ra),(pb,rb) in itertools.combinations(selected,2):
            a,b = read_meme(pa/'motifs.meme'),read_meme(pb/'motifs.meme')
            an,bn = list(a),list(b)
            matrix = np.array([[similarity(a[x],b[y],minimum_overlap=width)[0] for y in bn] for x in an])
            ai,bi = linear_sum_assignment(-matrix)
            ia = pd.read_csv(pa/'instances.tsv',sep='\t')
            ib = pd.read_csv(pb/'instances.tsv',sep='\t')
            for i,j in zip(ai,bi):
                sa = set(zip(ia[ia.motif==an[i]].sequence_id,ia[ia.motif==an[i]].start0))
                sb = set(zip(ib[ib.motif==bn[j]].sequence_id,ib[ib.motif==bn[j]].start0))
                stability.append({'source':source,'method':method,'layer':layer,'width':width,
                    'seed_a':ra['seed'],'seed_b':rb['seed'],'motif_a':an[i],'motif_b':bn[j],
                    'full_width_cosine':matrix[i,j], 'exact_instance_jaccard':len(sa&sb)/len(sa|sb),
                    'stable_pwm':bool(matrix[i,j]>=.8)})
    pd.DataFrame(stability).to_csv(out/'seed_stability.tsv',sep='\t',index=False)
    baseline = all_eval[(all_eval.layer==6)&(all_eval.width==10)&(all_eval.seed==20261005)].drop_duplicates(['source','method','motif','control'])
    counts = baseline.groupby(['source','method','control']).agg(candidates=('motif','size'),
                          bh_significant=('q_six_species',lambda v:int((v<=.05).sum()))).reset_index()
    counts.to_csv(out/'baseline_significant_counts.tsv',sep='\t',index=False)
    resource = {'completed_runs':len(runs),'total_serial_run_seconds':sum(r['seconds'] for r in runs),
                'branch_disk_bytes':sum(p.stat().st_size for p in BRANCH.rglob('*') if p.is_file()),
                'cuda':False,'gpu_memory':'not applicable','peak_cpu_memory':'not measured',
                'timing_limit':'concurrent CPU jobs; sums are not wall time; cached repeat seeds are cheaper'}
    save(out/'resources.json',resource)
    print(counts.to_string(index=False))
    print(json.dumps(resource))


if __name__ == '__main__':
    main()
