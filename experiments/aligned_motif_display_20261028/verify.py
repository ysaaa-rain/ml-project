"""Verify artifact provenance, unchanged history and common evaluation outputs."""
import numpy as np
import pandas as pd

from ..pr01_02_embedding_v2.discovery import max_scan, make_pwm
from .pipeline import (HERE, ROOT, OUT, BRANCH, DATA, read, save, sha, unchanged,
                       pools, load_split, background, rank, evaluate, sequence_windows, scan_windows)


def main():
    unchanged()
    manifest = read(OUT / 'manifest.json')
    for name, digest in manifest['outputs_sha256'].items():
        assert sha(OUT / name) == digest, name
    assert sha(HERE / 'report.py') == manifest['report_source_sha256']
    lock = read(OUT / 'fit_lock.json')
    assert lock['ranking_sha256'] == sha(OUT / 'discovery_rankings.tsv')
    for name, digest in lock['source_sha256'].items():
        assert sha(ROOT / name) == digest, name
    assert not lock['holdout_scored'] and not manifest['holdout_scored']
    pos, neg, _, mixed = load_split('discovery')
    dp, dn, shuffled, dmixed = load_split('development')
    bg = background(neg.sequence.tolist())
    family = pd.read_csv(OUT / 'family.tsv', sep='\t', dtype={'seed': str})
    ranking = pd.read_csv(OUT / 'discovery_rankings.tsv', sep='\t', dtype={'seed': str})
    candidates = pd.read_csv(OUT / 'candidates.tsv', sep='\t', dtype={'seed': str})
    positions = pd.read_csv(OUT / 'positions.tsv', sep='\t', dtype={'seed': str})
    old = pd.read_csv(BRANCH / 'real_bsub_comparison_20261027/summary/family.tsv', sep='\t', dtype={'seed': str})
    recomputed = 0
    scalar_checks = 0
    for method, seed, purpose, _, motifs in pools():
        computed = rank(motifs, pos, neg, bg, .9, mixed)
        table = ranking[(ranking.method == method) & (ranking.seed == seed)].sort_values('rank')
        assert table.candidate.tolist() == [r['candidate'] for r in computed]
        assert np.allclose(table.threshold, [r['threshold'] for r in computed])
        assert np.allclose(table.discovery_score, [r['discovery_score'] for r in computed])
        entries = [dict(row, rank=i) for i, row in enumerate(computed, 1)]
        cr, fr, loc = evaluate(motifs, entries, dp, dn, shuffled, dmixed, bg, positions=purpose == 'primary')
        actual_candidates = candidates[(candidates.method == method) & (candidates.seed == seed)]
        for field in ['positive_rate', 'natural_rate', 'shuffle_rate']:
            assert np.allclose(actual_candidates[field], [r[field] for r in cr])
        for expected in fr:
            row = family[(family.method == method) & (family.seed == seed) & (family.view == expected['view'])].iloc[0]
            for field, value in expected.items():
                if field == 'view':
                    assert row[field] == value
                else:
                    assert np.isclose(row[field], value), (method, seed, field)
            recomputed += 1
        if purpose == 'primary':
            actual = positions[(positions.method == method) & (positions.seed == seed)]
            assert actual[['candidate', 'rank', 'class', 'center_tss', 'groups']].to_dict('records') == loc
        # Independent per-sequence scalar scan also checks strand and coordinate mapping.
        sample = dp.sequence.head(8).tolist() + dn.sequence.head(8).tolist()
        for item in computed[:2]:
            pwm = motifs[item['candidate']]
            expected = max_scan(sample, pwm, bg)
            score, start, forward = scan_windows(sequence_windows(sample, len(pwm)), pwm, bg)
            assert np.allclose(score, [x[0] for x in expected])
            assert np.array_equal(start, [x[1] for x in expected])
            assert np.array_equal(forward, [x[2] == '+' for x in expected])
            scalar_checks += len(sample)
    # Existing seven pools must give exactly the previously published common-scan results.
    fields = ['candidates', 'positive_groups', 'natural_groups', 'positive_rate', 'natural_rate',
              'shuffle_rate', 'natural_difference', 'paired_positive_only', 'paired_shuffle_only', 'paired_groups']
    matched = old.merge(family, on=['method', 'seed', 'view'], suffixes=('_old', '_new'))
    assert len(matched) == len(old) == 14
    for field in fields:
        assert np.allclose(matched[field + '_old'], matched[field + '_new']), field
    # Synthetic coordinate fixture is a scanner check only, not an extra experiment.
    pattern = 'AACGTTTA'
    reverse = pattern.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    _, pwm = make_pwm([pattern] * 20)
    sample = ['CCCCC' + pattern + 'CCCC', 'CCCCC' + reverse + 'CCCC']
    _, start, forward = scan_windows(sequence_windows(sample, len(pwm)), pwm, np.full(4, .25))
    assert start.tolist() == [5, 5] and forward.tolist() == [True, False]
    primary = family[(family.purpose == 'primary') & (family.view == 'TOP5')]
    assert set(primary.method) == {'MEME', 'STREME', 'ONEHOT', 'PROKBERT', 'CNN_SHIFT'}
    assert primary.candidates.eq(5).all() and primary.positive_groups.eq(66).all() and primary.natural_groups.eq(77).all()
    alignment = pd.read_csv(OUT / 'alignment.tsv', sep='\t', dtype={'seed': str})
    assert len(alignment) == len(primary)
    assert alignment.evaluation_split.eq('development').all() and alignment.candidate_budget.eq(5).all()
    for column in alignment.columns:
        if column.endswith('_sha256'):
            filename = column.removesuffix('_sha256') + '.fasta'
            assert alignment[column].eq(sha(DATA / filename)).all(), column
    selected = pd.read_csv(OUT / 'figure_candidates.tsv', sep='\t', dtype={'seed': str})
    for row in selected.itertuples():
        expected = ranking[(ranking.method == row.method) & (ranking.seed == row.seed) & (ranking['rank'] == row.rank)]
        assert len(expected) == 1 and expected.iloc[0].candidate == row.candidate
    coverage = pd.read_csv(OUT / 'coverage.tsv', sep='\t')
    complete = coverage[(coverage.species == 'tjupan_bacillus_subtilis') & coverage.method.isin(primary.method)]
    assert complete.aligned_to_frozen_data.all()
    save(OUT / 'verification.json', {'status': 'PASSED', 'frozen_outputs_checked': 210,
                                     'historical_files_unchanged': len(read(OUT / 'protected_history.json')),
                                     'candidate_pools': len(pools()), 'family_views_recomputed': recomputed,
                                     'legacy_family_views_identical': len(matched), 'scalar_scan_checks': scalar_checks,
                                     'forward_reverse_coordinate_fixture': 'passed', 'figure_candidates_checked': len(selected),
                                     'holdout_scored': False, 'primary_method_count': 5})
    print('VERIFIED', recomputed, 'family views;', scalar_checks, 'scalar scans;', len(selected), 'figure candidates', flush=True)


if __name__ == '__main__':
    main()
