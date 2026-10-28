"""Check comparison units and scanner coordinates, independent of experiment data."""
import numpy as np
import pandas as pd

from experiments.aligned_motif_display_20261028.pipeline import evaluate, scan_windows, sequence_windows
from experiments.pr01_02_embedding_v2.discovery import make_pwm


def test_same_scanner_preserves_reverse_coordinates():
    pattern = 'AACGTTTA'
    reverse = pattern.translate(str.maketrans('ACGT', 'TGCA'))[::-1]
    _, pwm = make_pwm([pattern] * 20)
    scores, starts, forward = scan_windows(
        sequence_windows(['CCCCC' + pattern + 'CCCC', 'CCCCC' + reverse + 'CCCC'], 8), pwm, np.full(4, .25))
    assert starts.tolist() == [5, 5]
    assert forward.tolist() == [True, False]
    assert np.isclose(scores[0], scores[1])


def test_group_denominators_and_paired_shuffle_are_distinct():
    _, pwm = make_pwm(['AACGTTTA'] * 20)
    hit, blank = 'CCCCAACGTTTACCCCC', 'CCCCCCCCCCCCCCCCC'
    pos = pd.DataFrame({'sequence_id': ['p1', 'p2', 'q', 'mixedp'],
                        'sequence': [hit, hit, blank, hit], 'leakage_group': ['P', 'P', 'Q', 'M']})
    neg = pd.DataFrame({'sequence_id': ['n1', 'n2', 'mixedn'],
                        'sequence': [blank, blank, hit], 'leakage_group': ['N1', 'N2', 'M']})
    ranking = [{'candidate': 'fixture', 'rank': 1, 'threshold': 0., 'width': 8,
                'discovery_score': 0., 'information': 0.}]
    _, family, _ = evaluate({'fixture': pwm}, ranking, pos, neg, [blank] * 4, {'M'}, np.full(4, .25))
    for row in family:
        assert row['positive_groups'] == 2 and row['natural_groups'] == 2
        assert row['positive_rate'] == .5 and row['natural_rate'] == 0
        # Pairing includes all three positive groups; natural comparison excludes M.
        assert row['paired_groups'] == 3 and row['paired_positive_only'] == 2
        assert row['shuffle_rate'] == 0 and row['paired_positive_rate'] == 2 / 3


def test_positions_are_counted_once_per_group():
    _, pwm = make_pwm(['AACGTTTA'] * 20)
    hit, blank = 'CCCCAACGTTTACCCCC', 'CCCCCCCCCCCCCCCCC'
    pos = pd.DataFrame({'sequence_id': ['p1', 'p2'], 'sequence': [hit, hit], 'leakage_group': ['P', 'P']})
    neg = pd.DataFrame({'sequence_id': ['n'], 'sequence': [blank], 'leakage_group': ['N']})
    ranking = [{'candidate': 'fixture', 'rank': 1, 'threshold': 0., 'width': 8,
                'discovery_score': 0., 'information': 0.}]
    _, _, positions = evaluate({'fixture': pwm}, ranking, pos, neg, [blank] * 2, set(), np.full(4, .25), positions=True)
    assert sum(row['groups'] for row in positions if row['class'] == 'positive') == 1
