"""Boundary and evaluation safeguards for the scientific comparison."""
import numpy as np

from .pilot import (PROTOCOL, core_tokens, create_data, crop_spec, evaluate,
                    make_pwm, read, reverse_complement, scan)


def test_every_crop_preserves_core_and_token_coordinates():
    sequence = ''.join(np.random.default_rng(7).choice(list('ACGT'), 81))
    for start in range(70):
        left, offset = crop_spec(start, 12, 24)
        crop = sequence[left:left + 24]
        assert len(crop) == 24
        assert crop[offset:offset + 12] == sequence[start:start + 12]
        assert [crop[i:i + 6] for i in range(offset, offset + 7)] == [sequence[i:i + 6] for i in range(start, start + 7)]


def test_mean_and_order_use_identical_token_support():
    tokens = np.arange(2 * 76 * 3).reshape(2, 76, 3)
    cores = core_tokens(tokens, 12, 6)
    assert cores.shape == (2, 70, 7, 3)
    ordered = cores.reshape(2, 70, -1)
    assert np.array_equal(ordered.reshape(2, 70, 7, 3).mean(axis=2), cores.mean(axis=2))
    assert np.array_equal(cores[:, -1], tokens[:, 69:76])


def test_zero_reference_gate_does_not_hide_raw_false_positives():
    cfg = read(PROTOCOL)
    cfg['synthetic_pattern'] = 'CCCCCCCCCCCC'
    _, pwm = make_pwm(['AAAAAAAAAAAA'] * 32)
    data = {'scenario': {'mutation': 0.0}, 'discovery_control': ['ACGT' * 20 + 'A'],
        'calibration_control': ['C' * 81] * 8, 'test': ['A' * 81] * 2,
        'test_control': ['A' * 81] * 2,
        'test_truth': [{'start0': 0}, {'start0': 0}]}
    rows, candidates = evaluate({'bad': pwm}, data, cfg, 'test')
    assert all(row['f1'] == 0 and not row['reference_eligible'] for row in rows)
    assert candidates[0]['raw_blank_fpr'] == 1.0
    assert rows[0]['raw_blank_fpr'] == 1.0


def test_reverse_pwm_scanning_reports_orientation_and_start():
    sequence = 'TTGACATATAAT'
    _, pwm = make_pwm([sequence] * 256)
    scores, starts, forward = scan(['C' * 20 + sequence + 'C' * 49,
        'C' * 20 + reverse_complement(sequence) + 'C' * 49], pwm, np.full(4, .25))
    assert starts.tolist() == [20, 20]
    assert forward.tolist() == [True, False]
    assert np.isclose(scores[0], scores[1])


def test_empty_candidates_are_explicit_zero_recovery():
    cfg = read(PROTOCOL)
    data = {'scenario': {'mutation': .15}, 'discovery_control': ['ACGT' * 20 + 'A'], 'test': ['A' * 81]}
    rows, nulls = evaluate({}, data, cfg, 'test')
    assert len(rows) == 4 and not nulls
    assert all(row['f1'] == 0 and row['candidate_count'] == 0 for row in rows)


def test_split_streams_are_independent_of_other_split_sizes():
    cfg = read(PROTOCOL); smaller = {**cfg, 'discovery_n': 16}
    original = create_data(cfg, cfg['scenarios'][0], cfg['data_seeds'][0])
    changed = create_data(smaller, cfg['scenarios'][0], cfg['data_seeds'][0])
    assert original['test'] == changed['test']
    assert original['development'] == changed['development']
    for split in ['discovery', 'development', 'test']:
        for truth in original[split + '_truth']:
            assert original[split][truth['sequence_index']][truth['start0']:truth['end0']] == truth['fragment']
