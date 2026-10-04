from experiments.analyze_regulondb_sigma_fimo import (
    _completion_message,
    _nearest_pair_gap,
    _paired_discordance,
    _pair_spacing_rows,
)


def test_nearest_pair_gap_counts_unmatched_bases_between_adjacent_sites():
    first = {"start": 10, "stop": 15}
    second = {"start": 20, "stop": 24}
    assert _nearest_pair_gap(first, second) == (4, "a_before_b")


def test_nearest_pair_gap_marks_overlapping_site_intervals_negative():
    first = {"start": 10, "stop": 16}
    second = {"start": 15, "stop": 20}
    assert _nearest_pair_gap(first, second) == (-2, "a_before_b")


def test_nearest_pair_gap_tracks_reverse_input_order():
    first = {"start": 24, "stop": 29}
    second = {"start": 10, "stop": 15}
    assert _nearest_pair_gap(first, second) == (8, "b_before_a")


def test_paired_discordance_matches_primary_and_n1_by_original_record_id():
    metadata = {
        "Sigma70|positive|p1": {"pair_id": "p1"},
        "Sigma70|positive|p2": {"pair_id": "p2"},
        "Sigma70|positive|p3": {"pair_id": "p3"},
        "Sigma70|control|p1": {"pair_id": "p1"},
        "Sigma70|control|p2": {"pair_id": "p2"},
        "Sigma70|control|p3": {"pair_id": "p3"},
    }
    assert _paired_discordance(
        {"Sigma70|positive|p1", "Sigma70|positive|p2"},
        {"Sigma70|control|p2", "Sigma70|control|p3"},
        metadata,
    ) == (1, 1)


def test_pair_spacing_compares_positive_and_its_matched_n1_control():
    metadata = {
        "Sigma70|positive|p1": {"group": "Sigma70", "label": "positive", "pair_id": "p1"},
        "Sigma70|positive|p2": {"group": "Sigma70", "label": "positive", "pair_id": "p2"},
        "Sigma70|control|p1": {"group": "Sigma70", "label": "control", "pair_id": "p1"},
        "Sigma70|control|p2": {"group": "Sigma70", "label": "control", "pair_id": "p2"},
    }
    motifs = {
        "MEME-1": {"consensus": "AAAAAA", "ordinal": 1},
        "MEME-2": {"consensus": "CCCCCC", "ordinal": 2},
    }
    hits = {
        "MEME-1": {
            "Sigma70|positive|p1": [{"start": 10, "stop": 15, "p_value": 0.001}],
            "Sigma70|positive|p2": [{"start": 10, "stop": 15, "p_value": 0.001}],
            "Sigma70|control|p2": [{"start": 10, "stop": 15, "p_value": 0.001}],
        },
        "MEME-2": {
            "Sigma70|positive|p1": [{"start": 20, "stop": 25, "p_value": 0.001}],
            "Sigma70|control|p2": [{"start": 20, "stop": 25, "p_value": 0.001}],
        },
    }
    rows = _pair_spacing_rows(
        "Sigma70", motifs, hits,
        {"Sigma70|positive|p1", "Sigma70|positive|p2"}, metadata,
    )
    assert rows[0]["positive_sequences_with_both_motifs"] == 1
    assert rows[0]["N1_sequences_with_both_motifs"] == 1
    assert rows[0]["paired_positive_only"] == 1
    assert rows[0]["paired_N1_only"] == 1
    assert rows[0]["positive_median_nearest_signed_gap_bp"] == 4


def test_completion_message_sums_only_table_row_counts():
    result = {
        "rows": {"coverage.tsv": 360, "spacing.tsv": 810, "significant": {"0.01": 4}},
        "fimo_runs": {"Sigma70": {"site_rows_p_le_0.05": 100}},
    }
    assert _completion_message(result) == "Completed cross-sigma FIMO: 1170 summary rows, 100 site rows across six scans"
