import numpy as np
import pytest

from experiments.analyze_motif_strength import (
    _bh,
    _consume_fimo_stdout,
    _find_overlaps,
    _partial_spearman,
    _read_strength_table,
    _reverse_complement,
)


def test_strength_reader_strips_headers_cells_and_checks_unique_sequences(tmp_path):
    path = tmp_path / "strength.tsv"
    path.write_text(" strength\t promoter\n12.5\t" + "a" * 50 + " \n", encoding="utf-8")
    rows = _read_strength_table(path)
    assert rows == [{"sequence": "A" * 50, "strength": 12.5}]

    path.write_text("strength\tpromoter\n1\t" + "A" * 50 + "\n2\t" + "A" * 50 + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="duplicate sequences"):
        _read_strength_table(path)


def test_overlap_excludes_substrings_and_reverse_complements_from_both_labels():
    positive = "ACGT" * 20 + "A"
    control = "T" + "TGCA" * 20
    pos_fragment = positive[4:54]
    control_fragment = control[8:58]
    rows = [
        {"sequence": pos_fragment, "strength": 1},
        {"sequence": _reverse_complement(control_fragment), "strength": 2},
        {"sequence": "G" * 50, "strength": 3},
    ]
    excluded, counts = _find_overlaps(
        rows,
        [("p|label=1", positive, "positive"), ("n|label=0", control, "control")],
    )
    assert len(excluded) == 2
    assert counts["strength_sequences_excluded_for_overlap"] == 2
    assert counts["strength_sequences_matching_positive_any_50nt_window_or_reverse_complement"] == 1
    assert counts["strength_sequences_matching_control_any_50nt_window_or_reverse_complement"] == 1


def test_bh_and_partial_spearman_are_bounded_and_monotonic():
    assert _bh([0.01, 0.04, 0.03]) == pytest.approx([0.03, 0.04, 0.04])
    rho, p_value = _partial_spearman(
        np.asarray([1, 2, 3, 4, 5]),
        np.asarray([1, 2, 3, 4, 5]),
        np.asarray([5, 4, 3, 2, 1]),
    )
    assert -1 <= rho <= 1
    assert 0 <= p_value <= 1


def test_fimo_stream_aggregates_sites_without_preserving_sequence_text():
    lines = [
        "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\n",
        "m1\t\ts1\t1\t6\t+\t4.2\t0.02\n",
        "m1\t\ts1\t8\t13\t-\t7.1\t0.0005\n",
    ]
    state, site_count = _consume_fimo_stdout(lines, ["m1"], {"s1", "s2"})
    assert site_count == 2
    assert state["m1"]["s1"]["best_score"] == pytest.approx(7.1)
    assert state["m1"]["s1"]["min_p_value"] == pytest.approx(0.0005)
    assert state["m1"]["s1"]["hits_p_le_0.05"] == 2
    assert state["m1"]["s1"]["hits_p_le_0.001"] == 1
    assert state["m1"]["s2"]["min_p_value"] == 1.0
