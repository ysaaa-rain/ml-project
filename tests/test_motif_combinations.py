import pytest

from experiments.analyze_motif_combinations import (
    _nearest_site_pair,
    _parse_fimo_tsv,
    summarize_pair_presence,
)


def test_pair_summary_keeps_marginals_and_paired_cooccurrence_separate():
    result = summarize_pair_presence(
        positive_a={0, 1},
        positive_b={0, 2},
        null_a={0, 1},
        null_b={1, 2},
        sample_count=4,
    )
    assert result["positive_a_hits"] == 2
    assert result["positive_b_hits"] == 2
    assert result["positive_pair_cooccurrence"] == 1
    assert result["null_pair_cooccurrence"] == 1
    assert result["positive_only_pair_count"] == 1
    assert result["null_only_pair_count"] == 1
    assert result["paired_mcnemar_p_value"] == 1.0
    assert 0 <= result["positive_pair_fisher_independence_p_value"] <= 1


def test_pair_summary_rejects_invalid_paired_indexes():
    with pytest.raises(ValueError, match="nonnegative"):
        summarize_pair_presence({-1}, set(), set(), set(), 2)
    with pytest.raises(ValueError, match="exceeds"):
        summarize_pair_presence({2}, set(), set(), set(), 2)


def test_nearest_site_pair_reports_intervening_bases_and_overlap():
    assert _nearest_site_pair(
        [{"start": 1, "stop": 3}], [{"start": 6, "stop": 8}]
    ) == ("A_before_B", 2)
    assert _nearest_site_pair(
        [{"start": 1, "stop": 5}], [{"start": 4, "stop": 8}]
    ) == ("overlapping", -2)
    assert _nearest_site_pair([], [{"start": 1, "stop": 2}]) is None


def test_fimo_text_parser_reads_sites_and_checks_coordinates(tmp_path):
    output = tmp_path / "fimo.tsv"
    output.write_text(
        "# FIMO output\n"
        "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n"
        "species_MEME-1\t\tpositive_000001|label=1\t4\t9\t+\t8.2\t0.01\t0.2\tACGTAC\n",
        encoding="utf-8",
    )
    hit = _parse_fimo_tsv(output)[0]
    assert hit == {
        "motif_id": "species_MEME-1",
        "sequence_name": "positive_000001|label=1",
        "start": 4,
        "stop": 9,
        "strand": "+",
        "p_value": 0.01,
    }

    output.write_text(
        "motif_id\tsequence_name\tstart\tstop\tstrand\tp-value\n"
        "m1\ts1\t0\t3\t+\t0.01\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="invalid FIMO site"):
        _parse_fimo_tsv(output)
