from experiments.analyze_reference_pwm_sensitivity import (
    _read_raw_hits,
    summarize_ecoli_expected_windows,
    summarize_hits,
    summarize_position_bins,
)


def test_raw_fimo_parser_reads_pvalue_sites_without_requiring_qvalues(tmp_path):
    path = tmp_path / "fimo.tsv"
    path.write_text(
        "# FIMO text output\n"
        "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n"
        "Ecoli_sigma70_minus10\t\tp1|label=1\t20\t25\t+\t8.1\t0.01\t\t\n",
        encoding="utf-8",
    )
    hits = _read_raw_hits(path)
    assert len(hits) == 1
    assert hits[0]["p_value"] == 0.01
    assert hits[0]["start"] == 20


def test_sensitivity_summary_applies_cutoffs_deduplication_and_position_summary():
    sequences = {
        "p1|label=1": {"label": "1", "sequence_sha256": "pos-a"},
        "p2|label=1": {"label": "1", "sequence_sha256": "pos-a"},
        "c1|label=0": {"label": "0", "sequence_sha256": "ctl-a"},
        "shared-p|label=1": {"label": "1", "sequence_sha256": "shared"},
        "shared-c|label=0": {"label": "0", "sequence_sha256": "shared"},
    }
    hits = [
        {"motif_id": "Ecoli_sigma70_minus10", "sequence_name": "p1|label=1", "start": 20, "stop": 25, "strand": "+", "score": 8, "p_value": 0.02},
        {"motif_id": "Ecoli_sigma70_minus10", "sequence_name": "p2|label=1", "start": 21, "stop": 26, "strand": "-", "score": 8, "p_value": 0.03},
        {"motif_id": "Ecoli_sigma70_minus10", "sequence_name": "c1|label=0", "start": 30, "stop": 35, "strand": "-", "score": 8, "p_value": 0.04},
        {"motif_id": "Ecoli_sigma70_minus10", "sequence_name": "shared-p|label=1", "start": 45, "stop": 50, "strand": "+", "score": 8, "p_value": 0.001},
        {"motif_id": "Ecoli_sigma70_minus10", "sequence_name": "shared-c|label=0", "start": 45, "stop": 50, "strand": "+", "score": 8, "p_value": 0.001},
    ]
    rows = summarize_hits("escherichia_coli", sequences, hits, (0.05, 0.01))
    assert len(rows) == 8
    loose = next(row for row in rows if row["reference_motif"] == "Ecoli_sigma70_minus10" and row["site_p_cutoff"] == 0.05)
    strict = next(row for row in rows if row["reference_motif"] == "Ecoli_sigma70_minus10" and row["site_p_cutoff"] == 0.01)
    assert loose["positive_unique_hit_sequences"] == 2
    assert loose["positive_unique_sequences"] == 2
    assert loose["control_unique_hit_sequences"] == 2
    assert loose["cross_label_shared_sequences_excluded_from_fisher"] == 1
    assert loose["positive_best_hit_median_midpoint_1based"] == 35.0
    assert strict["positive_unique_hit_sequences"] == 1
    assert strict["control_unique_hit_sequences"] == 1
    assert all(row["fisher_bh_q_24_tests"] is not None for row in rows)


def test_position_distribution_counts_sequences_at_each_site_and_expected_window():
    sequences = {
        "p|label=1": {"label": "1", "sequence_sha256": "p", "gc_count": 40},
        "c|label=0": {"label": "0", "sequence_sha256": "c", "gc_count": 40},
    }
    hits = [
        {"motif_id": "Ecoli_sigma70_minus10", "sequence_name": "p|label=1", "start": 50, "stop": 55, "strand": "+", "p_value": 0.001},
        {"motif_id": "Ecoli_sigma70_minus10", "sequence_name": "p|label=1", "start": 50, "stop": 55, "strand": "+", "p_value": 0.002},
        {"motif_id": "Ecoli_sigma70_minus10", "sequence_name": "c|label=0", "start": 50, "stop": 55, "strand": "-", "p_value": 0.001},
    ]
    position = summarize_position_bins("escherichia_coli", sequences, hits, (0.001,))
    assert {row["unique_sequences_with_site"] for row in position} == {1}
    windows = summarize_ecoli_expected_windows(sequences, hits, (0.001,), permutations=200)
    plus = next(row for row in windows if row["reference_motif"].endswith("minus10") and row["strand"] == "+")
    minus = next(row for row in windows if row["reference_motif"].endswith("minus10") and row["strand"] == "-")
    assert plus["positive_unique_sequences_hit"] == 1
    assert plus["control_unique_sequences_hit"] == 0
    assert minus["control_unique_sequences_hit"] == 1
    assert all(row["gc_stratified_bh_q_12_tests"] is not None for row in windows)
