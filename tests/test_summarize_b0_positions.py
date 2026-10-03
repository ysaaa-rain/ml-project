import pandas as pd

from experiments.summarize_b0_positions import run


def test_summarize_b0_positions_writes_table_and_plot(tmp_path):
    input_path = tmp_path / "aligned.tsv"
    hits_path = tmp_path / "hits.tsv"
    output_dir = tmp_path / "positions"
    pd.DataFrame(
        {
            "sequence_id": ["a", "b"],
            "tss_offset_in_window": [60, 60],
            "sigma_factor_type": ["Sigma70", "Sigma24"],
        }
    ).to_csv(input_path, sep="\t", index=False)
    pd.DataFrame(
        {
            "motif_name": ["minus_35_box", "minus_10_box"],
            "sequence_id": ["a", "b"],
            "start": [27, 47],
            "end": [33, 53],
            "strand": ["+", "+"],
            "score": [10.0, 10.0],
        }
    ).to_csv(hits_path, sep="\t", index=False)

    result = run(input_path, hits_path, output_dir)
    assert result["hit_count"] == 2
    assert result["positive_strand_medians"]["minus_35_box"] == -33.0
    assert (output_dir / "position_summary.tsv").exists()
    assert (output_dir / "position_summary_by_sigma.tsv").exists()
    presence = pd.read_csv(output_dir / "presence_by_sigma.tsv", sep="\t")
    assert len(presence) == 4
    assert set(presence["sequence_hit_fraction"]) == {0.0, 1.0}
    assert (output_dir / "presence_overall.tsv").exists()
    by_sigma = pd.read_csv(output_dir / "position_summary_by_sigma.tsv", sep="\t")
    assert set(by_sigma["sigma_factor_type"]) == {"Sigma70", "Sigma24"}
    assert set(by_sigma["group_size"]) == {1}
    assert set(by_sigma["sequence_hit_fraction"]) == {1.0}
    assert (output_dir / "position_distribution.png").exists()
