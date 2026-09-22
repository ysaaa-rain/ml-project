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
    assert (output_dir / "position_distribution.png").exists()
