import pandas as pd

from experiments.eda import run_eda


def test_run_eda_writes_overall_and_sigma_position_statistics(tmp_path):
    input_path = tmp_path / "clean.tsv"
    output_dir = tmp_path / "eda"
    frame = pd.DataFrame(
        {
            "sequence_id": ["a", "b", "c"],
            "sequence": ["AAA", "ACA", "NCA"],
            "sequence_length": [3, 3, 3],
            "gc_fraction": [0.0, 1 / 3, 1 / 3],
            "species": ["E. coli"] * 3,
            "source_dataset": ["regulondb"] * 3,
            "sigma_factor_type": ["Sigma70", "Sigma70", None],
            "split": ["discovery"] * 3,
            "tss_offset_in_window": [1, 1, 1],
        }
    )
    frame.to_csv(input_path, sep="\t", index=False)

    summary = run_eda(input_path, output_dir)

    overall = pd.read_csv(output_dir / "position_composition.tsv", sep="\t")
    per_sigma = pd.read_csv(output_dir / "position_composition_by_sigma.tsv", sep="\t")
    missingness = pd.read_csv(output_dir / "missingness.tsv", sep="\t")
    assert summary["rows"] == 3
    assert summary["position_count"] == 3
    assert overall.loc[0, "relative_position_bp"] == -1
    assert overall.loc[0, "shannon_entropy_bits"] == 0
    assert overall.loc[0, "information_bits_vs_uniform"] == 2
    assert set(per_sigma["sigma_factor_type"]) == {"Sigma70", "unknown"}
    assert (output_dir / "position_information_overall.png").is_file()
    assert (output_dir / "position_information_by_sigma.png").is_file()
    assert missingness.loc[missingness["column"].eq("sigma_factor_type"), "missing_count"].item() == 1
