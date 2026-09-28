import pandas as pd
import pytest

from experiments.prepare_motif_inputs import run
from preprocessing.shuffle import dinucleotide_counts


def _records():
    return pd.DataFrame({
        "sequence_id": ["a", "b", "c"],
        "sequence": ["ACGTA", "TAGCT", "GATCA"],
        "source_dataset": ["regulondb"] * 3,
        "split": ["discovery"] * 3,
        "sigma_factor_type": ["Sigma70", "Sigma70", "unknown"],
        "tss_position": [3] * 3,
        "strand": ["+", "-", "+"],
        "window_upstream": [2] * 3,
        "window_downstream": [2] * 3,
        "tss_offset_in_window": [2] * 3,
        "window_strand": [1, -1, 1],
    })


def _fasta(path):
    lines = path.read_text(encoding="utf-8").splitlines()
    return dict(zip((line[1:] for line in lines[::2]), lines[1::2]))


def test_prepare_motif_inputs_keeps_groups_and_matched_n1_controls(tmp_path):
    source = tmp_path / "clean.tsv"
    _records().to_csv(source, sep="\t", index=False)
    manifest = run(source, tmp_path / "motifs", upstream=2, downstream=2, min_group_size=2, seed=7)
    assert set(manifest["outputs"]) == {"all", "Sigma70"}
    assert manifest["outputs"]["all"]["sequence_count"] == 3
    assert manifest["unknown_sigma_count"] == 1
    primary = _fasta(tmp_path / "motifs" / "all.fasta")
    control = _fasta(tmp_path / "motifs" / "all_n1.fasta")
    assert set(control) == {f"{identifier}__dinucleotide1" for identifier in primary}
    for identifier, sequence in primary.items():
        assert dinucleotide_counts(sequence) == dinucleotide_counts(control[f"{identifier}__dinucleotide1"])
    assert set(_fasta(tmp_path / "motifs" / "Sigma70.fasta")) == {"a", "b"}


def test_prepare_motif_inputs_rejects_mixed_sources_before_writing(tmp_path):
    source = tmp_path / "clean.tsv"
    records = _records()
    records.loc[2, "source_dataset"] = "dbtbs"
    records.to_csv(source, sep="\t", index=False)
    output = tmp_path / "motifs"
    with pytest.raises(ValueError, match="specified source"):
        run(source, output, upstream=2, downstream=2)
    assert not output.exists()
