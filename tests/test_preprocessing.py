from pathlib import Path

import pandas as pd

from preprocessing.pipeline import normalize_metadata, write_background_fasta
from preprocessing.schema import normalize_sequence, validate_metadata


FIXTURE = Path(__file__).parent / "fixtures" / "sample_promoters.tsv"


def test_sequence_normalization_and_validation():
    assert normalize_sequence(" au gc\n") == "ATGC"
    frame = pd.DataFrame(
        {
            "sequence_id": ["ok", "bad"],
            "species": ["A", "A"],
            "source_dataset": ["source", "source"],
            "sequence": ["ATGC", "ATGX"],
        }
    )
    report = validate_metadata(frame)
    assert report.invalid_row_count == 1
    assert report.invalid_sequence_ids == ["bad"]


def test_normalize_metadata_filters_quality_and_duplicates():
    original = pd.read_csv(FIXTURE, sep="\t")
    clean, report = normalize_metadata(
        original,
        max_n_fraction=0.1,
        validation_sources=["dbtbs"],
    )
    assert report.row_count == 8
    assert len(clean) == 6
    assert clean["sequence_id"].is_unique
    assert set(clean["split"]) == {"discovery", "validation"}
    assert set(clean.loc[clean["split"].eq("validation"), "source_dataset"]) == {"dbtbs"}


def test_shuffled_background_is_deterministic(tmp_path):
    original = pd.read_csv(FIXTURE, sep="\t")
    clean, _ = normalize_metadata(original)
    first = tmp_path / "first.fasta"
    second = tmp_path / "second.fasta"
    write_background_fasta(clean, first, replicates=2, seed=7)
    write_background_fasta(clean, second, replicates=2, seed=7)
    assert first.read_text(encoding="utf-8") == second.read_text(encoding="utf-8")
