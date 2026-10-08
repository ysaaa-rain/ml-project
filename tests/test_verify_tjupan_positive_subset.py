"""Tests for exact, ID-preserving TJU positive_samples redundancy audit."""
import csv
from pathlib import Path

from preprocessing.verify_tjupan_positive_subset import verify_pair


def _write(path: Path, rows: list[dict[str, str]]) -> None:
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=("seq_id", "seq_type", "seq", "label"))
        writer.writeheader()
        writer.writerows(rows)


def _paths(tmp_path):
    d = tmp_path / "Dataset.csv"
    p = tmp_path / "positive_samples.csv"
    first = {"seq_id": "1", "seq_type": "gene", "seq": "ATGC", "label": "1"}
    negative = {"seq_id": "2", "seq_type": "gene", "seq": "CGTA", "label": "0"}
    _write(d, [first, negative])
    _write(p, [first])
    return d, p, first


def test_exact_positive_subset_passes(tmp_path):
    d, p, _ = _paths(tmp_path)
    result = verify_pair(d, p)
    assert result["identical_positive_rows_by_id"]
    assert not result["errors"]


def test_positive_sequence_mismatch_is_rejected(tmp_path):
    d, p, row = _paths(tmp_path)
    _write(p, [{**row, "seq": "ATGT"}])
    assert not verify_pair(d, p)["identical_positive_rows_by_id"]


def test_missing_positive_is_rejected(tmp_path):
    d, p, _ = _paths(tmp_path)
    _write(p, [])
    assert not verify_pair(d, p)["identical_positive_rows_by_id"]


def test_extra_positive_is_rejected(tmp_path):
    d, p, row = _paths(tmp_path)
    _write(p, [row, {"seq_id": "100", "seq_type": "gene", "seq": "AAAA", "label": "1"}])
    assert not verify_pair(d, p)["identical_positive_rows_by_id"]


def test_duplicate_positive_id_is_rejected(tmp_path):
    d, p, row = _paths(tmp_path)
    _write(p, [row, row])
    result = verify_pair(d, p)
    assert not result["identical_positive_rows_by_id"]
    assert "duplicate" in " ".join(result["errors"])


def test_wrong_label_is_rejected(tmp_path):
    d, p, row = _paths(tmp_path)
    _write(p, [{**row, "label": "0"}])
    assert not verify_pair(d, p)["identical_positive_rows_by_id"]


def test_duplicate_dataset_id_is_rejected(tmp_path):
    d, p, row = _paths(tmp_path)
    _write(d, [row, row])
    assert not verify_pair(d, p)["identical_positive_rows_by_id"]
