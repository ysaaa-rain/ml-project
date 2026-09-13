"""Tests for TSS-anchored windowing and strand normalisation."""

from pathlib import Path

import pandas as pd
import pytest

from preprocessing.pipeline import normalize_metadata
from preprocessing.windowing import (
    apply_tss_window,
    parse_strand,
    parse_tss_position,
    reverse_complement,
)

FIXTURE = Path(__file__).parent / "fixtures" / "tss_promoters.tsv"


def _load() -> pd.DataFrame:
    return pd.read_csv(FIXTURE, sep="\t")


def test_reverse_complement_handles_n():
    assert reverse_complement("ACGTN") == "NACGT"


def test_parse_strand_accepts_common_spellings():
    assert parse_strand("+") == 1
    assert parse_strand("1") == 1
    assert parse_strand("forward") == 1
    assert parse_strand("plus") == 1
    assert parse_strand("-") == -1
    assert parse_strand("-1") == -1
    assert parse_strand("reverse") == -1
    assert parse_strand("crick") == -1
    assert parse_strand("") is None
    assert parse_strand(None) is None
    assert parse_strand(float("nan")) is None
    assert parse_strand("sideways") is None


def test_parse_tss_position_accepts_one_based_and_float_like_values():
    assert parse_tss_position(60) == 60
    assert parse_tss_position("60") == 60
    assert parse_tss_position("60.0") == 60
    assert parse_tss_position(0) is None
    assert parse_tss_position(-5) is None
    assert parse_tss_position("") is None
    assert parse_tss_position(None) is None
    assert parse_tss_position("abc") is None


def test_forward_strand_window_keeps_orientation_and_records_coordinates():
    frame = _load()
    windowed, report = apply_tss_window(frame, upstream=60, downstream=20)

    row = windowed.loc[windowed["sequence_id"].eq("tss_f01")].iloc[0]
    assert row["window_start_1based"] == 1
    assert row["window_end_1based"] == 80
    assert row["tss_offset_in_window"] == 60
    assert row["window_strand"] == 1
    assert len(row["sequence"]) == 80
    # The window is [tss-60, tss+20] = positions 1..80, so the window starts 60
    # bases upstream of the TSS. The fixture puts the -35 box at 1-based 28 and
    # the -10 box at 1-based 48, i.e. 33 and 13 bases upstream of the TSS.
    assert row["sequence"][60 - 33 : 60 - 27] == "TTGACA"
    assert row["sequence"][60 - 13 : 60 - 7] == "TATAAT"


def test_minus_strand_window_is_reverse_complemented_so_elements_read_forward():
    frame = _load()
    windowed, _ = apply_tss_window(frame, upstream=60, downstream=20)

    for identifier in ("tss_r01", "tss_r02"):
        row = windowed.loc[windowed["sequence_id"].eq(identifier)].iloc[0]
        assert row["window_strand"] == -1
        # Window [-60,+20] around TSS 80 -> positions 20..100, 81 bp.
        assert row["window_start_1based"] == 20
        assert row["window_end_1based"] == 100
        assert len(row["sequence"]) == 81
        # After reverse-complementing, the core elements must appear in the same
        # upstream-to-downstream order and at the same offsets as the forward
        # records, which is the whole point of strand normalisation.
        assert row["sequence"][27:33] == "TTGACA", row["sequence"]
        assert row["sequence"][47:53] == "TATAAT", row["sequence"]


def test_forward_and_minus_strand_windows_share_the_same_geometry():
    frame = _load()
    windowed, _ = apply_tss_window(frame, upstream=60, downstream=20)
    forward_window = windowed.loc[windowed["sequence_id"].eq("tss_f01"), "sequence"].iloc[0]
    minus_window = windowed.loc[windowed["sequence_id"].eq("tss_r01"), "sequence"].iloc[0]
    # Same relative offsets even though the two records have different TSS
    # coordinates and opposite strands.
    assert forward_window[27:33] == minus_window[27:33] == "TTGACA"
    assert forward_window[47:53] == minus_window[47:53] == "TATAAT"


def test_strand_token_spellings_agree():
    frame = _load()
    windowed, _ = apply_tss_window(frame, upstream=60, downstream=20)
    minus = windowed.loc[windowed["sequence_id"].eq("tss_r01"), "sequence"].iloc[0]
    minus_word = windowed.loc[windowed["sequence_id"].eq("tss_r02"), "sequence"].iloc[0]
    assert minus == minus_word


def test_records_with_unusable_annotations_are_dropped_and_reported():
    frame = _load()
    windowed, report = apply_tss_window(frame, upstream=60, downstream=20)

    assert report["input_rows"] == 8
    # tss_bad_pos (out of range), tss_no_tss, tss_no_strand are unusable.
    assert report["dropped_rows"] == 3
    assert report["dropped_by_reason"]["window_out_of_range"] == 1
    assert report["dropped_by_reason"]["missing_tss"] == 1
    assert report["dropped_by_reason"]["missing_strand"] == 1
    assert "tss_bad_pos" in report["dropped_sequence_ids"]["window_out_of_range"]
    assert "tss_no_tss" in report["dropped_sequence_ids"]["missing_tss"]
    assert "tss_no_strand" in report["dropped_sequence_ids"]["missing_strand"]

    assert set(windowed["sequence_id"]) == {"tss_f01", "tss_f02", "tss_f03", "tss_r01", "tss_r02"}
    # Every kept row must carry complete window coordinates.
    assert windowed["window_start_1based"].notna().all()
    assert windowed["window_end_1based"].notna().all()
    assert windowed["window_strand"].notna().all()


def test_optional_annotations_keep_rows_when_not_required():
    frame = _load()
    windowed, report = apply_tss_window(
        frame, upstream=60, downstream=20, require_annotations=False
    )
    assert report["dropped_rows"] == 0
    assert len(windowed) == 8
    # Rows without usable annotations keep their original sequence untouched.
    missing = windowed.loc[windowed["sequence_id"].eq("tss_no_tss")].iloc[0]
    assert len(missing["sequence"]) == 100
    assert pd.isna(missing["window_start_1based"])


def test_window_out_of_range_is_dropped_not_clipped():
    frame = pd.DataFrame(
        {
            "sequence_id": ["a"],
            "species": ["E"],
            "source_dataset": ["s"],
            "sequence": ["ACGT" * 5],
            "tss_position": [3],
            "strand": ["+"],
        }
    )
    windowed, report = apply_tss_window(frame, upstream=60, downstream=20)
    assert len(windowed) == 0
    assert report["dropped_by_reason"]["window_out_of_range"] == 1


def test_window_size_is_validated():
    frame = _load()
    with pytest.raises(ValueError):
        apply_tss_window(frame, upstream=-1, downstream=20)
    with pytest.raises(ValueError):
        apply_tss_window(frame, upstream=60, downstream=-1)


def test_windowing_requires_the_sequence_column():
    with pytest.raises(ValueError):
        apply_tss_window(pd.DataFrame({"sequence_id": ["a"]}))


def test_windowing_through_normalize_metadata_records_the_report():
    frame = _load()
    clean, _ = normalize_metadata(frame, tss_window=(60, 20), validation_sources=["dbtbs"])

    assert "window_start_1based" in clean.columns
    # Forward records are TSS 60 with window [0,80]; the TSS 62 record yields
    # [2,82]. The minus-strand records use TSS 80 and yield [20,100].
    assert set(clean["sequence_length"]) == {80, 81}
    report = clean.attrs["tss_windowing"]
    assert report is not None
    assert report["window_upstream"] == 60
    assert report["window_downstream"] == 20
    # The windowed table feeds the existing split logic unchanged.
    assert set(clean["split"]) <= {"discovery", "validation"}


def test_normalize_metadata_rejects_windowing_without_annotation_columns():
    frame = pd.DataFrame(
        {
            "sequence_id": ["a"],
            "species": ["E"],
            "source_dataset": ["s"],
            "sequence": ["ACGT" * 20],
        }
    )
    with pytest.raises(ValueError, match="annotation columns"):
        normalize_metadata(frame, tss_window=(60, 20))


def test_windowing_defaults_to_disabled():
    frame = pd.DataFrame(
        {
            "sequence_id": ["a"],
            "species": ["E"],
            "source_dataset": ["s"],
            "sequence": ["ACGT" * 20],
        }
    )
    clean, _ = normalize_metadata(frame)
    assert "window_start_1based" not in clean.columns
    assert clean.attrs.get("tss_windowing") is None
    assert len(clean.iloc[0]["sequence"]) == 80
