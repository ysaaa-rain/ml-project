"""TSS-anchored windowing and strand normalisation.

Motif positions are only comparable across records when every sequence uses the
same coordinate frame. This module puts each promoter into a window relative to
its transcription start site (TSS):

- the window is ``[-upstream, +downstream]`` around the TSS;
- coordinates are 1-based and inclusive, because the source databases annotate
  TSS positions that way;
- on the minus strand the extracted region is reverse-complemented, so
  "upstream" always means "towards the promoter" and a relative coordinate has
  the same meaning on both strands.

Records whose window would extend past the available sequence are reported and
dropped instead of being silently truncated or shifted.
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from .schema import normalize_sequence

DEFAULT_WINDOW_UPSTREAM = 60
DEFAULT_WINDOW_DOWNSTREAM = 20

FORWARD_STRAND_TOKENS = frozenset({"+", "1", "+1", "f", "forward", "plus", "watson"})
REVERSE_STRAND_TOKENS = frozenset({"-", "-1", "r", "reverse", "minus", "crick"})

WINDOW_COLUMNS = (
    "window_upstream",
    "window_downstream",
    "window_start_1based",
    "window_end_1based",
    "window_strand",
    "tss_offset_in_window",
)


def reverse_complement(sequence: str) -> str:
    """Return the reverse complement of a DNA sequence."""

    return sequence.upper().translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def parse_strand(value: Any) -> int | None:
    """Map a strand annotation to ``+1``/``-1``, or ``None`` if unrecognised."""

    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    token = str(value).strip().lower()
    if not token:
        return None
    if token in FORWARD_STRAND_TOKENS:
        return 1
    if token in REVERSE_STRAND_TOKENS:
        return -1
    return None


def parse_tss_position(value: Any) -> int | None:
    """Parse a 1-based TSS coordinate, returning ``None`` when unusable."""

    if value is None or (isinstance(value, float) and pd.isna(value)):
        return None
    try:
        position = int(float(str(value).strip()))
    except (TypeError, ValueError):
        return None
    return position if position > 0 else None


def apply_tss_window(
    frame: pd.DataFrame,
    *,
    upstream: int = DEFAULT_WINDOW_UPSTREAM,
    downstream: int = DEFAULT_WINDOW_DOWNSTREAM,
    sequence_column: str = "sequence",
    tss_column: str = "tss_position",
    strand_column: str = "strand",
    require_annotations: bool = True,
) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Extract a TSS-anchored, strand-normalised window for every row.

    Returns the windowed frame (with the ``window_*`` columns added) and a
    report describing what was dropped and why. When ``require_annotations`` is
    false, rows without usable TSS/strand annotations are passed through
    unchanged instead of being dropped.
    """

    if upstream < 0 or downstream < 0:
        raise ValueError("window sizes must not be negative")
    if sequence_column not in frame:
        raise ValueError(f"frame is missing the sequence column: {sequence_column}")

    work = frame.copy()
    sequences = work[sequence_column].map(normalize_sequence)

    windowed_sequences: list[str] = []
    starts: list[int | None] = []
    ends: list[int | None] = []
    offsets: list[int | None] = []
    strands: list[int | None] = []
    keep: list[bool] = []

    dropped: dict[str, int] = {
        "missing_tss": 0,
        "missing_strand": 0,
        "window_out_of_range": 0,
        "empty_window": 0,
    }
    dropped_ids: dict[str, list[str]] = {reason: [] for reason in dropped}

    def record(
        sequence: str,
        *,
        keep_row: bool,
        reason: str | None = None,
        window: str | None = None,
        start: int | None = None,
        end: int | None = None,
        offset: int | None = None,
        strand: int | None = None,
        identifier: str = "",
    ) -> None:
        """Append one row's outcome, counting the drop reason when there is one."""

        if reason is not None:
            dropped[reason] += 1
            dropped_ids[reason].append(identifier)
        keep.append(keep_row)
        windowed_sequences.append(sequence if window is None else window)
        starts.append(start)
        ends.append(end)
        offsets.append(offset)
        strands.append(strand)

    for index, sequence in sequences.items():
        identifier = str(work.at[index, "sequence_id"]) if "sequence_id" in work else str(index)
        tss = parse_tss_position(work.at[index, tss_column]) if tss_column in work else None
        strand = parse_strand(work.at[index, strand_column]) if strand_column in work else None

        if tss is None:
            record(sequence, keep_row=not require_annotations, reason="missing_tss", identifier=identifier)
            continue
        if strand is None:
            record(sequence, keep_row=not require_annotations, reason="missing_strand", identifier=identifier)
            continue

        # Window around the TSS: [tss - upstream, tss + downstream], 1-based and
        # inclusive. ``raw_start`` may be 0 when the TSS sits exactly
        # ``upstream`` bases from the start of the sequence, which is a valid
        # although tight window rather than an error. Anything that would reach
        # past either end of the sequence is rejected instead of being clipped.
        raw_start = tss - upstream
        raw_end = tss + downstream
        if raw_start < 0 or raw_end > len(sequence):
            record(
                sequence,
                keep_row=not require_annotations,
                reason="window_out_of_range",
                identifier=identifier,
            )
            continue

        window = sequence[max(0, raw_start - 1) : raw_end]
        if not window:
            record(sequence, keep_row=False, reason="empty_window", identifier=identifier)
            continue

        if strand == -1:
            window = reverse_complement(window)
        # Position 0 of the window is the base at -upstream relative to the TSS.
        record(
            sequence,
            keep_row=True,
            window=window,
            start=max(1, raw_start),
            end=raw_end,
            offset=upstream,
            strand=strand,
        )

    work[sequence_column] = windowed_sequences
    work["window_start_1based"] = starts
    work["window_end_1based"] = ends
    work["tss_offset_in_window"] = offsets
    work["window_strand"] = strands
    work["window_upstream"] = upstream
    work["window_downstream"] = downstream
    work = work.loc[keep].copy()

    report = {
        "window_upstream": upstream,
        "window_downstream": downstream,
        "input_rows": int(len(frame)),
        "output_rows": int(len(work)),
        "dropped_rows": int(len(frame) - len(work)),
        "dropped_by_reason": dropped,
        "dropped_sequence_ids": dropped_ids,
        "require_annotations": require_annotations,
    }
    return work, report
