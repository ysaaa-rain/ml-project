"""Recover DBTBS 4.1 promoter strand/TSS coordinates from saved official pages.

This is a data-repair step only. It does not run motif discovery.

DBTBS 4.1 defines promoter ``Location`` relative to the transcription start
site and reports ``Absolute position`` against accession NC_000964. The operon
page also reports gene ``Direction``. We use those explicit fields; missing or
inconsistent coordinates remain unresolved instead of being guessed.
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path
from typing import Any

import pandas as pd

from .download_dbtbs import (
    DNA_RE,
    POSITION_RE,
    _HTMLTableParser,
    _clean_sequence,
    _clean_text,
    sha256_bytes,
)


LOCATION_RE = re.compile(r"^([+-]?\d+):([+-]?\d+)$")
DBTBS_COORDINATE_ACCESSION = "NC_000964"


def _normalised_header(row: list[str]) -> list[str]:
    return [re.sub(r"[^a-z]", "", cell.lower()) for cell in row]


def _find_promoter_header(rows: list[list[str]]) -> int | None:
    for index, row in enumerate(rows):
        normalised = _normalised_header(row)
        if any("bindingfactor" in cell for cell in normalised) and any(
            "bindingseq" in cell for cell in normalised
        ):
            return index
    return None


def _page_gene_direction(rows: list[list[str]]) -> tuple[str | None, list[str]]:
    """Return a unique operon gene direction when the page supports one.

    DBTBS regulated-operon pages list gene direction in the Genes table. An
    operon should have one transcription direction; if a page exposes mixed or
    missing directions we refuse to infer a strand for its promoter records.
    """

    directions: list[str] = []
    for row in rows:
        if len(row) < 4:
            continue
        direction = _clean_text(row[2])
        genome_position = _clean_text(row[3])
        if direction in {"+", "-"} and POSITION_RE.fullmatch(genome_position):
            directions.append(direction)
    unique = sorted(set(directions))
    return (unique[0] if len(unique) == 1 else None), unique


def _recover_coordinate_fields(
    *,
    location: str,
    absolute_position: str,
    sequence: str,
    strand: str | None,
) -> dict[str, Any]:
    """Recover local TSS and genomic TSS using only explicit DBTBS fields.

    In DBTBS, a location such as ``-43:+15`` spans 58 displayed nucleotides;
    the zero point is the TSS. Therefore the TSS has zero-based offset 43 in
    the displayed sequence (one-based ``tss_position`` 44). This is accepted
    only when the relative span and the inclusive absolute span both match the
    actual displayed sequence length.
    """

    result: dict[str, Any] = {
        "strand": strand if strand in {"+", "-"} else pd.NA,
        "promoter_location_start": pd.NA,
        "promoter_location_end": pd.NA,
        "source_window_start": pd.NA,
        "source_window_end": pd.NA,
        "tss_position": pd.NA,
        "source_tss_coordinate": pd.NA,
        "coordinate_status": "unresolved",
        "coordinate_issue": "",
    }

    if location.upper() == "ND" or not location:
        result["coordinate_issue"] = "location_nd_or_missing"
        return result
    location_match = LOCATION_RE.fullmatch(location)
    if location_match is None:
        result["coordinate_issue"] = "location_unparseable"
        return result
    relative_start = int(location_match.group(1))
    relative_end = int(location_match.group(2))
    result["promoter_location_start"] = relative_start
    result["promoter_location_end"] = relative_end

    position_match = POSITION_RE.fullmatch(absolute_position)
    if position_match is None:
        result["coordinate_issue"] = "absolute_position_nd_or_unparseable"
        return result
    absolute_start = int(position_match.group(1))
    absolute_end = int(position_match.group(2))
    result["source_window_start"] = absolute_start
    result["source_window_end"] = absolute_end

    if strand not in {"+", "-"}:
        result["coordinate_issue"] = "gene_direction_missing_or_mixed"
        return result
    if not (relative_start <= 0 < relative_end):
        result["coordinate_issue"] = "location_does_not_span_tss"
        return result

    relative_span = relative_end - relative_start
    absolute_span = absolute_end - absolute_start + 1
    if relative_span != len(sequence):
        result["coordinate_issue"] = "location_length_mismatch"
        return result
    if absolute_span != len(sequence):
        result["coordinate_issue"] = "absolute_length_mismatch"
        return result

    tss_offset = -relative_start
    if not (0 <= tss_offset < len(sequence)):
        result["coordinate_issue"] = "tss_offset_out_of_range"
        return result

    result["tss_position"] = tss_offset + 1  # sequence-local, 1-based
    result["source_tss_coordinate"] = (
        absolute_start + tss_offset
        if strand == "+"
        else absolute_end - tss_offset
    )
    result["coordinate_status"] = "recovered_from_dbtbs"
    result["coordinate_issue"] = ""
    return result


def extract_recovered_promoter_rows(
    html: str, *, page_name: str
) -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Extract sigma-promoter rows and recover coordinate fields when possible."""

    parser = _HTMLTableParser()
    parser.feed(html)
    rows = parser.rows
    header_index = _find_promoter_header(rows)
    page_strand, page_directions = _page_gene_direction(rows)
    counts = {
        "promoter_rows_seen": 0,
        "sigma_promoters_kept": 0,
        "non_sigma_promoters_skipped": 0,
        "invalid_or_missing_sequence_skipped": 0,
        "coordinates_recovered": 0,
        "coordinates_unresolved": 0,
        "location_nd_or_missing": 0,
        "mixed_or_missing_gene_direction_pages": int(page_strand is None),
    }
    if header_index is None:
        return [], counts

    records: list[dict[str, Any]] = []
    for row_index, row in enumerate(rows[header_index + 1 :], start=1):
        if len(row) < 6:
            continue
        binding_factor = _clean_text(row[0])
        regulation = _clean_text(row[1])
        if regulation.lower() != "promoter":
            continue
        counts["promoter_rows_seen"] += 1
        if not re.match(r"^sig(?:ma)?[a-z0-9]+$", binding_factor, re.IGNORECASE):
            counts["non_sigma_promoters_skipped"] += 1
            continue

        location = _clean_text(row[2])
        absolute_position = _clean_text(row[3])
        sequence = _clean_sequence(row[4])
        evidence = _clean_text(row[5])
        if not sequence or not DNA_RE.fullmatch(sequence):
            counts["invalid_or_missing_sequence_skipped"] += 1
            continue

        coordinates = _recover_coordinate_fields(
            location=location,
            absolute_position=absolute_position,
            sequence=sequence,
            strand=page_strand,
        )
        if coordinates["coordinate_status"] == "recovered_from_dbtbs":
            counts["coordinates_recovered"] += 1
        else:
            counts["coordinates_unresolved"] += 1
            if coordinates["coordinate_issue"] == "location_nd_or_missing":
                counts["location_nd_or_missing"] += 1

        source_record_id = f"dbtbs_v4.1:{page_name}:promoter:{row_index}"
        records.append(
            {
                "sequence_id": source_record_id,
                "species": "Bacillus subtilis",
                "taxon_id": 224308,
                "assembly_accession": DBTBS_COORDINATE_ACCESSION,
                "source_dataset": "dbtbs_v4.1",
                "sequence": sequence,
                "sequence_length": len(sequence),
                "tss_position": coordinates["tss_position"],
                "strand": coordinates["strand"],
                "sigma_factor_type": binding_factor,
                "known_element_annotations": pd.NA,
                "promoter_strength": pd.NA,
                "evidence_level": "experimental",
                "source_record_id": source_record_id,
                "source_page": page_name,
                "promoter_location": location,
                "promoter_location_start": coordinates["promoter_location_start"],
                "promoter_location_end": coordinates["promoter_location_end"],
                "source_window_start": coordinates["source_window_start"],
                "source_window_end": coordinates["source_window_end"],
                "source_tss_coordinate": coordinates["source_tss_coordinate"],
                "coordinate_status": coordinates["coordinate_status"],
                "coordinate_issue": coordinates["coordinate_issue"],
                "coordinate_reference": DBTBS_COORDINATE_ACCESSION,
                "page_gene_directions": ",".join(page_directions),
                "source_evidence": evidence,
                "source_binding_factor": binding_factor,
            }
        )
        counts["sigma_promoters_kept"] += 1
    return records, counts


def recover_saved_pages(
    *,
    pages_dir: str | Path,
    metadata_output: str | Path,
    report_output: str | Path,
) -> dict[str, Any]:
    """Reparse a frozen DBTBS page snapshot without downloading anything."""

    pages_path = Path(pages_dir)
    page_paths = sorted(pages_path.glob("*.html"))
    if not page_paths:
        raise FileNotFoundError(f"no DBTBS HTML pages found under {pages_path}")

    all_records: list[dict[str, Any]] = []
    aggregate: dict[str, int] = {}
    page_hashes: dict[str, str] = {}
    for page_path in page_paths:
        payload = page_path.read_bytes()
        page_hashes[page_path.name] = sha256_bytes(payload)
        records, counts = extract_recovered_promoter_rows(
            payload.decode("utf-8", errors="replace"), page_name=page_path.name
        )
        all_records.extend(records)
        for key, value in counts.items():
            aggregate[key] = aggregate.get(key, 0) + value

    frame = pd.DataFrame(all_records).sort_values("sequence_id").reset_index(drop=True)
    if frame.empty:
        raise RuntimeError("DBTBS coordinate recovery produced no usable sigma-promoter rows")

    output_path = Path(metadata_output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(output_path, sep="\t", index=False)

    status_counts = frame["coordinate_status"].value_counts(dropna=False).to_dict()
    issue_counts = frame.loc[
        frame["coordinate_status"] != "recovered_from_dbtbs", "coordinate_issue"
    ].value_counts(dropna=False).to_dict()
    report = {
        "source_id": "dbtbs_v4.1",
        "operation": "coordinate_recovery_from_frozen_pages",
        "pages_dir": str(pages_path),
        "pages_seen": len(page_paths),
        "page_hashes": page_hashes,
        "coordinate_reference": DBTBS_COORDINATE_ACCESSION,
        "coordinate_reference_note": (
            "DBTBS terms state that Absolute position was recalculated using "
            "NCBI accession NC_000964. Do not silently substitute a newer accession "
            "version without sequence/coordinate validation."
        ),
        "metadata_rows": len(frame),
        "status_counts": status_counts,
        "issue_counts": issue_counts,
        "parse_counts": aggregate,
        "metadata_output": str(output_path),
        "metadata_sha256": sha256_bytes(output_path.read_bytes()),
    }
    report_path = Path(report_output)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pages-dir", default="data/raw/dbtbs_v4.1_20261005/pages"
    )
    parser.add_argument(
        "--metadata-output",
        default="data/interim/dbtbs_metadata_recovered_20261005.tsv",
    )
    parser.add_argument(
        "--report-output",
        default="data/interim/dbtbs_coordinate_recovery_20261005.json",
    )
    args = parser.parse_args()
    report = recover_saved_pages(
        pages_dir=args.pages_dir,
        metadata_output=args.metadata_output,
        report_output=args.report_output,
    )
    print(f"metadata_rows={report['metadata_rows']}")
    print(f"status_counts={report['status_counts']}")
    print(f"issue_counts={report['issue_counts']}")
    print(f"metadata_sha256={report['metadata_sha256']}")


if __name__ == "__main__":
    main()
