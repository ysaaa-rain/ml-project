"""Map source-specific promoter tables into the project metadata contract.

The adapters deliberately only rename fields and add the explicitly supplied
source identifier. They do not infer sigma factors, TSS coordinates, strand,
or promoter sequences. Missing optional annotations remain missing so the
downstream quality report can expose the limitation.
"""

from __future__ import annotations

import json
import re
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Mapping
from urllib.parse import unquote

import pandas as pd

from .schema import normalize_sequence


CANONICAL_COLUMNS = (
    "sequence_id",
    "species",
    "taxon_id",
    "assembly_accession",
    "source_dataset",
    "sequence",
    "sequence_length",
    "tss_position",
    "strand",
    "sigma_factor_type",
    "known_element_annotations",
    "promoter_strength",
    "evidence_level",
    "source_record_id",
)

REQUIRED_SOURCE_FIELDS = ("sequence_id", "species", "sequence")

OPTIONAL_SOURCE_ALIASES: dict[str, tuple[str, ...]] = {
    "taxon_id": ("taxon_id", "taxid", "ncbi_taxon_id"),
    "assembly_accession": (
        "assembly_accession",
        "assembly",
        "genome_accession",
        "accession",
    ),
    "sequence_length": ("sequence_length", "length", "promoter_length"),
    "tss_position": (
        "tss_position",
        "tss",
        "tss_coordinate",
        "tss_coord",
        "transcription_start_site",
    ),
    "strand": ("strand", "orientation", "direction", "gene_strand"),
    "sigma_factor_type": (
        "sigma_factor_type",
        "sigma_factor",
        "sigma",
        "sigma_type",
        "factor_sigma",
    ),
    "known_element_annotations": (
        "known_element_annotations",
        "known_elements",
        "element_annotations",
        "promoter_elements",
    ),
    "promoter_strength": (
        "promoter_strength",
        "strength",
        "expression",
        "activity",
    ),
    "evidence_level": (
        "evidence_level",
        "evidence",
        "evidence_type",
        "support_level",
    ),
    "source_record_id": (
        "source_record_id",
        "record_id",
        "promoter_id",
        "promoter_record_id",
        "regulon_record_id",
        "db_record_id",
    ),
}

REQUIRED_SOURCE_ALIASES: dict[str, tuple[str, ...]] = {
    "sequence_id": (
        "sequence_id",
        "promoter_id",
        "promoter_record_id",
        "record_id",
        "id",
        "name",
    ),
    "species": (
        "species",
        "organism",
        "organism_name",
        "species_name",
        "taxon_name",
    ),
    "sequence": (
        "sequence",
        "promoter_sequence",
        "dna_sequence",
        "sequence_dna",
        "seq",
    ),
}


class SourceMappingError(ValueError):
    """Raised when a source table cannot satisfy the required field contract."""


@dataclass
class AdapterReport:
    """Auditable summary of one source-table adaptation."""

    source_id: str
    input_rows: int
    output_rows: int
    mapped_columns: dict[str, str]
    missing_required_fields: list[str]
    missing_optional_fields: list[str]
    derived_fields: list[str]
    warnings: list[str]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def _normalise_column_name(value: Any) -> str:
    """Make column matching tolerant of case, spaces, hyphens, and punctuation."""

    text = str(value).strip().lower()
    return re.sub(r"_+", "_", re.sub(r"[^a-z0-9]+", "_", text)).strip("_")


def _column_lookup(frame: pd.DataFrame) -> dict[str, list[str]]:
    lookup: dict[str, list[str]] = {}
    for column in frame.columns:
        lookup.setdefault(_normalise_column_name(column), []).append(str(column))
    return lookup


def _values_conflict(frame: pd.DataFrame, first: str, second: str) -> bool:
    """Return whether two candidate aliases carry different non-empty values."""

    left = frame[first].astype("string").fillna("").str.strip()
    right = frame[second].astype("string").fillna("").str.strip()
    comparable = left.ne("") & right.ne("")
    return bool((left[comparable] != right[comparable]).any())


def _resolve_column(
    frame: pd.DataFrame,
    lookup: dict[str, list[str]],
    canonical: str,
    aliases: tuple[str, ...],
    explicit_map: Mapping[str, str],
    warnings: list[str],
) -> str | None:
    if canonical in explicit_map:
        raw_name = explicit_map[canonical]
        if raw_name not in frame.columns:
            raise SourceMappingError(
                f"explicit mapping for {canonical!r} points to missing column {raw_name!r}"
            )
        return raw_name

    alias_matches: list[tuple[str, list[str]]] = []
    for alias in aliases:
        candidates = lookup.get(_normalise_column_name(alias), [])
        if candidates:
            alias_matches.append((alias, candidates))
    if not alias_matches:
        return None

    selected_alias, selected_candidates = alias_matches[0]
    selected = selected_candidates[0]
    for duplicate in selected_candidates[1:]:
        if _values_conflict(frame, selected, duplicate):
            raise SourceMappingError(
                f"conflicting columns for {canonical!r}: {selected!r} and {duplicate!r}"
            )
        warnings.append(
            f"multiple equivalent columns for {canonical!r}; kept {selected!r}"
        )
    if len(alias_matches) > 1:
        alternatives = [alias for alias, _ in alias_matches[1:]]
        warnings.append(
            f"multiple alias candidates for {canonical!r}; kept {selected_alias!r} "
            f"over {', '.join(alternatives)}"
        )
    return selected


def adapt_source_table(
    frame: pd.DataFrame,
    *,
    source_id: str,
    default_species: str | None = None,
    column_map: Mapping[str, str] | None = None,
) -> tuple[pd.DataFrame, AdapterReport]:
    """Adapt one source table to the project's canonical metadata schema.

    ``source_id`` is required and is written explicitly to every output row.
    ``default_species`` is optional and must be supplied by the caller; the
    adapter never guesses a species merely from a database name. ``column_map``
    can override alias matching when a downloaded file uses an unfamiliar
    header. Required fields raise :class:`SourceMappingError` when unavailable.
    """

    source_id = str(source_id).strip()
    if not source_id:
        raise SourceMappingError("source_id must not be empty")
    if default_species is not None and not str(default_species).strip():
        raise SourceMappingError("default_species must not be empty when supplied")
    explicit_map = dict(column_map or {})
    unknown_map_keys = sorted(set(explicit_map) - set(CANONICAL_COLUMNS))
    if unknown_map_keys:
        raise SourceMappingError(
            "column_map contains unknown canonical fields: " + ", ".join(unknown_map_keys)
        )

    lookup = _column_lookup(frame)
    warnings: list[str] = []
    mapped: dict[str, str] = {}
    for canonical, aliases in REQUIRED_SOURCE_ALIASES.items():
        raw_name = _resolve_column(
            frame, lookup, canonical, aliases, explicit_map, warnings
        )
        if raw_name is not None:
            mapped[canonical] = raw_name

    raw_source_record_id = _resolve_column(
        frame,
        lookup,
        "source_record_id",
        OPTIONAL_SOURCE_ALIASES["source_record_id"],
        explicit_map,
        warnings,
    )
    if raw_source_record_id is not None:
        mapped["source_record_id"] = raw_source_record_id

    missing_required = [
        field
        for field in REQUIRED_SOURCE_FIELDS
        if field not in mapped and (field != "species" or default_species is None)
    ]
    if missing_required:
        raise SourceMappingError(
            "source table is missing required fields: " + ", ".join(missing_required)
        )

    output = pd.DataFrame(index=frame.index)
    output["sequence_id"] = frame[mapped["sequence_id"]]
    if "species" in mapped:
        output["species"] = frame[mapped["species"]]
    else:
        output["species"] = default_species
        warnings.append("species was filled from the explicit default_species argument")
    output["source_dataset"] = source_id
    output["sequence"] = frame[mapped["sequence"]]

    for canonical, aliases in OPTIONAL_SOURCE_ALIASES.items():
        if canonical == "source_record_id":
            raw_name = raw_source_record_id
        else:
            raw_name = _resolve_column(
                frame, lookup, canonical, aliases, explicit_map, warnings
            )
            if raw_name is not None:
                mapped[canonical] = raw_name
        output[canonical] = frame[raw_name] if raw_name is not None else pd.NA

    derived_fields: list[str] = []
    if "source_record_id" not in mapped and "sequence_id" in mapped:
        warnings.append("source_record_id is unavailable; sequence_id was not relabeled as provenance")
    if "sequence_length" not in mapped:
        output["sequence_length"] = output["sequence"].map(normalize_sequence).str.len()
        derived_fields.append("sequence_length")

    output = output.loc[:, CANONICAL_COLUMNS].reset_index(drop=True)
    missing_optional = [
        field
        for field in CANONICAL_COLUMNS
        if field not in REQUIRED_SOURCE_FIELDS
        and field not in {"source_dataset", "sequence_length"}
        and field not in mapped
    ]
    report = AdapterReport(
        source_id=source_id,
        input_rows=len(frame),
        output_rows=len(output),
        mapped_columns=mapped,
        missing_required_fields=[],
        missing_optional_fields=missing_optional,
        derived_fields=derived_fields,
        warnings=warnings,
    )
    return output, report


def write_adapter_report(report: AdapterReport, output_path: str | Path) -> None:
    """Write an adapter report as UTF-8 JSON."""

    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(report.to_dict(), ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )


def parse_gff3_attributes(value: str) -> dict[str, str]:
    """Parse the simple ``key=value`` attributes used by source GFF3 files."""

    attributes: dict[str, str] = {}
    for item in str(value).split(";"):
        if not item.strip() or "=" not in item:
            continue
        key, raw_value = item.split("=", 1)
        attributes[unquote(key.strip())] = unquote(raw_value.strip())
    return attributes


def adapt_regulondb_gff3(
    path: str | Path,
    *,
    source_id: str = "regulondb",
    species: str,
) -> tuple[pd.DataFrame, AdapterReport]:
    """Adapt RegulonDB ``PromoterSet.gff3`` records to project metadata.

    RegulonDB encodes the local TSS base as the single uppercase nucleotide in
    its ``Sequence`` attribute. If that marker is absent or ambiguous, the
    local ``tss_position`` is left missing rather than inferred from a guessed
    window. The original GFF3 file remains the authoritative provenance record.
    """

    species = str(species).strip()
    if not species:
        raise SourceMappingError("species is required when adapting RegulonDB GFF3")
    rows: list[dict[str, Any]] = []
    seen_record_ids: dict[str, int] = {}
    input_path = Path(path)
    with input_path.open("r", encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n\r")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 9:
                raise SourceMappingError(
                    f"invalid GFF3 row at line {line_number}: expected 9 columns"
                )
            seqid, _source, feature_type, start, end, score, strand, phase, raw_attributes = fields
            attributes = parse_gff3_attributes(raw_attributes)
            if feature_type not in {"transcription_start_site", "promoter"}:
                continue
            raw_sequence = attributes.get("Sequence", "")
            uppercase_positions = [
                index + 1
                for index, nucleotide in enumerate(raw_sequence)
                if nucleotide in "ACGT"
            ]
            local_tss = uppercase_positions[0] if len(uppercase_positions) == 1 else pd.NA
            record_id = attributes.get("name") or f"{seqid}:{start}-{end}:{strand}:{line_number}"
            occurrence = seen_record_ids.get(record_id, 0) + 1
            seen_record_ids[record_id] = occurrence
            sequence_id = (
                record_id if occurrence == 1 else f"{record_id}__duplicate_{occurrence}"
            )
            rows.append(
                {
                    "sequence_id": sequence_id,
                    "species": species,
                    "sequence": raw_sequence.upper(),
                    "taxon_id": pd.NA,
                    "assembly_accession": seqid,
                    "sequence_length": len(raw_sequence) if raw_sequence else pd.NA,
                    "tss_position": local_tss,
                    "strand": strand,
                    "sigma_factor_type": attributes.get("SigmaFactor") or pd.NA,
                    "known_element_annotations": pd.NA,
                    "promoter_strength": pd.NA,
                    "evidence_level": attributes.get("Confidence") or pd.NA,
                    "source_record_id": record_id,
                    "source_tss_coordinate": start,
                    "source_tss_end_coordinate": end,
                    "source_evidence": attributes.get("Evidence") or pd.NA,
                    "source_feature_type": feature_type,
                    "source_score": score,
                    "source_phase": phase,
                }
            )
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise SourceMappingError(f"no promoter records found in GFF3: {input_path}")
    adapted, report = adapt_source_table(
        frame,
        source_id=source_id,
        default_species=species,
        column_map={
            "sequence_id": "sequence_id",
            "species": "species",
            "sequence": "sequence",
            "taxon_id": "taxon_id",
            "assembly_accession": "assembly_accession",
            "sequence_length": "sequence_length",
            "tss_position": "tss_position",
            "strand": "strand",
            "sigma_factor_type": "sigma_factor_type",
            "known_element_annotations": "known_element_annotations",
            "promoter_strength": "promoter_strength",
            "evidence_level": "evidence_level",
            "source_record_id": "source_record_id",
        },
    )
    for extra_column in (
        "source_tss_coordinate",
        "source_tss_end_coordinate",
        "source_evidence",
        "source_feature_type",
        "source_score",
        "source_phase",
    ):
        adapted[extra_column] = frame[extra_column].reset_index(drop=True)
    return adapted, report


def main() -> None:
    """Adapt a CSV/TSV source file from the command line."""

    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, help="Source CSV/TSV file")
    parser.add_argument("--output", required=True, help="Canonical metadata TSV")
    parser.add_argument("--source-id", required=True, help="Audited source identifier")
    parser.add_argument(
        "--default-species",
        default=None,
        help="Explicit species value only when the source file has no species column",
    )
    parser.add_argument(
        "--column-map",
        default=None,
        help='JSON object mapping canonical fields to raw columns, e.g. {"sequence":"DNA"}',
    )
    parser.add_argument(
        "--report",
        default=None,
        help="Adapter report path; defaults to <output>.adapter_report.json",
    )
    args = parser.parse_args()
    input_path = Path(args.input)
    if input_path.suffix.lower() == ".gff3":
        if not args.default_species:
            raise SourceMappingError(
                "--default-species is required for GFF3 because species must be explicit"
            )
        adapted, report = adapt_regulondb_gff3(
            input_path,
            source_id=args.source_id,
            species=args.default_species,
        )
    else:
        frame = pd.read_csv(input_path, sep="," if input_path.suffix.lower() == ".csv" else "\t")
        column_map = json.loads(args.column_map) if args.column_map else None
        adapted, report = adapt_source_table(
            frame,
            source_id=args.source_id,
            default_species=args.default_species,
            column_map=column_map,
        )
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    adapted.to_csv(output_path, sep="\t", index=False)
    report_path = Path(args.report) if args.report else output_path.with_suffix(".adapter_report.json")
    write_adapter_report(report, report_path)
    print(f"input_rows={report.input_rows}")
    print(f"output_rows={report.output_rows}")
    print(f"output={output_path.resolve()}")
    print(f"report={report_path.resolve()}")


if __name__ == "__main__":
    main()
