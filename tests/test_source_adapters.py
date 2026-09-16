import pandas as pd
import pytest

from preprocessing.pipeline import normalize_metadata
from preprocessing.source_adapters import SourceMappingError, adapt_source_table


def test_adapter_maps_common_source_headers_and_preserves_record_id():
    source = pd.DataFrame(
        {
            "Promoter ID": ["p1", "p2"],
            "Organism": ["Escherichia coli", "Escherichia coli"],
            "DNA sequence": [" au gc ", "TTGACA"],
            "Sigma": ["RpoD", "RpoS"],
            "Evidence": ["experimental", "inferred"],
            "TSS coordinate": [5, 5],
            "Orientation": ["+", "-"],
        }
    )

    adapted, report = adapt_source_table(source, source_id="regulondb")

    assert adapted["sequence_id"].tolist() == ["p1", "p2"]
    assert adapted["source_dataset"].tolist() == ["regulondb", "regulondb"]
    assert adapted["sequence"].tolist() == [" au gc ", "TTGACA"]
    assert adapted["sequence_length"].tolist() == [4, 6]
    assert adapted["sigma_factor_type"].tolist() == ["RpoD", "RpoS"]
    assert adapted["tss_position"].tolist() == [5, 5]
    assert adapted["source_record_id"].tolist() == ["p1", "p2"]
    assert report.missing_optional_fields == [
        "taxon_id",
        "assembly_accession",
        "known_element_annotations",
        "promoter_strength",
    ]


def test_adapter_uses_explicit_default_species_and_keeps_missing_fields_missing():
    source = pd.DataFrame(
        {
            "record_id": ["dbtbs-1"],
            "seq": ["ACGT"],
        }
    )

    adapted, report = adapt_source_table(
        source,
        source_id="dbtbs",
        default_species="Bacillus subtilis",
    )

    assert adapted.loc[0, "sequence_id"] == "dbtbs-1"
    assert adapted.loc[0, "source_record_id"] == "dbtbs-1"
    assert adapted.loc[0, "species"] == "Bacillus subtilis"
    assert pd.isna(adapted.loc[0, "sigma_factor_type"])
    assert "species was filled" in report.warnings[0]


def test_adapter_output_can_enter_existing_normalization_pipeline():
    source = pd.DataFrame(
        {
            "promoter_id": ["p1", "p2"],
            "species_name": ["Escherichia coli", "Escherichia coli"],
            "promoter_sequence": ["AUGC", "ATGX"],
            "record_id": ["r1", "r2"],
            "evidence_type": ["experimental", "inferred"],
        }
    )
    adapted, _ = adapt_source_table(source, source_id="regulondb")
    clean, report = normalize_metadata(adapted, max_n_fraction=0.1)

    assert len(clean) == 1
    assert clean.loc[0, "sequence"] == "ATGC"
    assert clean.loc[0, "source_record_id"] == "r1"
    assert report.invalid_row_count == 1


def test_adapter_rejects_missing_required_field_without_guessing():
    source = pd.DataFrame({"id": ["p1"], "sequence": ["ACGT"]})

    with pytest.raises(SourceMappingError, match="species"):
        adapt_source_table(source, source_id="regulondb")


def test_adapter_rejects_conflicting_alias_columns():
    source = pd.DataFrame(
        {
            "id": ["p1"],
            "sequence_id": ["p1"],
            "Sequence ID": ["different"],
            "species": ["Escherichia coli"],
            "sequence": ["ACGT"],
        }
    )

    with pytest.raises(SourceMappingError, match="conflicting columns"):
        adapt_source_table(source, source_id="regulondb")
