import pandas as pd
import pytest

from preprocessing.pipeline import normalize_metadata
from preprocessing.assemble_metadata import assemble_source_tables
from preprocessing.download_dbtbs import extract_promoter_rows, parse_index_links
from preprocessing.source_adapters import (
    SourceMappingError,
    adapt_regulondb_gff3,
    adapt_source_table,
)


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


def test_regulondb_gff3_adapter_reads_sequence_sigma_evidence_and_tss(tmp_path):
    gff3 = tmp_path / "PromoterSet.gff3"
    gff3.write_text(
        "NC_000913.3\tRegulonDB\ttranscription_start_site\t101\t101\t.\t+\t.\t"
        "name=geneP;SigmaFactor=Sigma70;Sequence=acgtA;Evidence=[TIM|S|mapping];Confidence=Strong\n",
        encoding="utf-8",
    )

    adapted, report = adapt_regulondb_gff3(
        gff3,
        species="Escherichia coli",
    )

    assert len(adapted) == 1
    assert adapted.loc[0, "sequence_id"] == "geneP"
    assert adapted.loc[0, "sequence"] == "ACGTA"
    assert adapted.loc[0, "sigma_factor_type"] == "Sigma70"
    assert adapted.loc[0, "evidence_level"] == "Strong"
    assert adapted.loc[0, "tss_position"] == 5
    assert adapted.loc[0, "source_tss_coordinate"] == "101"
    assert report.source_id == "regulondb"


def test_regulondb_gff3_adapter_makes_sequence_ids_unique_without_losing_source_ids(tmp_path):
    gff3 = tmp_path / "PromoterSet.gff3"
    row = (
        "NC_000913.3\tRegulonDB\ttranscription_start_site\t101\t101\t.\t+\t.\t"
        "name=geneP;SigmaFactor=Sigma70;Sequence=acgtA;Confidence=Strong\n"
    )
    gff3.write_text(row + row, encoding="utf-8")

    adapted, _ = adapt_regulondb_gff3(gff3, species="Escherichia coli")

    assert adapted["sequence_id"].tolist() == ["geneP", "geneP__duplicate_2"]
    assert adapted["source_record_id"].tolist() == ["geneP", "geneP"]


def test_dbtbs_parser_keeps_sigma_promoters_and_skips_tf_binding_rows():
    html = """
    <h2>Promoters</h2>
    <table>
      <tr><th>Binding<br>factor</th><th>Regulation</th><th>Location</th>
          <th>Absolute position</th><th>Binding seq.(cis-element)</th>
          <th>Experimental evidence</th></tr>
      <tr><td>SigA</td><td>Promoter</td><td><tt>-42:+5</tt></td>
          <td>100..146</td><td><tt>TTGACA<font color=red>AA</font>TATAAT</tt></td>
          <td>Paper A: PE<br>Paper B: S1</td></tr>
          <tr><td>LexA</td><td>Promoter</td><td>ND</td><td>200..210</td>
          <td><tt>ACGT</tt></td><td>Paper C: DB</td></tr>
      <tr><td>SigB</td><td>Promoter</td><td>ND</td><td>300..303</td>
          <td><tt>ND</tt></td><td>Paper D: PE</td></tr>
    </table>
    """

    records, counts = extract_promoter_rows(html, page_name="operon.html")

    assert len(records) == 1
    assert records[0]["sigma_factor_type"] == "SigA"
    assert records[0]["sequence"] == "TTGACAAATATAAT"
    assert records[0]["evidence_level"] == "experimental"
    assert records[0]["source_tss_coordinate"] == 100
    assert counts["promoter_rows_seen"] == 3
    assert counts["non_sigma_promoters_skipped"] == 1
    assert counts["invalid_or_missing_sequence_skipped"] == 1


def test_dbtbs_index_parser_deduplicates_operon_links():
    html = """
    <a href="COG/prom/a.html">a</a>
    <a href="COG/prom/a.html">a again</a>
    <a href="COG/prom/b.html">b</a>
    """

    assert parse_index_links(html, base_url="https://dbtbs.hgc.jp/ver4/") == [
        "https://dbtbs.hgc.jp/ver4/COG/prom/a.html",
        "https://dbtbs.hgc.jp/ver4/COG/prom/b.html",
    ]


def test_assemble_source_tables_tracks_hashes_and_rejects_duplicate_ids(tmp_path):
    first = tmp_path / "first.tsv"
    second = tmp_path / "second.tsv"
    frame = pd.DataFrame(
        {
            "sequence_id": ["a"],
            "species": ["A"],
            "source_dataset": ["source_a"],
            "sequence": ["ACGT"],
        }
    )
    frame.to_csv(first, sep="\t", index=False)
    frame.assign(sequence_id="b", source_dataset="source_b").to_csv(second, sep="\t", index=False)
    output = tmp_path / "combined.tsv"
    report_path = tmp_path / "combined.json"

    report = assemble_source_tables(
        [first, second],
        output_path=output,
        report_path=report_path,
    )

    assert report["output_rows"] == 2
    assert report["output_sha256"]
    assert len(report["input_tables"]) == 2
    assert pd.read_csv(output, sep="\t")["sequence_id"].tolist() == ["a", "b"]

    frame.assign(sequence_id="a", source_dataset="source_c").to_csv(second, sep="\t", index=False)
    with pytest.raises(ValueError, match="globally unique"):
        assemble_source_tables([first, second], output_path=output, report_path=report_path)
