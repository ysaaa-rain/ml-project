import pandas as pd

from preprocessing.recover_dbtbs_coordinates import (
    _recover_coordinate_fields,
    extract_recovered_promoter_rows,
)


def test_recover_negative_strand_tss_from_dbtbs_location_and_absolute_position():
    recovered = _recover_coordinate_fields(
        location="-43:+15",
        absolute_position="2518244..2518301",
        sequence="A" * 58,
        strand="-",
    )

    assert recovered["coordinate_status"] == "recovered_from_dbtbs"
    assert recovered["tss_position"] == 44
    assert recovered["source_tss_coordinate"] == 2518258
    assert recovered["source_window_start"] == 2518244
    assert recovered["source_window_end"] == 2518301


def test_recover_positive_strand_tss_uses_forward_coordinate():
    recovered = _recover_coordinate_fields(
        location="-10:+5",
        absolute_position="100..114",
        sequence="A" * 15,
        strand="+",
    )

    assert recovered["coordinate_status"] == "recovered_from_dbtbs"
    assert recovered["tss_position"] == 11
    assert recovered["source_tss_coordinate"] == 110


def test_location_nd_is_not_guessed():
    recovered = _recover_coordinate_fields(
        location="ND",
        absolute_position="100..120",
        sequence="A" * 21,
        strand="+",
    )

    assert recovered["coordinate_status"] == "unresolved"
    assert recovered["coordinate_issue"] == "location_nd_or_missing"
    assert pd.isna(recovered["tss_position"])
    assert pd.isna(recovered["source_tss_coordinate"])


def test_inconsistent_spans_are_rejected():
    recovered = _recover_coordinate_fields(
        location="-10:+5",
        absolute_position="100..114",
        sequence="A" * 14,
        strand="+",
    )

    assert recovered["coordinate_status"] == "unresolved"
    assert recovered["coordinate_issue"] == "location_length_mismatch"


def test_page_parser_recovers_gene_direction_and_tss():
    html = """
    <div class="table_title">Genes</div>
    <table>
      <tr><th>Genes</th><th>Synonyms</th><th>Direction</th><th>Genome position</th><th>Function</th></tr>
      <tr><td>geneA</td><td></td><td>-</td><td>200..300</td><td>x</td></tr>
    </table>
    <h2>Promoters</h2>
    <table>
      <tr><th>Binding<br>factor</th><th>Regulation</th><th>Location</th>
          <th>Absolute position</th><th>Binding seq.(cis-element)</th>
          <th>Experimental evidence</th></tr>
      <tr><td>SigA</td><td>Promoter</td><td><tt>-4:+4</tt></td>
          <td>400..407</td><td><tt>AAA<b>T</b>CCCC</tt></td><td>Paper: S1</td></tr>
    </table>
    """

    records, counts = extract_recovered_promoter_rows(html, page_name="geneA.html")

    assert len(records) == 1
    record = records[0]
    assert record["strand"] == "-"
    assert record["tss_position"] == 5
    assert record["source_tss_coordinate"] == 403
    assert record["coordinate_status"] == "recovered_from_dbtbs"
    assert record["assembly_accession"] == "NC_000964"
    assert counts["coordinates_recovered"] == 1
