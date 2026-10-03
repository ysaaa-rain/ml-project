import pytest

from experiments.summarize_motif_comparison import parse_discovered_motifs, read_tomtom, validate_reference_consensus


def test_parse_discovered_motifs_extracts_consensus_and_statistics(tmp_path):
    path = tmp_path / "meme.txt"
    path.write_text(
        "MOTIF TATAAT MEME-1\twidth =  6 sites = 123 llr = 200 p-value = 1e-8 E-value = 2.5e-6\n",
        encoding="utf-8",
    )

    assert parse_discovered_motifs(path) == {
        "TATAAT": {"meme_id": "MEME-1", "width": 6, "sites": 123, "meme_evalue": 2.5e-6}
    }


def test_read_tomtom_skips_footer_comments(tmp_path):
    path = tmp_path / "tomtom.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        handle.write("Query_ID\tTarget_ID\tOptimal_offset\tp-value\tE-value\tq-value\tOverlap\tQuery_consensus\tTarget_consensus\tOrientation\n")
        handle.write("TATAAT\tminus10\t0\t0.1\t0.4\t0.2\t6\tTATAAT\tTATAAT\t+\n")
        handle.write("\n# Tomtom version footer\n# command footer\n")

    rows = read_tomtom(path)

    assert len(rows) == 1
    assert rows[0]["Query_ID"] == "TATAAT"
    assert rows[0]["Target_ID"] == "minus10"


def test_reference_validator_checks_named_core_box_sequences():
    from pathlib import Path

    path = Path(__file__).parents[1] / "baselines" / "known_promoter_elements.meme"

    assert validate_reference_consensus(path) == {
        "Ecoli_sigma70_minus10": "TATAAT",
        "Ecoli_sigma70_minus35": "TTGACA",
    }


def test_reference_validator_rejects_mislabeled_consensus(tmp_path):
    path = tmp_path / "wrong.meme"
    path.write_text(
        "MOTIF Ecoli_sigma70_minus10\nletter-probability matrix: alength= 4 w= 6\n"
        "1 0 0 0\n0 0 0 1\n1 0 0 0\n0 0 0 1\n0 0 0 1\n0 1 0 0\n"
        "MOTIF Ecoli_sigma70_minus35\nletter-probability matrix: alength= 4 w= 6\n"
        "0 0 0 1\n0 0 0 1\n0 0 0 1\n1 0 0 0\n0 0 1 0\n1 0 0 0\n",
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="expected consensus TATAAT, got ATATTC"):
        validate_reference_consensus(path)
