import json

from experiments.render_meme_logos import parse_meme_xml, render


def _write_meme_xml(path):
    path.write_text(
        """<MEME version="5.5.9"><motifs><motif id="motif_1" name="AT" alt="MEME-1" width="2" sites="12" e_value="1e-5">
        <probabilities><alphabet_matrix>
        <alphabet_array><value letter_id="A">0.8</value><value letter_id="C">0.1</value><value letter_id="G">0.05</value><value letter_id="T">0.05</value></alphabet_array>
        <alphabet_array><value letter_id="A">0.05</value><value letter_id="C">0.05</value><value letter_id="G">0.1</value><value letter_id="T">0.8</value></alphabet_array>
        </alphabet_matrix></probabilities></motif></motifs></MEME>""",
        encoding="utf-8",
    )


def test_parse_meme_xml_and_render_report_logos(tmp_path):
    meme_root = tmp_path / "meme"
    (meme_root / "bacillus_subtilis").mkdir(parents=True)
    _write_meme_xml(meme_root / "bacillus_subtilis" / "meme.xml")

    motifs = parse_meme_xml(meme_root / "bacillus_subtilis" / "meme.xml")
    manifest = render(meme_root, tmp_path / "results", ["bacillus_subtilis"])

    assert len(motifs) == 1
    assert motifs[0]["consensus"] == "AT"
    assert manifest["species_motif_counts"] == {"bacillus_subtilis": 1}
    assert manifest["status"] == "partial_validation_subset"
    assert (tmp_path / "results" / "bacillus_subtilis_motif_logos.png").stat().st_size > 1000
    assert (tmp_path / "results" / "six_species_motif_logos.pdf").stat().st_size > 1000
    saved = json.loads((tmp_path / "results" / "logo_manifest.json").read_text(encoding="utf-8"))
    assert "input_xml_sha256" in saved
