import hashlib
import json

from experiments.summarize_tjupan_meme import export


def _write_xml(path):
    path.write_text(
        """<MEME version="5.5.9"><training_set><background_frequencies><alphabet_array>
        <value letter_id="A">0.3</value><value letter_id="C">0.2</value><value letter_id="G">0.2</value><value letter_id="T">0.3</value>
        </alphabet_array></background_frequencies></training_set><motifs><motif id="motif_1" name="AT" alt="MEME-1" width="2" sites="12" p_value="1e-5" e_value="2e-5">
        <probabilities><alphabet_matrix>
        <alphabet_array><value letter_id="A">0.8</value><value letter_id="C">0.1</value><value letter_id="G">0.05</value><value letter_id="T">0.05</value></alphabet_array>
        <alphabet_array><value letter_id="A">0.05</value><value letter_id="C">0.05</value><value letter_id="G">0.1</value><value letter_id="T">0.8</value></alphabet_array>
        </alphabet_matrix></probabilities></motif></motifs></MEME>""",
        encoding="utf-8",
    )


def test_export_writes_sequence_free_matrices_and_hashed_summary(tmp_path):
    meme_root = tmp_path / "meme"
    group = meme_root / "bacillus_subtilis"
    group.mkdir(parents=True)
    _write_xml(group / "meme.xml")
    (group / "meme.txt").write_text("MOTIF AT MEME-1 width = 2\n", encoding="utf-8")
    output_hashes = {
        name: hashlib.sha256((group / name).read_bytes()).hexdigest()
        for name in ("meme.xml", "meme.txt")
    }
    run_manifest = tmp_path / "run_manifest.json"
    run_manifest.write_text(json.dumps({
        "meme_version": "MEME version 5.5.9",
        "parameters": {"-dna": True},
        "input_manifest_sha256": "input-hash",
        "results": {"bacillus_subtilis": {
            "status": "success", "return_code": 0, "positive_count": 691, "control_count": 809,
            "positive_sha256": "positive-hash", "control_sha256": "control-hash",
            "output_files": {name: {"sha256": digest} for name, digest in output_hashes.items()},
        }},
    }), encoding="utf-8")

    manifest = export(meme_root, run_manifest, tmp_path / "summary", ["bacillus_subtilis"])

    matrix = (tmp_path / "summary" / "matrices" / "bacillus_subtilis.meme").read_text(encoding="utf-8")
    summary = (tmp_path / "summary" / "motif_summary.tsv").read_text(encoding="utf-8")
    assert "0.8 0.1 0.05 0.05" in matrix
    assert "seq_id" not in matrix and "ACGTAC" not in matrix
    assert "mean_information_bits" in summary
    assert manifest["species_motif_counts"] == {"bacillus_subtilis": 1}
    assert manifest["status"] == "partial_validation_subset"
    assert "matrices/bacillus_subtilis.meme" in manifest["output_sha256"]
