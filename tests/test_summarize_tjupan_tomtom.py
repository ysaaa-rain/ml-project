import csv
import json
from pathlib import Path

import pytest

from experiments.summarize_tjupan_tomtom import REFERENCE_IDS, TOMTOM_COLUMNS, summarize


def _prepare(tmp_path, *, omit_target=None):
    matrix_dir = tmp_path / "matrix_summary"
    (matrix_dir / "matrices").mkdir(parents=True)
    (matrix_dir / "matrices" / "bacillus_subtilis.meme").write_text("MEME matrix", encoding="utf-8")
    (matrix_dir / "motif_summary.tsv").write_text(
        "species\tmotif_id\tconsensus\twidth_bp\tsites_reported\te_value\n"
        "bacillus_subtilis\tMEME-1\tATCGHC\t6\t100\t0.01\n",
        encoding="utf-8",
    )
    (matrix_dir / "summary_manifest.json").write_text(json.dumps({
        "included_species": ["bacillus_subtilis"],
        "species_motif_counts": {"bacillus_subtilis": 1},
    }), encoding="utf-8")

    tomtom_dir = tmp_path / "tomtom" / "bacillus_subtilis"
    tomtom_dir.mkdir(parents=True)
    tsv = tomtom_dir / "tomtom.tsv"
    with tsv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=TOMTOM_COLUMNS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for target in sorted(REFERENCE_IDS - ({omit_target} if omit_target else set())):
            writer.writerow({
                "Query_ID": "bacillus_subtilis_MEME-1", "Target_ID": target,
                "Optimal_offset": "0", "p-value": "0.1", "E-value": "0.4", "q-value": "0.4",
                "Overlap": "6", "Query_consensus": "ATCGHC", "Target_consensus": "TATAAT", "Orientation": "+",
            })
    return matrix_dir, tmp_path / "tomtom"


def test_summarize_tomtom_marks_partial_and_exports_all_pairs(tmp_path):
    matrix_dir, tomtom_root = _prepare(tmp_path)
    reference = Path(__file__).parents[1] / "baselines" / "known_promoter_elements.meme"

    manifest = summarize(tomtom_root, matrix_dir, reference, tmp_path / "results", ["bacillus_subtilis"])

    assert manifest["status"] == "partial_validation_subset"
    assert manifest["pair_count"] == 4
    best = (tmp_path / "results" / "best_reference_match_per_motif.tsv").read_text(encoding="utf-8")
    assert "best_reference_match_descriptive_only" in best
    assert "TATAAT" not in best  # target consensus is not copied into the compact best-match table


def test_summarize_tomtom_rejects_missing_reference_match(tmp_path):
    matrix_dir, tomtom_root = _prepare(tmp_path, omit_target="Ecoli_UP_distal")
    reference = Path(__file__).parents[1] / "baselines" / "known_promoter_elements.meme"

    with pytest.raises(ValueError, match="expected exactly 4 query-reference pairs"):
        summarize(tomtom_root, matrix_dir, reference, tmp_path / "results", ["bacillus_subtilis"])
