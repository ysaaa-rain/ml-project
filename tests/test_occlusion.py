import json

import pandas as pd

from experiments.select_occlusion_representatives import run as select_representatives
from experiments.summarize_occlusion import run as summarize_occlusion


def test_select_representatives_is_balanced_and_deterministic(tmp_path):
    input_path = tmp_path / "input.tsv"
    output_path = tmp_path / "selected.tsv"
    frame = pd.DataFrame(
        {
            "sequence_id": ["b2", "a2", "b1", "a1", "a3", "b3"],
            "sequence": ["CCCC", "AAAA", "CCCA", "AAAC", "AAAG", "CCCG"],
            "sigma_factor_type": ["B", "A", "B", "A", "A", "B"],
        }
    )
    frame.to_csv(input_path, sep="\t", index=False)

    manifest = select_representatives(input_path, output_path, per_group=2)
    selected = pd.read_csv(output_path, sep="\t")

    assert manifest["sequence_count"] == 4
    assert selected["sequence_id"].tolist() == ["a1", "a2", "b1", "b2"]
    assert selected["sigma_factor_type"].value_counts().to_dict() == {"A": 2, "B": 2}


def test_summarize_occlusion_uses_tss_relative_top_k(tmp_path):
    input_path = tmp_path / "input.tsv"
    windows_path = tmp_path / "windows.tsv"
    output_dir = tmp_path / "summary"
    metadata = pd.DataFrame(
        {
            "sequence_id": ["a", "b"],
            "sigma_factor_type": ["SigmaA", "SigmaB"],
            "tss_offset_in_window": [5, 5],
        }
    )
    windows = pd.DataFrame(
        {
            "sequence_id": ["a", "a", "a", "b", "b"],
            "start": [8, 2, 4, 7, 1],
            "end": [10, 4, 6, 9, 3],
            "importance": [0.2, 0.8, 0.4, 0.3, 0.1],
            "original_fragment": ["CC", "AA", "GG", "TT", "AC"],
        }
    )
    metadata.to_csv(input_path, sep="\t", index=False)
    windows.to_csv(windows_path, sep="\t", index=False)

    result = summarize_occlusion(input_path, windows_path, output_dir, top_k=2)
    top = pd.read_csv(output_dir / "top_occlusion_windows.tsv", sep="\t")
    payload = json.loads((output_dir / "occlusion_summary.json").read_text(encoding="utf-8"))

    assert result["top_window_count"] == 4
    assert top.loc[top["sequence_id"].eq("a"), "relative_start"].tolist() == [-3, -1]
    assert payload["top_k_per_sequence"] == 2
    assert "No resampling stability" in payload["boundary"]
