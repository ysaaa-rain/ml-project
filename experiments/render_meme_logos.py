"""Render MEME XML probability matrices as information-content sequence logos."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET
from typing import Any

from experiments.prepare_tjupan_motif_inputs import SPECIES

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.font_manager import FontProperties, findfont
from matplotlib.patches import PathPatch
from matplotlib.textpath import TextPath
from matplotlib.transforms import Affine2D


LETTERS = ("A", "C", "G", "T")
COLORS = {"A": "#2a9d55", "C": "#2878b5", "G": "#e49b25", "T": "#d1495b"}
GLYPH_FONT = FontProperties(family="DejaVu Sans", weight="bold")


def _has_font(family: str) -> bool:
    try:
        findfont(FontProperties(family=family), fallback_to_default=False)
        return True
    except ValueError:
        return False


CJK_FONT = next(
    (family for family in ("Hiragino Sans GB", "PingFang SC", "Noto Sans CJK SC", "Microsoft YaHei", "SimHei") if _has_font(family)),
    None,
)
TEXT_FONT = FontProperties(family=CJK_FONT or "DejaVu Sans")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_meme_xml(path: Path) -> list[dict[str, Any]]:
    root = ET.parse(path).getroot()
    motifs: list[dict[str, Any]] = []
    for node in root.findall(".//motif"):
        matrix = node.find("./probabilities/alphabet_matrix")
        if matrix is None:
            continue
        columns: list[dict[str, float]] = []
        for array in matrix.findall("./alphabet_array"):
            column = {value.attrib["letter_id"]: float(value.text or "nan") for value in array.findall("./value")}
            if set(column) != set(LETTERS) or any(not math.isfinite(value) or value < 0 or value > 1 for value in column.values()):
                raise ValueError(f"invalid probability column in {path}: {column}")
            if not math.isclose(sum(column.values()), 1.0, rel_tol=0, abs_tol=2e-5):
                raise ValueError(f"probabilities do not sum to one in {path}: {column}")
            columns.append(column)
        width = int(node.attrib["width"])
        if len(columns) != width:
            raise ValueError(f"motif width {width} != probability columns {len(columns)} in {path}")
        motifs.append({
            "id": node.attrib.get("alt", node.attrib["id"]),
            "consensus": node.attrib.get("name", ""),
            "width": width,
            "sites": int(node.attrib.get("sites", "0")),
            "p_value": float(node.attrib.get("p_value", "nan")),
            "e_value": float(node.attrib.get("e_value", "nan")),
            "columns": columns,
        })
    if not motifs:
        raise ValueError(f"no probability matrices found in {path}")
    return motifs


def draw_logo(axis: Any, motif: dict[str, Any], *, show_ylabel: bool) -> None:
    for position, column in enumerate(motif["columns"]):
        entropy = -sum(p * math.log2(p) for p in column.values() if p > 0)
        information = max(0.0, 2.0 - entropy)
        baseline = 0.0
        for letter in sorted(LETTERS, key=lambda item: column[item]):
            height = column[letter] * information
            if height <= 0.002:
                continue
            glyph = TextPath((0, 0), letter, size=1, prop=GLYPH_FONT)
            bounds = glyph.get_extents()
            width = 0.82
            scale_x = width / bounds.width
            scale_y = height / bounds.height
            transform = (
                Affine2D()
                .scale(scale_x, scale_y)
                .translate(position + 0.5 - (bounds.x0 + bounds.width / 2) * scale_x,
                           baseline - bounds.y0 * scale_y)
            )
            axis.add_patch(PathPatch(glyph.transformed(transform), facecolor=COLORS[letter], edgecolor="none"))
            baseline += height

    axis.set_xlim(0, motif["width"])
    axis.set_ylim(0, 2.08)
    axis.set_xticks(range(1, motif["width"] + 1, 2))
    axis.tick_params(axis="both", labelsize=7, length=2)
    axis.set_xlabel("Motif 内位置（bp）" if CJK_FONT else "Motif position (bp)", fontsize=8, fontproperties=TEXT_FONT)
    if show_ylabel:
        axis.set_ylabel("信息量（bit）" if CJK_FONT else "Information (bits)", fontsize=8, fontproperties=TEXT_FONT)
        axis.set_yticks((0, 1, 2))
    else:
        axis.set_yticks((0, 1, 2))
        axis.set_yticklabels([])
    axis.set_title(
        f"{motif['id']}  {motif['consensus']}  |  sites={motif['sites']}  |  E={motif['e_value']:.2g}",
        fontsize=9,
        loc="left",
        fontproperties=TEXT_FONT,
    )
    axis.spines["top"].set_visible(False)
    axis.spines["right"].set_visible(False)


def render(meme_root: Path, output_dir: Path, species: list[str]) -> dict[str, Any]:
    if not species or len(set(species)) != len(species) or any(slug not in SPECIES for slug in species):
        raise ValueError(f"invalid or repeated species list: {species}")
    output_dir.mkdir(parents=True, exist_ok=True)
    pdf_path = output_dir / "six_species_motif_logos.pdf"
    sources: dict[str, str] = {}
    summaries: dict[str, int] = {}
    png_paths: list[Path] = []

    with PdfPages(pdf_path) as pdf:
        for slug in species:
            xml_path = meme_root / slug / "meme.xml"
            motifs = parse_meme_xml(xml_path)
            sources[f"{slug}/meme.xml"] = sha256_file(xml_path)
            summaries[slug] = len(motifs)
            rows = math.ceil(len(motifs) / 2)
            fig, axes = plt.subplots(rows, 2, figsize=(13.5, max(3.2, rows * 1.85)), squeeze=False)
            flat_axes = list(axes.flat)
            for index, motif in enumerate(motifs):
                draw_logo(flat_axes[index], motif, show_ylabel=index % 2 == 0)
            for axis in flat_axes[len(motifs):]:
                axis.axis("off")
            title = (
                f"{slug}：MEME motif 序列 logo（信息量刻度）"
                if CJK_FONT
                else f"{slug}: MEME motif sequence logos (information content)"
            )
            fig.suptitle(title, fontsize=14, y=0.998, fontproperties=TEXT_FONT)
            fig.tight_layout(rect=(0, 0, 1, 0.98))
            png_path = output_dir / f"{slug}_motif_logos.png"
            fig.savefig(png_path, dpi=180, bbox_inches="tight")
            pdf.savefig(fig, bbox_inches="tight")
            plt.close(fig)
            png_paths.append(png_path)

    outputs = [pdf_path, *png_paths]
    manifest = {
        "purpose": "从 MEME XML 概率矩阵绘制 motif 序列 logo；未读取或输出原始序列",
        "method": "每列信息量 = 2 - Shannon entropy；字母高度 = 碱基概率 × 列信息量；未作有限样本修正",
        "status": "complete" if species == list(SPECIES) else "partial_validation_subset",
        "expected_species": list(SPECIES),
        "included_species": species,
        "species_motif_counts": summaries,
        "input_xml_sha256": sources,
        "outputs": {path.name: {"bytes": path.stat().st_size, "sha256": sha256_file(path)} for path in outputs},
    }
    (output_dir / "logo_manifest.json").write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--meme-root", type=Path, required=True, help="包含六个物种子目录和 meme.xml 的目录")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--species", nargs="+", required=True)
    args = parser.parse_args()
    result = render(args.meme_root, args.output_dir, args.species)
    print(f"rendered {sum(result['species_motif_counts'].values())} motif logos for {len(result['species_motif_counts'])} groups")


if __name__ == "__main__":
    main()
