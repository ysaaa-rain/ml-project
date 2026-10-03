import math

from experiments.run_tjupan_fimo import _motif_pi0_estimates, _prepare_scan_inputs
from experiments.summarize_tjupan_fimo import _read_hits, benjamini_hochberg, fisher_exact_two_sided


def test_prepare_scan_inputs_uses_controls_for_zero_order_background(tmp_path):
    positive = tmp_path / "positive.fasta"
    control = tmp_path / "control.fasta"
    combined = tmp_path / "combined.fasta"
    background = tmp_path / "background"
    sequence = "ACGT" * 20 + "A"
    positive.write_text(f">p1|label=1\n{sequence}\n", encoding="ascii")
    control.write_text(f">n1|label=0\n{sequence}\n>n2|label=0\n{sequence}\n", encoding="ascii")

    counts = _prepare_scan_inputs(positive, control, combined, background)

    assert counts["positive_records"] == 1
    assert counts["control_records"] == 2
    assert counts["positive_unique_sequences"] == 1
    assert counts["control_unique_sequences"] == 1
    assert counts["cross_label_exact_sequence_overlap"] == 1
    assert math.isclose(sum(counts["background_frequencies"].values()), 1.0)
    assert combined.read_text(encoding="ascii").count(">") == 3
    assert background.read_text(encoding="ascii").splitlines() == [
        f"{base} {counts['background_frequencies'][base]:.12g}" for base in "ACGT"
    ]


def test_fisher_exact_two_sided_known_example():
    assert math.isclose(fisher_exact_two_sided(1, 9, 11, 3), 0.0027594561852200836, rel_tol=1e-10)


def test_benjamini_hochberg_preserves_input_order_and_monotonicity():
    assert benjamini_hochberg([0.01, 0.04, 0.03, 0.002]) == [0.02, 0.04, 0.04, 0.008]


def test_read_hits_accepts_fimo_output_with_no_significant_matches(tmp_path):
    fimo_tsv = tmp_path / "fimo.tsv"
    fimo_tsv.write_text(
        "\n# FIMO (Find Individual Motif Occurrences): Version 5.5.9\n",
        encoding="utf-8",
    )

    assert _read_hits(fimo_tsv, "baumannii") == []


def test_motif_pi0_estimates_are_associated_with_the_current_motif(tmp_path):
    log = tmp_path / "fimo.log"
    log.write_text(
        "Using motif +species_MEME-1 of width 6.\n"
        "Computing q-values.\n#   Estimated pi_0=0.912345\n"
        "Using motif -species_MEME-1 of width 6.\n"
        "Using motif +species_MEME-2 of width 9.\n"
        "#   Estimated pi_0=1\n",
        encoding="utf-8",
    )

    assert _motif_pi0_estimates(log, "species") == {"MEME-1": 0.912345, "MEME-2": 1.0}
