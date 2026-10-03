from experiments.run_tjupan_reference_fimo import _read_reference_names
from experiments.summarize_tjupan_reference_overlap import (
    _read_fasta,
    _read_reference_hits,
    intervals_overlap,
)


def test_reference_matrix_contains_the_four_audited_templates():
    from pathlib import Path

    assert _read_reference_names(Path("baselines/known_promoter_elements.meme")) == {
        "Ecoli_sigma70_minus10",
        "Ecoli_sigma70_minus35",
        "Ecoli_UP_proximal",
        "Ecoli_UP_distal",
    }


def test_inclusive_interval_overlap_handles_shared_boundary():
    assert intervals_overlap({"start": 1, "stop": 6}, {"start": 6, "stop": 11})
    assert not intervals_overlap({"start": 1, "stop": 5}, {"start": 6, "stop": 11})


def test_reference_fimo_parser_filters_by_q_and_checks_template_id(tmp_path):
    path = tmp_path / "fimo.tsv"
    path.write_text(
        "# FIMO\n"
        "motif_id\tmotif_alt_id\tsequence_name\tstart\tstop\tstrand\tscore\tp-value\tq-value\tmatched_sequence\n"
        "Ecoli_sigma70_minus10\t.\tp1|label=1\t1\t6\t+\t8\t0.001\t0.01\tTATAAT\n"
        "Ecoli_sigma70_minus10\t.\tp2|label=0\t5\t10\t-\t8\t0.1\t0.2\tTATAAT\n",
        encoding="utf-8",
    )
    result = _read_reference_hits(path, {"Ecoli_sigma70_minus10"}, 0.05)
    assert len(result) == 1
    assert result[0]["sequence_name"] == "p1|label=1"


def test_fasta_parser_exports_only_label_and_sequence_hash(tmp_path):
    path = tmp_path / "input.fasta"
    sequence = "ACGT" * 20 + "A"
    path.write_text(f">p1|label=1\n{sequence}\n", encoding="ascii")
    records = _read_fasta(path)
    assert set(records["p1|label=1"]) == {"label", "sequence_sha256"}
    assert records["p1|label=1"]["label"] == "1"
