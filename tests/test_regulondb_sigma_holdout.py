from experiments.validate_regulondb_sigma_holdout import (
    ROOT,
    _require_project_output,
    _split_paired_records,
    _write_fimo_matrix,
)


def _example_records():
    primary = [(f"p{i}", "A" * 81) for i in range(10)]
    controls = [(f"p{i}__dinucleotide1", "C" * 81) for i in range(10)]
    return primary, controls


def test_split_is_reproducible_and_keeps_matched_controls_together():
    primary, controls = _example_records()
    first = _split_paired_records(primary, controls, fraction=0.2, seed=31)
    second = _split_paired_records(primary, controls, fraction=0.2, seed=31)
    assert first == second
    assert len(first["holdout"]["primary"]) == 2
    for partition in ("train", "holdout"):
        primary_ids = {identifier for identifier, _ in first[partition]["primary"]}
        control_ids = {identifier.removesuffix("__dinucleotide1") for identifier, _ in first[partition]["control"]}
        assert primary_ids == control_ids


def test_split_partitions_are_disjoint_and_recover_all_records():
    primary, controls = _example_records()
    split = _split_paired_records(primary, controls, fraction=0.3, seed=9)
    train_ids = {identifier for identifier, _ in split["train"]["primary"]}
    holdout_ids = {identifier for identifier, _ in split["holdout"]["primary"]}
    assert train_ids.isdisjoint(holdout_ids)
    assert train_ids | holdout_ids == {identifier for identifier, _ in primary}


def test_split_rejects_unpaired_control_ids():
    primary, controls = _example_records()
    with_error = controls[:-1] + [("other__dinucleotide1", "C" * 81)]
    try:
        _split_paired_records(primary, with_error)
    except ValueError as error:
        assert "not paired" in str(error)
    else:
        raise AssertionError("unpaired controls must be rejected")


def test_fimo_matrix_headers_use_xml_motif_aliases(tmp_path):
    source = tmp_path / "meme.txt"
    target = tmp_path / "normalized.meme"
    source.write_text(
        f"PRIMARY SEQUENCES= {ROOT}/tmp/train.fasta\n"
        "MOTIF ACGT MEME-1 width = 4\nletter-probability matrix\n"
        "MOTIF TGCA MEME-2 width = 4\nletter-probability matrix\n",
        encoding="utf-8",
    )
    _write_fimo_matrix(source, target, ["MEME-1", "MEME-2"])
    headers = [line for line in target.read_text(encoding="utf-8").splitlines() if line.startswith("MOTIF ")]
    assert headers == ["MOTIF 1 MEME-1", "MOTIF 2 MEME-2"]
    normalized = target.read_text(encoding="utf-8")
    assert str(ROOT) not in normalized
    assert all(line == line.rstrip() for line in normalized.splitlines())


def test_experiment_outputs_must_stay_inside_project():
    _require_project_output(ROOT / "tmp" / "example")
    try:
        _require_project_output(ROOT.parent / "outside-project")
    except ValueError as error:
        assert "inside the project directory" in str(error)
    else:
        raise AssertionError("outputs outside the project must be rejected")
