from experiments.analyze_tjupan_motif_pairs import _quartiles, _read_fasta, site_relation


def test_site_relation_reports_order_and_uncovered_gap():
    assert site_relation({"start": 1, "stop": 6}, {"start": 10, "stop": 15}) == ("A_before_B", 3)
    assert site_relation({"start": 10, "stop": 15}, {"start": 1, "stop": 6}) == ("B_before_A", 3)


def test_site_relation_encodes_overlapping_base_count_as_negative_gap():
    assert site_relation({"start": 1, "stop": 8}, {"start": 6, "stop": 12}) == ("overlapping", -3)


def test_spacing_quartiles_handle_empty_single_and_multiple_values():
    assert _quartiles([]) == (None, None)
    assert _quartiles([4]) == (4.0, 4.0)
    assert _quartiles([1, 2, 3, 4]) == (1.75, 3.25)


def test_fasta_reader_retains_hash_not_raw_sequence(tmp_path):
    fasta = tmp_path / "records.fasta"
    fasta.write_text(f">one|label=1\n{'A' * 81}\n", encoding="ascii")
    record = _read_fasta(fasta)["one|label=1"]
    assert set(record) == {"label", "sequence_sha256"}
    assert record["label"] == "1"
