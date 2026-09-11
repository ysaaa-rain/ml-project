from models.pwm import consensus_pwm, reverse_complement, scan_pwm


def test_consensus_pwm_prefers_consensus_sequence():
    pwm = consensus_pwm("minus_10", "TATAAT")
    assert pwm.score("TATAAT") > pwm.score("CCCCCC")


def test_reverse_complement_and_both_strand_scan():
    pwm = consensus_pwm("minus_35", "TTGACA")
    sequence = "AAAA" + reverse_complement("TTGACA") + "CCC"
    hits = scan_pwm(pwm, "seq1", sequence, min_score=2.0, both_strands=True)
    assert hits
    target_hits = [hit for hit in hits if hit.strand == "-" and hit.start == 4]
    assert target_hits
    assert {hit.end for hit in target_hits} == {10}
