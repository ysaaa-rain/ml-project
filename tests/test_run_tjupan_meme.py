from experiments.run_tjupan_meme import _attempt_history


def test_attempt_history_is_empty_for_first_run():
    assert _attempt_history(None) == []


def test_attempt_history_preserves_interrupted_attempts_on_retry():
    prior = {
        "status": "interrupted",
        "partial_motifs_written": 3,
        "attempt_history": [{"status": "failed", "return_code": 2}],
    }

    history = _attempt_history(prior)

    assert history == [
        {"status": "failed", "return_code": 2},
        prior,
    ]
    assert history[1] is not prior
