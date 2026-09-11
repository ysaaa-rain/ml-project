import hashlib

from experiments.preflight import MEME_TOOLS, collect_environment
from preprocessing.provenance import sha256_file


def test_sha256_file_matches_standard_library(tmp_path):
    payload = b"PR01-02 provenance test\n"
    target = tmp_path / "input.tsv"
    target.write_bytes(payload)
    assert sha256_file(target) == hashlib.sha256(payload).hexdigest()


def test_environment_snapshot_has_required_shape():
    snapshot = collect_environment()
    assert snapshot["python"]["version"]
    assert set(snapshot["meme_suite"]["tools"]) == set(MEME_TOOLS)
    assert isinstance(snapshot["meme_suite"]["ready"], bool)
