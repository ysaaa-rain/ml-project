from pathlib import Path


ALPHABET = "ACGT"
REFERENCE = Path(__file__).parents[1] / "baselines" / "known_promoter_elements.meme"


def _consensus_by_unique_max(matrix):
    consensus = []
    for row in matrix:
        maximum = max(row)
        winners = [letter for letter, probability in zip(ALPHABET, row) if probability == maximum]
        consensus.append(winners[0] if len(winners) == 1 else "N")
    return "".join(consensus)


def _read_motifs(path):
    motifs = {}
    lines = path.read_text(encoding="utf-8").splitlines()
    index = 0
    while index < len(lines):
        if not lines[index].startswith("MOTIF "):
            index += 1
            continue
        name = lines[index].split()[1]
        index += 1
        while index < len(lines) and not lines[index].startswith("letter-probability matrix:"):
            index += 1
        if index == len(lines):
            raise AssertionError(f"missing matrix for {name}")
        width = int(lines[index].split("w=", 1)[1].split()[0])
        index += 1
        matrix = []
        for _ in range(width):
            matrix.append([float(value) for value in lines[index].split()])
            index += 1
        motifs[name] = matrix
    return motifs


def test_sigma70_reference_matrices_encode_the_named_consensus_boxes():
    motifs = _read_motifs(REFERENCE)

    assert _consensus_by_unique_max(motifs["Ecoli_sigma70_minus10"]) == "TATAAT"
    assert _consensus_by_unique_max(motifs["Ecoli_sigma70_minus35"]) == "TTGACA"
