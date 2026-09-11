"""Small, dependency-light PWM implementation for known-element baselines."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable


BASES = "ACGT"
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")


def reverse_complement(sequence: str) -> str:
    return sequence.upper().translate(COMPLEMENT)[::-1]


@dataclass(frozen=True)
class PWMHit:
    motif_name: str
    sequence_id: str
    start: int
    end: int
    strand: str
    score: float
    p_value_proxy: float
    matched_sequence: str


@dataclass(frozen=True)
class PWM:
    name: str
    probabilities: tuple[dict[str, float], ...]
    background: dict[str, float]

    @property
    def width(self) -> int:
        return len(self.probabilities)

    def score(self, sequence: str) -> float:
        """Return the log2 likelihood ratio; N contributes zero information."""

        if len(sequence) != self.width:
            raise ValueError(f"expected width {self.width}, got {len(sequence)}")
        score = 0.0
        for position, base in enumerate(sequence.upper()):
            if base == "N":
                continue
            if base not in BASES:
                raise ValueError(f"invalid DNA base: {base}")
            score += math.log2(self.probabilities[position][base] / self.background[base])
        return score


def consensus_pwm(
    name: str,
    consensus: str,
    *,
    match_probability: float = 0.85,
    pseudocount: float = 0.01,
    background: dict[str, float] | None = None,
) -> PWM:
    """Create a PWM from an IUPAC-like consensus using a transparent prior.

    Supported symbols are A/C/G/T and N. Ambiguous IUPAC symbols are not
    silently guessed; they should be expanded into aligned training sequences
    before a PWM is created.
    """

    if not 0 < match_probability < 1:
        raise ValueError("match_probability must be in (0, 1)")
    background = background or {base: 0.25 for base in BASES}
    if set(background) != set(BASES) or not math.isclose(sum(background.values()), 1.0):
        raise ValueError("background must contain A/C/G/T probabilities summing to 1")
    probabilities = []
    for symbol in consensus.upper():
        if symbol not in set(BASES) | {"N"}:
            raise ValueError(f"unsupported consensus symbol: {symbol}")
        favored = set(BASES) if symbol == "N" else {symbol}
        other_probability = (1.0 - match_probability) / (len(BASES) - 1)
        row = {base: (match_probability if base in favored else other_probability) for base in BASES}
        row = {base: (value + pseudocount) for base, value in row.items()}
        normalizer = sum(row.values())
        probabilities.append({base: value / normalizer for base, value in row.items()})
    return PWM(name=name, probabilities=tuple(probabilities), background=background)


def scan_pwm(
    pwm: PWM,
    sequence_id: str,
    sequence: str,
    *,
    min_score: float = 2.0,
    both_strands: bool = True,
) -> list[PWMHit]:
    """Scan a sequence and return deterministic hits above ``min_score``."""

    sequence = sequence.upper()
    hits: list[PWMHit] = []
    strands = [("+", sequence)]
    if both_strands:
        strands.append(("-", reverse_complement(sequence)))
    for strand, oriented in strands:
        for start in range(0, len(oriented) - pwm.width + 1):
            window = oriented[start : start + pwm.width]
            if set(window) - set(BASES) - {"N"}:
                continue
            score = pwm.score(window)
            if score < min_score:
                continue
            hits.append(
                PWMHit(
                    motif_name=pwm.name,
                    sequence_id=sequence_id,
                    start=start,
                    end=start + pwm.width,
                    strand=strand,
                    score=score,
                    p_value_proxy=2 ** (-max(score, 0.0)),
                    matched_sequence=window,
                )
            )
    return sorted(hits, key=lambda hit: (hit.sequence_id, hit.start, hit.strand, -hit.score))


def scan_many(
    motifs: Iterable[PWM],
    records: Iterable[tuple[str, str]],
    *,
    min_score: float = 2.0,
    both_strands: bool = True,
) -> list[PWMHit]:
    hits = []
    for motif in motifs:
        for sequence_id, sequence in records:
            hits.extend(
                scan_pwm(
                    motif,
                    sequence_id,
                    sequence,
                    min_score=min_score,
                    both_strands=both_strands,
                )
            )
    return hits
