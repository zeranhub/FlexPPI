"""Residue correspondence for a partner-stripped chain and its complex."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Sequence, Tuple

from flexppi.data.structures import ProteinStructure


@dataclass
class Alignment:
    indices_a: List[int]
    indices_b: List[int]
    identity: float
    coverage: float
    score: int


def global_sequence_alignment(
    sequence_a: str,
    sequence_b: str,
    match_score: int = 2,
    mismatch_score: int = -1,
    gap_score: int = -2,
) -> Alignment:
    """Needleman-Wunsch alignment returning only exact residue correspondences."""

    n, m = len(sequence_a), len(sequence_b)
    scores = [[0] * (m + 1) for _ in range(n + 1)]
    trace = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        scores[i][0] = i * gap_score
        trace[i][0] = 1
    for j in range(1, m + 1):
        scores[0][j] = j * gap_score
        trace[0][j] = 2
    for i in range(1, n + 1):
        aa = sequence_a[i - 1]
        for j in range(1, m + 1):
            diagonal = scores[i - 1][j - 1] + (
                match_score if aa == sequence_b[j - 1] else mismatch_score
            )
            up = scores[i - 1][j] + gap_score
            left = scores[i][j - 1] + gap_score
            best = max(diagonal, up, left)
            scores[i][j] = best
            trace[i][j] = 0 if best == diagonal else (1 if best == up else 2)

    exact_a: List[int] = []
    exact_b: List[int] = []
    aligned_positions = 0
    exact_matches = 0
    i, j = n, m
    while i > 0 or j > 0:
        direction = trace[i][j]
        if i > 0 and j > 0 and direction == 0:
            aligned_positions += 1
            if sequence_a[i - 1] == sequence_b[j - 1]:
                exact_a.append(i - 1)
                exact_b.append(j - 1)
                exact_matches += 1
            i -= 1
            j -= 1
        elif i > 0 and (j == 0 or direction == 1):
            i -= 1
        else:
            j -= 1
    exact_a.reverse()
    exact_b.reverse()
    identity = exact_matches / max(aligned_positions, 1)
    coverage = exact_matches / max(n, 1)
    return Alignment(exact_a, exact_b, identity, coverage, scores[n][m])


def match_apo_to_holo(
    apo: ProteinStructure,
    holo: ProteinStructure,
) -> Tuple[List[int], List[int], List[str], float, float]:
    """Map every apo chain to the most similar unused chain in the complex."""

    apo_indices: List[int] = []
    holo_indices: List[int] = []
    matched_holo_chains: List[str] = []
    weighted_identity = 0.0
    weighted_coverage = 0.0
    total_weight = 0
    unused_holo = set(holo.unique_chains())

    for apo_chain in apo.unique_chains():
        source_global = apo.indices_for_chain(apo_chain).tolist()
        source_sequence = apo.chain_sequence(apo_chain)
        candidates = []
        for holo_chain in sorted(unused_holo):
            alignment = global_sequence_alignment(source_sequence, holo.chain_sequence(holo_chain))
            candidates.append((alignment.coverage, alignment.identity, alignment.score, holo_chain, alignment))
        if not candidates:
            continue
        _, _, _, best_chain, best = max(candidates, key=lambda value: value[:3])
        target_global = holo.indices_for_chain(best_chain).tolist()
        apo_indices.extend(source_global[i] for i in best.indices_a)
        holo_indices.extend(target_global[j] for j in best.indices_b)
        matched_holo_chains.append(best_chain)
        unused_holo.remove(best_chain)
        weight = len(source_sequence)
        weighted_identity += best.identity * weight
        weighted_coverage += best.coverage * weight
        total_weight += weight

    if not apo_indices:
        raise ValueError("No exact residue correspondence found between apo and holo structures")
    return (
        apo_indices,
        holo_indices,
        matched_holo_chains,
        weighted_identity / max(total_weight, 1),
        weighted_coverage / max(total_weight, 1),
    )

