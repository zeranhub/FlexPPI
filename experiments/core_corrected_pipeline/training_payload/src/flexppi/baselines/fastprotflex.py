"""Python-compatible execution of the published FastProtFlex linear model.

The algorithm, coefficients, cutoffs, and preprocessing follow ``predict.r``
and ``FUNCTION_GDV.r`` from https://github.com/jure-praznikar/FastProtFlex.
FastProtFlex describes its distributed predictor as CC BY-NC 4.0 in the
README, so this adapter is intended for non-commercial research evaluation.
Orbit counts are computed with ``orca-graphlets``, which uses the same ORCA
node-orbit numbering 0--14 as the R implementation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import time
from typing import List, Tuple

import numpy as np

from flexppi.constants import THREE_TO_ONE


FASTPROTFLEX_BETA = np.asarray(
    [
        -1.9063749,
        -0.9075780,
        2.9751331,
        2.9623554,
        -0.9232815,
        -0.1610720,
        0.1171952,
        -0.4793347,
        0.2411250,
        0.3110215,
        0.4226776,
        -0.7755721,
        0.5147505,
        -1.7265063,
        -1.4770230,
    ],
    dtype=np.float64,
)


@dataclass
class HeavyAtomStructure:
    coordinates: np.ndarray
    atom_names: List[str]
    residue_keys: List[Tuple[str, str, str]]
    residue_names: List[str]
    ca_atom_indices: np.ndarray
    ca_residue_keys: List[Tuple[str, str, str]]


@dataclass
class FastProtFlexPrediction:
    atom_scores: np.ndarray
    ca_scores: np.ndarray
    ca_residue_keys: List[Tuple[str, str, str]]
    atom_count: int
    residue_count: int
    edge_count: int
    runtime_seconds: float


def _element(atom_name: str, element_column: str) -> str:
    explicit = element_column.strip().upper()
    if explicit:
        return explicit
    letters = "".join(character for character in atom_name if character.isalpha())
    return letters[:1].upper()


def parse_heavy_protein_atoms(path: str) -> HeavyAtomStructure:
    coordinates: List[np.ndarray] = []
    atom_names: List[str] = []
    atom_residue_keys: List[Tuple[str, str, str]] = []
    residue_order: List[Tuple[str, str, str]] = []
    residue_name_by_key = {}
    ca_by_key = {}
    saw_model = False
    with Path(path).open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            record = line[:6].strip()
            if record == "MODEL":
                if saw_model:
                    break
                saw_model = True
                continue
            if record == "ENDMDL":
                break
            if record not in {"ATOM", "HETATM"} or len(line) < 54:
                continue
            if line[16:17] not in {" ", "A", "1"}:
                continue
            residue_name = line[17:20].strip().upper()
            if residue_name not in THREE_TO_ONE:
                continue
            atom_name = line[12:16].strip().upper()
            if _element(atom_name, line[76:78] if len(line) >= 78 else "") in {"H", "D"}:
                continue
            key = (
                line[21:22].strip() or "_",
                line[22:26].strip(),
                line[26:27].strip(),
            )
            try:
                xyz = np.asarray(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=np.float64,
                )
            except ValueError:
                continue
            if key not in residue_name_by_key:
                residue_order.append(key)
                residue_name_by_key[key] = residue_name
            atom_index = len(coordinates)
            coordinates.append(xyz)
            atom_names.append(atom_name)
            atom_residue_keys.append(key)
            if atom_name == "CA" and key not in ca_by_key:
                ca_by_key[key] = atom_index
    if not coordinates:
        raise ValueError("No heavy protein atoms found in {}".format(path))
    ca_keys = [key for key in residue_order if key in ca_by_key]
    return HeavyAtomStructure(
        coordinates=np.stack(coordinates),
        atom_names=atom_names,
        residue_keys=atom_residue_keys,
        residue_names=[residue_name_by_key[key] for key in residue_order],
        ca_atom_indices=np.asarray([ca_by_key[key] for key in ca_keys], dtype=np.int64),
        ca_residue_keys=ca_keys,
    )


def _column_zscore(values: np.ndarray) -> np.ndarray:
    mean = np.mean(values, axis=0)
    standard_deviation = np.std(values, axis=0, ddof=1)
    if np.any(standard_deviation == 0):
        zero_columns = np.flatnonzero(standard_deviation == 0).tolist()
        raise ValueError("FastProtFlex orbit columns have zero variance: {}".format(zero_columns))
    return (values - mean) / standard_deviation


def predict_fastprotflex(path: str, maximum_atoms: int = 10000) -> FastProtFlexPrediction:
    """Predict per-atom normalized RMSF using the official all-at-once path."""

    try:
        from orca import orca_nodes
        from scipy.spatial import cKDTree
    except ImportError as exc:
        raise RuntimeError(
            "FastProtFlex benchmark requires: pip install -e '.[benchmark]'"
        ) from exc

    started = time.perf_counter()
    structure = parse_heavy_protein_atoms(path)
    atom_count = int(len(structure.coordinates))
    if atom_count > maximum_atoms:
        raise ValueError(
            "{} heavy atoms exceed the all-at-once limit {}; the published R code would use "
            "its by-parts path".format(atom_count, maximum_atoms)
        )
    tree = cKDTree(structure.coordinates)
    graph_edges = tree.query_pairs(7.0, output_type="ndarray").astype(np.int64, copy=False)
    degrees = np.bincount(graph_edges.reshape(-1), minlength=atom_count)
    if np.any(degrees == 0):
        raise ValueError("FastProtFlex graph contains isolated heavy atoms")
    orbit_counts = orca_nodes(graph_edges, num_nodes=atom_count, graphlet_size=4).astype(
        np.float64
    )

    covalent_edges = tree.query_pairs(2.1, output_type="ndarray").astype(np.int64, copy=False)
    neighbour_sum = np.zeros_like(orbit_counts)
    neighbour_count = np.zeros(atom_count, dtype=np.int64)
    if len(covalent_edges):
        left, right = covalent_edges[:, 0], covalent_edges[:, 1]
        np.add.at(neighbour_sum, left, orbit_counts[right])
        np.add.at(neighbour_sum, right, orbit_counts[left])
        np.add.at(neighbour_count, left, 1)
        np.add.at(neighbour_count, right, 1)
    smoothed = orbit_counts.copy()
    has_neighbours = neighbour_count > 0
    smoothed[has_neighbours] = 0.5 * (
        orbit_counts[has_neighbours]
        + neighbour_sum[has_neighbours] / neighbour_count[has_neighbours, None]
    )
    features = _column_zscore(np.log1p(smoothed))
    raw_prediction = features @ FASTPROTFLEX_BETA
    prediction_std = float(np.std(raw_prediction, ddof=1))
    if prediction_std == 0:
        raise ValueError("FastProtFlex produced a constant prediction")
    atom_scores = (raw_prediction - float(np.mean(raw_prediction))) / prediction_std
    return FastProtFlexPrediction(
        atom_scores=atom_scores,
        ca_scores=atom_scores[structure.ca_atom_indices],
        ca_residue_keys=structure.ca_residue_keys,
        atom_count=atom_count,
        residue_count=len(structure.ca_atom_indices),
        edge_count=int(len(graph_edges)),
        runtime_seconds=time.perf_counter() - started,
    )

