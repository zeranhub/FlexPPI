"""Extract and audit dynamic PPI-interface descriptors from DynaRepo pairs.

This script deliberately does not import or modify the FlexPPI model.  Its job
is to answer a prior scientific question: do trajectory-derived interface
descriptors contain reproducible information about apo-to-bound dynamical
redistribution beyond the static distance-to-partner gradient?

The analysis unit is a PDB complex (not a residue).  Every trajectory is split
into non-overlapping 100 ns blocks.  Dynamic descriptors are accepted as useful
only when they are reproducible between blocks and improve leave-one-PDB-out
explanation of the observed RMSF response beyond static geometry.
"""

from __future__ import annotations

import argparse
import json
import math
import traceback
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.cluster.hierarchy import fcluster, linkage
from scipy.spatial.distance import pdist
from scipy.stats import rankdata, spearmanr

from flexppi.constants import THREE_TO_ONE
from flexppi.data.labels import (
    _chainwise_rmsf,
    _shrake_rupley_residue_sasa,
    _sphere_points,
    load_label_file,
)
from flexppi.data.structures import parse_pdb


MAX_ASA = {
    "ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0,
    "CYS": 167.0, "GLN": 225.0, "GLU": 223.0, "GLY": 104.0,
    "HIS": 224.0, "HID": 224.0, "HIE": 224.0, "HIP": 224.0,
    "ILE": 197.0, "LEU": 201.0, "LYS": 236.0, "MET": 224.0,
    "MSE": 224.0, "PHE": 240.0, "PRO": 159.0, "SER": 155.0,
    "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0,
}

DYNAMIC_RESIDUE_FEATURES = [
    "interface_persistence_7a",
    "interface_persistence_8a",
    "interface_persistence_9a",
    "distance_mean",
    "distance_std",
    "contact_count_mean",
    "contact_count_std",
    "edge_entropy_mean",
    "edge_entropy_sum",
    "switching_edge_fraction",
    "stable_edge_fraction",
    "partner_identity_entropy",
    "initial_contact_retention",
]

STATIC_FEATURES = [
    "static_distance_exp",
    "static_burial_proxy",
    "static_partner_count",
    "static_local_ca_contacts",
]

GLOBAL_DYNAMIC_FEATURES = [
    "state_entropy_025_norm",
    "state_entropy_035_norm",
    "state_entropy_045_norm",
    "contact_turnover",
    "consensus_dispersion",
    "total_contacts_cv",
    "pose_translation_rms",
    "pose_rotation_rms_deg",
    "pose_rmsd_rms",
]


def _read_jsonl(path: Path) -> List[dict]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def _local_state_paths(record: Mapping[str, object], root: Path, label_dir: Path) -> Tuple[Path, Path, Path]:
    sample_id = str(record["sample_id"])
    sample_root = root / sample_id
    structure = sample_root / "structure.pdb"
    trajectory = sample_root / "trajectory.xtc"
    label = label_dir / (sample_id + ".npz")
    for path in (structure, trajectory, label):
        if not path.exists():
            raise FileNotFoundError(path)
    return structure, trajectory, label


def _standard_residues(universe: object) -> Tuple[List[object], List[str]]:
    residues = []
    chains = []
    for residue in universe.residues:
        name = str(residue.resname).upper()
        atom_names = {str(atom.name).upper() for atom in residue.atoms}
        if name not in THREE_TO_ONE or "CA" not in atom_names:
            continue
        residues.append(residue)
        chain = str(residue.atoms[0].chainID).strip()
        if not chain:
            chain = str(residue.segid).strip() or "_"
        chains.append(chain)
    return residues, chains


def _atom_lookup(residue: object) -> Dict[str, int]:
    return {str(atom.name).upper(): int(atom.index) for atom in residue.atoms}


def _is_heavy_atom(atom: object) -> bool:
    name = str(atom.name).strip().upper()
    return bool(name) and not name.startswith("H")


def _kabsch_transform(mobile: np.ndarray, reference: np.ndarray) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    centered_mobile = mobile - mobile_center
    centered_reference = reference - reference_center
    covariance = centered_mobile.T @ centered_reference
    u, _, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    aligned = centered_mobile @ rotation + reference_center
    translation = reference_center - mobile_center @ rotation
    return aligned, rotation, translation


def _rotation_angle_degrees(rotation: np.ndarray) -> float:
    cosine = np.clip((np.trace(rotation) - 1.0) / 2.0, -1.0, 1.0)
    return float(np.degrees(np.arccos(cosine)))


def _entropy(probabilities: np.ndarray) -> float:
    p = np.asarray(probabilities, dtype=np.float64)
    p = p[np.isfinite(p) & (p > 0)]
    if not len(p):
        return 0.0
    p = p / p.sum()
    return float(-np.sum(p * np.log(p)))


def _binary_entropy(probabilities: np.ndarray) -> np.ndarray:
    p = np.clip(np.asarray(probabilities, dtype=np.float64), 1e-12, 1.0 - 1e-12)
    result = -(p * np.log(p) + (1.0 - p) * np.log(1.0 - p))
    result[(probabilities <= 0) | (probabilities >= 1)] = 0.0
    return result


def _jaccard_turnover(first: set[int], second: set[int]) -> float:
    union = first | second
    if not union:
        return 0.0
    return 1.0 - len(first & second) / len(union)


def _contact_state_statistics(contact_sets: Sequence[set[int]]) -> Dict[str, float]:
    n_frames = len(contact_sets)
    if n_frames == 0:
        return {name: float("nan") for name in (
            "contact_turnover", "consensus_dispersion", "total_contacts_mean",
            "total_contacts_std", "total_contacts_cv", "unique_edges",
            "state_entropy_025", "state_entropy_025_norm", "state_clusters_025",
            "state_entropy_035", "state_entropy_035_norm", "state_clusters_035",
            "state_entropy_045", "state_entropy_045_norm", "state_clusters_045",
        )}
    all_edges = sorted(set().union(*contact_sets)) if contact_sets else []
    sizes = np.asarray([len(value) for value in contact_sets], dtype=np.float64)
    turnover = np.asarray(
        [_jaccard_turnover(contact_sets[i - 1], contact_sets[i]) for i in range(1, n_frames)],
        dtype=np.float64,
    )
    counts: Dict[int, int] = {}
    for values in contact_sets:
        for edge in values:
            counts[edge] = counts.get(edge, 0) + 1
    consensus = {edge for edge, count in counts.items() if count / n_frames >= 0.5}
    dispersion = np.asarray([_jaccard_turnover(value, consensus) for value in contact_sets])
    result: Dict[str, float] = {
        "contact_turnover": float(np.mean(turnover)) if len(turnover) else 0.0,
        "consensus_dispersion": float(np.mean(dispersion)),
        "total_contacts_mean": float(np.mean(sizes)),
        "total_contacts_std": float(np.std(sizes, ddof=1)) if n_frames > 1 else 0.0,
        "total_contacts_cv": float(np.std(sizes, ddof=1) / max(np.mean(sizes), 1e-8)) if n_frames > 1 else 0.0,
        "unique_edges": float(len(all_edges)),
    }
    if len(all_edges) == 0 or n_frames < 3:
        for threshold in (0.25, 0.35, 0.45):
            suffix = "%03d" % int(round(threshold * 100))
            result["state_entropy_" + suffix] = 0.0
            result["state_entropy_" + suffix + "_norm"] = 0.0
            result["state_clusters_" + suffix] = 1.0
        return result
    edge_to_column = {edge: index for index, edge in enumerate(all_edges)}
    matrix = np.zeros((n_frames, len(all_edges)), dtype=bool)
    for frame, values in enumerate(contact_sets):
        for edge in values:
            matrix[frame, edge_to_column[edge]] = True
    distances = pdist(matrix, metric="jaccard")
    distances = np.nan_to_num(distances, nan=0.0, posinf=1.0, neginf=0.0)
    hierarchy = linkage(distances, method="average")
    for threshold in (0.25, 0.35, 0.45):
        suffix = "%03d" % int(round(threshold * 100))
        labels = fcluster(hierarchy, t=threshold, criterion="distance")
        _, cluster_counts = np.unique(labels, return_counts=True)
        state_entropy = _entropy(cluster_counts.astype(np.float64))
        result["state_entropy_" + suffix] = state_entropy
        result["state_entropy_" + suffix + "_norm"] = state_entropy / max(math.log(n_frames), 1e-8)
        result["state_clusters_" + suffix] = float(len(cluster_counts))
    return result


def _residue_contact_statistics(
    contact_sets: Sequence[set[int]],
    min_distances: np.ndarray,
    contact_counts: np.ndarray,
    n_target: int,
    n_partner: int,
) -> pd.DataFrame:
    n_frames = len(contact_sets)
    edge_counts: Dict[int, int] = {}
    for values in contact_sets:
        for edge in values:
            edge_counts[edge] = edge_counts.get(edge, 0) + 1
    edges_by_target: List[List[Tuple[int, float]]] = [[] for _ in range(n_target)]
    for edge, count in edge_counts.items():
        target = edge // n_partner
        partner = edge % n_partner
        edges_by_target[target].append((partner, count / max(n_frames, 1)))
    initial = contact_sets[0] if contact_sets else set()
    rows = []
    for target in range(n_target):
        entries = edges_by_target[target]
        occupancies = np.asarray([value for _, value in entries], dtype=np.float64)
        entropies = _binary_entropy(occupancies) if len(occupancies) else np.asarray([], dtype=np.float64)
        initial_occupancies = [
            edge_counts.get(edge, 0) / max(n_frames, 1)
            for edge in initial if edge // n_partner == target
        ]
        rows.append({
            "target_local_index": target,
            "interface_persistence": float(np.mean(contact_counts[:, target] > 0)),
            "interface_persistence_7a": float(np.mean(min_distances[:, target] <= 7.0)),
            "interface_persistence_8a": float(np.mean(min_distances[:, target] <= 8.0)),
            "interface_persistence_9a": float(np.mean(min_distances[:, target] <= 9.0)),
            "distance_mean": float(np.mean(min_distances[:, target])),
            "distance_std": float(np.std(min_distances[:, target], ddof=1)) if n_frames > 1 else 0.0,
            "contact_count_mean": float(np.mean(contact_counts[:, target])),
            "contact_count_std": float(np.std(contact_counts[:, target], ddof=1)) if n_frames > 1 else 0.0,
            "unique_partner_count": float(len(entries)),
            "edge_entropy_mean": float(np.mean(entropies)) if len(entropies) else 0.0,
            "edge_entropy_sum": float(np.sum(entropies)) if len(entropies) else 0.0,
            "edge_entropy_max": float(np.max(entropies)) if len(entropies) else 0.0,
            "switching_edge_fraction": float(np.mean((occupancies >= 0.2) & (occupancies <= 0.8))) if len(occupancies) else 0.0,
            "stable_edge_fraction": float(np.mean(occupancies >= 0.8)) if len(occupancies) else 0.0,
            "partner_identity_entropy": _entropy(occupancies),
            "initial_contact_retention": float(np.mean(initial_occupancies)) if initial_occupancies else 0.0,
        })
    return pd.DataFrame(rows)


def _static_sasa_features(
    universe: object,
    standard_residues: Sequence[object],
    target_global: np.ndarray,
    sphere_points: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    sphere = _sphere_points(sphere_points)
    atom_indices = []
    radii = []
    atom_to_residue = []
    for residue_index, residue in enumerate(standard_residues):
        for atom in residue.atoms:
            if not _is_heavy_atom(atom):
                continue
            name = str(atom.name).strip().upper()
            element = next((character for character in name if character.isalpha()), "C")
            radius = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}.get(element, 1.70)
            atom_indices.append(int(atom.index))
            radii.append(radius + 1.40)
            atom_to_residue.append(residue_index)
    positions = universe.atoms.positions[np.asarray(atom_indices, dtype=np.int64)]
    radii_array = np.asarray(radii, dtype=np.float32)
    residue_array = np.asarray(atom_to_residue, dtype=np.int64)
    complex_sasa = _shrake_rupley_residue_sasa(
        positions, radii_array, residue_array, len(standard_residues), sphere
    )
    global_to_target = {int(value): index for index, value in enumerate(target_global)}
    keep = np.asarray([int(value) in global_to_target for value in residue_array], dtype=bool)
    target_positions = positions[keep]
    target_radii = radii_array[keep]
    target_atom_residues = np.asarray(
        [global_to_target[int(value)] for value in residue_array[keep]], dtype=np.int64
    )
    alone_sasa = _shrake_rupley_residue_sasa(
        target_positions, target_radii, target_atom_residues, len(target_global), sphere
    )
    complex_target = complex_sasa[target_global]
    delta_sasa = alone_sasa - complex_target
    maximum = np.asarray(
        [MAX_ASA.get(str(standard_residues[int(index)].resname).upper(), 220.0) for index in target_global],
        dtype=np.float32,
    )
    rasa_alone = alone_sasa / maximum
    rasa_complex = complex_target / maximum
    return delta_sasa, rasa_alone, rasa_complex, complex_target


def _load_all_ca_frames(
    structure_path: Path,
    trajectory_path: Path,
    frame_count: int,
) -> Tuple[np.ndarray, np.ndarray]:
    import MDAnalysis as mda

    universe = mda.Universe(str(structure_path), str(trajectory_path))
    residues, chains = _standard_residues(universe)
    ca_indices = [_atom_lookup(residue)["CA"] for residue in residues]
    chain_to_index: Dict[str, int] = {}
    chain_index = np.zeros(len(chains), dtype=np.int64)
    for index, chain in enumerate(chains):
        if chain not in chain_to_index:
            chain_to_index[chain] = len(chain_to_index)
        chain_index[index] = chain_to_index[chain]
    frames = np.empty((frame_count, len(ca_indices), 3), dtype=np.float32)
    for frame, timestep in enumerate(universe.trajectory[:frame_count]):
        frames[frame] = timestep.positions[ca_indices]
    if frame + 1 != frame_count:
        raise ValueError("Expected {} frames, loaded {}".format(frame_count, frame + 1))
    return frames, chain_index


def _extract_pair(
    pair: Mapping[str, object],
    root: Path,
    label_dir: Path,
    frame_count: int,
    block_frames: int,
    contact_stride: int,
    contact_cutoff: float,
    distance_cap: float,
    sasa_points: int,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    import MDAnalysis as mda
    from MDAnalysis.lib.distances import capped_distance

    apo = pair["apo"]
    holo = pair["holo"]
    apo_structure, apo_trajectory, apo_label = _local_state_paths(apo, root, label_dir)
    holo_structure, holo_trajectory, holo_label = _local_state_paths(holo, root, label_dir)
    target_global = np.asarray(pair["holo_indices"], dtype=np.int64)
    apo_global = np.asarray(pair["apo_indices"], dtype=np.int64)
    if len(target_global) != len(apo_global):
        raise ValueError("Apo/holo mapped residue counts differ")
    universe = mda.Universe(str(holo_structure), str(holo_trajectory))
    standard_residues, standard_chains = _standard_residues(universe)
    parsed = parse_pdb(str(holo_structure))
    if len(standard_residues) != parsed.num_residues:
        raise ValueError("MDAnalysis/parser residue mismatch: {} versus {}".format(
            len(standard_residues), parsed.num_residues
        ))
    if standard_chains != parsed.chain_ids:
        raise ValueError("MDAnalysis/parser chain ordering mismatch")
    all_global = np.arange(len(standard_residues), dtype=np.int64)
    target_set = {int(value) for value in target_global}
    partner_global = np.asarray([value for value in all_global if int(value) not in target_set], dtype=np.int64)
    if len(partner_global) == 0:
        raise ValueError("No partner residues remain after target selection")

    target_ca_indices = np.asarray(
        [_atom_lookup(standard_residues[int(index)])["CA"] for index in target_global], dtype=np.int64
    )
    all_ca_indices = np.asarray(
        [_atom_lookup(residue)["CA"] for residue in standard_residues], dtype=np.int64
    )
    partner_ca_indices = np.asarray(
        [_atom_lookup(standard_residues[int(index)])["CA"] for index in partner_global], dtype=np.int64
    )
    n_target = len(target_global)
    n_partner = len(partner_global)
    sample_frames = list(range(0, frame_count, contact_stride))
    contact_sets: List[set[int]] = []
    min_distances = np.full((len(sample_frames), n_target), distance_cap, dtype=np.float32)
    contact_counts = np.zeros((len(sample_frames), n_target), dtype=np.float32)
    holo_ca_all = np.empty((frame_count, len(standard_residues), 3), dtype=np.float32)
    pose_translation = []
    pose_rotation = []
    pose_rmsd = []
    target_fit_rmsd = []
    static_burial_proxy = np.zeros(n_target, dtype=np.float32)
    reference_target = None
    reference_partner = None
    sample_lookup = {frame: index for index, frame in enumerate(sample_frames)}
    for frame, timestep in enumerate(universe.trajectory[:frame_count]):
        target_ca = timestep.positions[target_ca_indices].copy()
        holo_ca_all[frame] = timestep.positions[all_ca_indices]
        if frame not in sample_lookup:
            continue
        sample_index = sample_lookup[frame]
        partner_ca = timestep.positions[partner_ca_indices].copy()
        pairs, distances = capped_distance(
            target_ca,
            partner_ca,
            max_cutoff=distance_cap,
            box=timestep.dimensions,
            return_distances=True,
            method="nsgrid",
        )
        if len(pairs):
            target_local = pairs[:, 0]
            np.minimum.at(min_distances[sample_index], target_local, distances.astype(np.float32))
            if sample_index == 0:
                burial_pairs = distances <= 12.0
                if np.any(burial_pairs):
                    static_burial_proxy = np.bincount(
                        pairs[burial_pairs, 0], minlength=n_target
                    ).astype(np.float32)
            close = distances <= contact_cutoff
            if np.any(close):
                edge_codes = np.unique(
                    pairs[close, 0] * n_partner + pairs[close, 1]
                )
            else:
                edge_codes = np.asarray([], dtype=np.int64)
        else:
            edge_codes = np.asarray([], dtype=np.int64)
        values = {int(value) for value in edge_codes.tolist()}
        contact_sets.append(values)
        if values:
            targets = np.fromiter((value // n_partner for value in values), dtype=np.int64)
            contact_counts[sample_index] = np.bincount(targets, minlength=n_target)
        if reference_target is None:
            reference_target = target_ca.copy()
            reference_partner = partner_ca.copy()
        aligned_target, target_rotation, target_translation = _kabsch_transform(
            target_ca, reference_target
        )
        aligned_partner = partner_ca @ target_rotation + target_translation
        target_fit_rmsd.append(float(np.sqrt(np.mean(np.sum((aligned_target - reference_target) ** 2, axis=1)))))
        pose_translation.append(float(np.linalg.norm(aligned_partner.mean(axis=0) - reference_partner.mean(axis=0))))
        centered_aligned = aligned_partner - aligned_partner.mean(axis=0)
        centered_reference = reference_partner - reference_partner.mean(axis=0)
        _, relative_rotation, _ = _kabsch_transform(centered_aligned, centered_reference)
        pose_rotation.append(_rotation_angle_degrees(relative_rotation))
        pose_rmsd.append(float(np.sqrt(np.mean(np.sum((aligned_partner - reference_partner) ** 2, axis=1)))))

    if len(contact_sets) != len(sample_frames):
        raise ValueError("Contact frame count mismatch")
    apo_ca_all, apo_chain_index = _load_all_ca_frames(
        apo_structure, apo_trajectory, frame_count
    )
    n_blocks = frame_count // block_frames
    if n_blocks < 2 or n_blocks * block_frames != frame_count:
        raise ValueError("frame_count must be an exact multiple of block_frames")
    block_delta = []
    block_apo_rmsf = []
    block_holo_rmsf = []
    for block in range(n_blocks):
        start = block * block_frames
        stop = start + block_frames
        apo_rmsf = _chainwise_rmsf(apo_ca_all[start:stop], apo_chain_index)[apo_global]
        holo_rmsf = _chainwise_rmsf(
            holo_ca_all[start:stop], parsed.chain_index
        )[target_global]
        block_apo_rmsf.append(apo_rmsf)
        block_holo_rmsf.append(holo_rmsf)
        block_delta.append(holo_rmsf - apo_rmsf)
    block_delta_array = np.stack(block_delta)
    block_apo_array = np.stack(block_apo_rmsf)
    block_holo_array = np.stack(block_holo_rmsf)

    universe.trajectory[0]
    if sasa_points > 0:
        static_delta_sasa, rasa_alone, rasa_complex, complex_sasa = _static_sasa_features(
            universe, standard_residues, target_global, sasa_points
        )
    else:
        static_delta_sasa = np.full(n_target, np.nan, dtype=np.float32)
        rasa_alone = np.full(n_target, np.nan, dtype=np.float32)
        rasa_complex = np.full(n_target, np.nan, dtype=np.float32)
        complex_sasa = np.full(n_target, np.nan, dtype=np.float32)
    target_ca0 = universe.atoms.positions[target_ca_indices]
    ca_distances = np.sqrt(np.sum((target_ca0[:, None, :] - target_ca0[None, :, :]) ** 2, axis=-1))
    local_ca_contacts = np.sum((ca_distances < 10.0) & (ca_distances > 1e-6), axis=1)

    apo_values, apo_variances, apo_mask = load_label_file(str(apo_label))
    holo_values, holo_variances, holo_mask = load_label_file(str(holo_label))
    delta_rmsf_label = holo_values[target_global, 0] - apo_values[apo_global, 0]
    delta_rmsf_variance = holo_variances[target_global, 0] + apo_variances[apo_global, 0]
    delta_chi = holo_values[target_global, 1] - apo_values[apo_global, 1]
    delta_chi_mask = holo_mask[target_global, 1] & apo_mask[apo_global, 1]

    full_residue = _residue_contact_statistics(
        contact_sets, min_distances, contact_counts, n_target, n_partner
    )
    full_global = _contact_state_statistics(contact_sets)
    full_global.update({
        "pose_translation_mean": float(np.mean(pose_translation)),
        "pose_translation_rms": float(np.sqrt(np.mean(np.square(pose_translation)))),
        "pose_translation_std": float(np.std(pose_translation, ddof=1)),
        "pose_rotation_mean_deg": float(np.mean(pose_rotation)),
        "pose_rotation_rms_deg": float(np.sqrt(np.mean(np.square(pose_rotation)))),
        "pose_rotation_std_deg": float(np.std(pose_rotation, ddof=1)),
        "pose_rmsd_mean": float(np.mean(pose_rmsd)),
        "pose_rmsd_rms": float(np.sqrt(np.mean(np.square(pose_rmsd)))),
        "target_fit_rmsd_mean": float(np.mean(target_fit_rmsd)),
    })
    static_distance = min_distances[0].copy()
    full_residue["pair_id"] = str(pair["pair_id"])
    full_residue["pdb_id"] = str(pair["metadata"]["pdb_id"])
    full_residue["apo_sample_id"] = str(apo["sample_id"])
    full_residue["holo_sample_id"] = str(holo["sample_id"])
    full_residue["target_global_index"] = target_global
    full_residue["chain_id"] = [parsed.chain_ids[int(index)] for index in target_global]
    full_residue["author_residue_id"] = [parsed.author_residue_ids[int(index)] for index in target_global]
    full_residue["residue_name"] = [parsed.residue_names[int(index)] for index in target_global]
    full_residue["static_distance"] = static_distance
    full_residue["static_distance_exp"] = -np.exp(-static_distance / 5.0)
    full_residue["static_partner_count"] = contact_counts[0]
    full_residue["static_burial_proxy"] = static_burial_proxy
    full_residue["static_delta_sasa"] = static_delta_sasa
    full_residue["static_rasa_alone"] = rasa_alone
    full_residue["static_rasa_complex"] = rasa_complex
    full_residue["static_complex_sasa"] = complex_sasa
    full_residue["static_local_ca_contacts"] = local_ca_contacts
    if sasa_points > 0:
        full_residue["interface_region"] = np.where(
            (rasa_alone < 0.25) & (rasa_complex < 0.25), "support",
            np.where((rasa_alone >= 0.25) & (rasa_complex < 0.25), "core",
                     np.where((rasa_alone >= 0.25) & (rasa_complex >= 0.25), "rim", "other")),
        )
    else:
        full_residue["interface_region"] = "not_computed"
    full_residue["delta_rmsf"] = delta_rmsf_label
    full_residue["delta_rmsf_variance"] = delta_rmsf_variance
    full_residue["delta_rmsf_z"] = delta_rmsf_label / np.sqrt(np.maximum(delta_rmsf_variance, 1e-6))
    full_residue["delta_chi_entropy"] = np.where(delta_chi_mask, delta_chi, np.nan)
    full_residue["delta_rmsf_block_mean"] = block_delta_array.mean(axis=0)
    full_residue["delta_rmsf_block_variance"] = block_delta_array.var(axis=0, ddof=1)

    block_residue_rows = []
    block_complex_rows = []
    sampled_per_block = block_frames // contact_stride
    for block in range(n_blocks):
        sample_start = block * sampled_per_block
        sample_stop = sample_start + sampled_per_block
        block_sets = contact_sets[sample_start:sample_stop]
        block_distances = min_distances[sample_start:sample_stop]
        block_counts = contact_counts[sample_start:sample_stop]
        table = _residue_contact_statistics(
            block_sets, block_distances, block_counts, n_target, n_partner
        )
        table["pair_id"] = str(pair["pair_id"])
        table["pdb_id"] = str(pair["metadata"]["pdb_id"])
        table["block"] = block
        table["delta_rmsf_block"] = block_delta_array[block]
        table["apo_rmsf_block"] = block_apo_array[block]
        table["holo_rmsf_block"] = block_holo_array[block]
        table["static_distance"] = static_distance
        table["static_delta_sasa"] = static_delta_sasa
        table["static_burial_proxy"] = static_burial_proxy
        block_residue_rows.append(table)
        global_row = _contact_state_statistics(block_sets)
        global_row.update({
            "pair_id": str(pair["pair_id"]),
            "pdb_id": str(pair["metadata"]["pdb_id"]),
            "block": block,
            "pose_translation_rms": float(np.sqrt(np.mean(np.square(pose_translation[sample_start:sample_stop])))),
            "pose_rotation_rms_deg": float(np.sqrt(np.mean(np.square(pose_rotation[sample_start:sample_stop])))),
            "pose_rmsd_rms": float(np.sqrt(np.mean(np.square(pose_rmsd[sample_start:sample_stop])))),
            "mean_abs_delta_rmsf_block": float(np.mean(np.abs(block_delta_array[block]))),
            "delta_rmsf_sd_block": float(np.std(block_delta_array[block], ddof=1)),
        })
        block_complex_rows.append(global_row)

    full_global.update({
        "pair_id": str(pair["pair_id"]),
        "pdb_id": str(pair["metadata"]["pdb_id"]),
        "apo_sample_id": str(apo["sample_id"]),
        "holo_sample_id": str(holo["sample_id"]),
        "target_residues": n_target,
        "partner_residues": n_partner,
        "mean_abs_delta_rmsf": float(np.mean(np.abs(delta_rmsf_label))),
        "delta_rmsf_sd": float(np.std(delta_rmsf_label, ddof=1)),
        "interface_fraction": float(np.mean(static_distance <= contact_cutoff)),
        "dynamic_interface_fraction": float(np.mean(full_residue["interface_persistence"] >= 0.5)),
        "mean_interface_persistence": float(np.mean(full_residue["interface_persistence"])),
    })
    qc = {
        "pair_id": str(pair["pair_id"]),
        "pdb_id": str(pair["metadata"]["pdb_id"]),
        "target_residues": n_target,
        "partner_residues": n_partner,
        "contact_frames": len(contact_sets),
        "raw_frames": frame_count,
        "label_vs_recomputed_delta_spearman": _safe_spearman(delta_rmsf_label, block_delta_array.mean(axis=0)),
        "label_vs_recomputed_delta_mae": float(np.mean(np.abs(delta_rmsf_label - block_delta_array.mean(axis=0)))),
        "target_fit_rmsd_mean": float(np.mean(target_fit_rmsd)),
        "pose_translation_max": float(np.max(pose_translation)),
        "pose_rotation_max_deg": float(np.max(pose_rotation)),
    }
    return (
        full_residue,
        pd.concat(block_residue_rows, ignore_index=True),
        pd.DataFrame([full_global]),
        pd.DataFrame(block_complex_rows),
        qc,
    )


def _safe_spearman(first: Sequence[float], second: Sequence[float]) -> float:
    x = np.asarray(first, dtype=np.float64)
    y = np.asarray(second, dtype=np.float64)
    valid = np.isfinite(x) & np.isfinite(y)
    if np.sum(valid) < 3 or np.std(x[valid]) <= 1e-12 or np.std(y[valid]) <= 1e-12:
        return float("nan")
    return float(spearmanr(x[valid], y[valid]).statistic)


def _partial_spearman(x: np.ndarray, y: np.ndarray, controls: np.ndarray) -> float:
    valid = np.isfinite(x) & np.isfinite(y) & np.all(np.isfinite(controls), axis=1)
    if np.sum(valid) < max(8, controls.shape[1] + 4):
        return float("nan")
    xr = rankdata(x[valid])
    yr = rankdata(y[valid])
    cr = np.column_stack([rankdata(controls[valid, column]) for column in range(controls.shape[1])])
    design = np.column_stack((np.ones(len(cr)), cr))
    x_residual = xr - design @ np.linalg.lstsq(design, xr, rcond=None)[0]
    y_residual = yr - design @ np.linalg.lstsq(design, yr, rcond=None)[0]
    return _safe_spearman(x_residual, y_residual)


def _icc_oneway(matrix: np.ndarray) -> float:
    values = np.asarray(matrix, dtype=np.float64)
    valid_rows = np.all(np.isfinite(values), axis=1)
    values = values[valid_rows]
    n, k = values.shape if values.ndim == 2 else (0, 0)
    if n < 3 or k < 2:
        return float("nan")
    grand = values.mean()
    row_means = values.mean(axis=1)
    between = k * np.sum((row_means - grand) ** 2) / max(n - 1, 1)
    within = np.sum((values - row_means[:, None]) ** 2) / max(n * (k - 1), 1)
    denominator = between + (k - 1) * within
    return float((between - within) / denominator) if denominator > 1e-12 else float("nan")


def _bootstrap_complex_mean(values: pd.DataFrame, column: str, repetitions: int = 2000, seed: int = 7) -> Tuple[float, float, float]:
    grouped = values.groupby("pdb_id", sort=False)[column].mean().dropna()
    array = grouped.to_numpy(dtype=np.float64)
    if len(array) == 0:
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    samples = np.asarray([
        np.mean(rng.choice(array, size=len(array), replace=True)) for _ in range(repetitions)
    ])
    return float(np.mean(array)), float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def _feature_reliability(block_residue: pd.DataFrame) -> pd.DataFrame:
    rows = []
    feature_names = DYNAMIC_RESIDUE_FEATURES + ["delta_rmsf_block"]
    for feature in feature_names:
        pair_rows = []
        for (pdb_id, pair_id), group in block_residue.groupby(["pdb_id", "pair_id"], sort=False):
            pivot = group.pivot(index="target_local_index", columns="block", values=feature)
            if pivot.shape[1] < 2:
                continue
            correlations = []
            columns = list(pivot.columns)
            for i, first in enumerate(columns):
                for second in columns[i + 1:]:
                    correlations.append(_safe_spearman(pivot[first], pivot[second]))
            pair_rows.append({
                "pdb_id": pdb_id,
                "pair_id": pair_id,
                "block_spearman": float(np.nanmedian(correlations)),
                "icc": _icc_oneway(pivot.to_numpy()),
            })
        pair_table = pd.DataFrame(pair_rows)
        for metric in ("block_spearman", "icc"):
            mean, low, high = _bootstrap_complex_mean(pair_table, metric)
            rows.append({
                "feature": feature,
                "metric": metric,
                "macro_mean": mean,
                "ci_low": low,
                "ci_high": high,
                "evaluable_pairs": len(pair_table),
            })
    return pd.DataFrame(rows)


def _global_feature_reliability(block_complex: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feature in GLOBAL_DYNAMIC_FEATURES:
        pivot = block_complex.pivot(index="pair_id", columns="block", values=feature)
        pair_to_pdb = block_complex.drop_duplicates("pair_id").set_index("pair_id")["pdb_id"]
        pivot = pivot.loc[pivot.index.intersection(pair_to_pdb.index)]
        block_correlations = []
        columns = list(pivot.columns)
        for i, first in enumerate(columns):
            for second in columns[i + 1:]:
                block_correlations.append(_safe_spearman(pivot[first], pivot[second]))
        rows.append({
            "feature": feature,
            "block_spearman_median": float(np.nanmedian(block_correlations)),
            "block_spearman_min": float(np.nanmin(block_correlations)),
            "icc_across_complexes": _icc_oneway(pivot.to_numpy()),
            "pairs": len(pivot),
        })
    return pd.DataFrame(rows)


def _cross_half_associations(block_residue: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feature in DYNAMIC_RESIDUE_FEATURES:
        for outcome_name in ("delta_rmsf", "abs_delta_rmsf"):
            pair_rows = []
            for (pdb_id, pair_id), group in block_residue.groupby(["pdb_id", "pair_id"], sort=False):
                early = group[group["block"].isin([0, 1])].groupby("target_local_index", sort=True).mean(numeric_only=True)
                late = group[group["block"].isin([3, 4])].groupby("target_local_index", sort=True).mean(numeric_only=True)
                common = early.index.intersection(late.index)
                if len(common) < 8:
                    continue
                x = early.loc[common, feature].to_numpy(dtype=np.float64)
                y = late.loc[common, "delta_rmsf_block"].to_numpy(dtype=np.float64)
                if outcome_name == "abs_delta_rmsf":
                    y = np.abs(y)
                controls = early.loc[common, ["static_distance", "static_burial_proxy"]].to_numpy(dtype=np.float64)
                pair_rows.append({
                    "pdb_id": pdb_id,
                    "pair_id": pair_id,
                    "spearman": _safe_spearman(x, y),
                    "partial_spearman": _partial_spearman(x, y, controls),
                })
            pair_table = pd.DataFrame(pair_rows)
            for metric in ("spearman", "partial_spearman"):
                mean, low, high = _bootstrap_complex_mean(pair_table, metric)
                rows.append({
                    "feature": feature,
                    "outcome": outcome_name,
                    "metric": metric,
                    "macro_mean": mean,
                    "ci_low": low,
                    "ci_high": high,
                    "evaluable_pairs": len(pair_table),
                })
    return pd.DataFrame(rows)


def _weighted_ridge_fit(
    train: pd.DataFrame,
    test: pd.DataFrame,
    feature_names: Sequence[str],
    alpha: float,
) -> np.ndarray:
    x_train = train[list(feature_names)].to_numpy(dtype=np.float64)
    y_train = train["delta_rmsf"].to_numpy(dtype=np.float64)
    x_test = test[list(feature_names)].to_numpy(dtype=np.float64)
    group_sizes = train.groupby("pdb_id")["pdb_id"].transform("size").to_numpy(dtype=np.float64)
    weights = 1.0 / np.maximum(group_sizes, 1.0)
    weights = weights / weights.sum()
    mean_x = np.sum(x_train * weights[:, None], axis=0)
    var_x = np.sum((x_train - mean_x) ** 2 * weights[:, None], axis=0)
    scale_x = np.sqrt(np.maximum(var_x, 1e-8))
    mean_y = float(np.sum(y_train * weights))
    x_centered = (x_train - mean_x) / scale_x
    test_centered = (x_test - mean_x) / scale_x
    y_centered = y_train - mean_y
    sqrt_weight = np.sqrt(weights * len(weights))
    design = x_centered * sqrt_weight[:, None]
    response = y_centered * sqrt_weight
    penalty = alpha * np.eye(len(feature_names))
    coefficients = np.linalg.solve(design.T @ design + penalty, design.T @ response)
    return mean_y + test_centered @ coefficients


def _loco_models(residue: pd.DataFrame, alphas: Sequence[float]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    feature_sets = {
        "distance_only": ["static_distance_exp"],
        "static": STATIC_FEATURES,
        "dynamic": DYNAMIC_RESIDUE_FEATURES,
        "static_dynamic": STATIC_FEATURES + DYNAMIC_RESIDUE_FEATURES,
    }
    required = sorted(set(sum(feature_sets.values(), [])) | {"delta_rmsf"})
    data = residue[np.all(np.isfinite(residue[required].to_numpy(dtype=np.float64)), axis=1)].copy()
    predictions = []
    for alpha in alphas:
        for feature_set, names in feature_sets.items():
            for pdb_id in data["pdb_id"].drop_duplicates():
                train = data[data["pdb_id"] != pdb_id]
                test = data[data["pdb_id"] == pdb_id]
                prediction = _weighted_ridge_fit(train, test, names, alpha)
                table = test[["pdb_id", "pair_id", "target_local_index", "delta_rmsf"]].copy()
                table["feature_set"] = feature_set
                table["alpha"] = alpha
                table["prediction"] = prediction
                predictions.append(table)
    prediction_table = pd.concat(predictions, ignore_index=True)
    summary_rows = []
    for (feature_set, alpha), group in prediction_table.groupby(["feature_set", "alpha"], sort=False):
        per_pdb = []
        for pdb_id, pdb_group in group.groupby("pdb_id", sort=False):
            per_pdb.append({
                "pdb_id": pdb_id,
                "spearman": _safe_spearman(pdb_group["prediction"], pdb_group["delta_rmsf"]),
                "mae": float(np.mean(np.abs(pdb_group["prediction"] - pdb_group["delta_rmsf"]))),
            })
        per_pdb_table = pd.DataFrame(per_pdb)
        row = {"feature_set": feature_set, "alpha": alpha, "pdb_count": len(per_pdb_table)}
        for metric in ("spearman", "mae"):
            mean, low, high = _bootstrap_complex_mean(per_pdb_table, metric)
            row[metric] = mean
            row[metric + "_ci_low"] = low
            row[metric + "_ci_high"] = high
        summary_rows.append(row)
    return prediction_table, pd.DataFrame(summary_rows)


def _loco_paired_comparisons(predictions: pd.DataFrame) -> pd.DataFrame:
    metric_rows = []
    for (alpha, feature_set, pdb_id), group in predictions.groupby(
        ["alpha", "feature_set", "pdb_id"], sort=False
    ):
        metric_rows.append({
            "alpha": alpha,
            "feature_set": feature_set,
            "pdb_id": pdb_id,
            "spearman": _safe_spearman(group["prediction"], group["delta_rmsf"]),
            "mae": float(np.mean(np.abs(group["prediction"] - group["delta_rmsf"]))),
        })
    metric_table = pd.DataFrame(metric_rows)
    comparisons = [
        ("dynamic", "distance_only"),
        ("static_dynamic", "static"),
        ("static_dynamic", "distance_only"),
    ]
    rows = []
    for alpha in sorted(metric_table["alpha"].unique()):
        alpha_table = metric_table[metric_table["alpha"] == alpha]
        for candidate, reference in comparisons:
            candidate_table = alpha_table[alpha_table["feature_set"] == candidate].set_index("pdb_id")
            reference_table = alpha_table[alpha_table["feature_set"] == reference].set_index("pdb_id")
            common = candidate_table.index.intersection(reference_table.index)
            paired = pd.DataFrame({
                "pdb_id": common,
                "spearman_gain": candidate_table.loc[common, "spearman"].to_numpy()
                - reference_table.loc[common, "spearman"].to_numpy(),
                "mae_reduction": reference_table.loc[common, "mae"].to_numpy()
                - candidate_table.loc[common, "mae"].to_numpy(),
            })
            row = {
                "alpha": alpha,
                "candidate": candidate,
                "reference": reference,
                "pdb_count": len(common),
            }
            for metric in ("spearman_gain", "mae_reduction"):
                mean, low, high = _bootstrap_complex_mean(paired, metric)
                row[metric] = mean
                row[metric + "_ci_low"] = low
                row[metric + "_ci_high"] = high
            rows.append(row)
    return pd.DataFrame(rows)


def _bootstrap_spearman_rows(
    table: pd.DataFrame,
    feature: str,
    outcome: str,
    repetitions: int = 2000,
    seed: int = 7,
) -> Tuple[float, float, float, int]:
    grouped = table.groupby("pdb_id", as_index=False)[[feature, outcome]].mean(numeric_only=True)
    valid = np.isfinite(grouped[feature]) & np.isfinite(grouped[outcome])
    grouped = grouped[valid]
    observed = _safe_spearman(grouped[feature], grouped[outcome])
    if len(grouped) < 4:
        return observed, float("nan"), float("nan"), len(grouped)
    rng = np.random.default_rng(seed)
    values = []
    for _ in range(repetitions):
        indices = rng.integers(0, len(grouped), size=len(grouped))
        sampled = grouped.iloc[indices]
        value = _safe_spearman(sampled[feature], sampled[outcome])
        if np.isfinite(value):
            values.append(value)
    return observed, float(np.quantile(values, 0.025)), float(np.quantile(values, 0.975)), len(grouped)


def _global_associations(complex_table: pd.DataFrame, residue: pd.DataFrame) -> pd.DataFrame:
    outcomes = residue.groupby(["pdb_id", "pair_id"], as_index=False).agg(
        mean_abs_delta_rmsf=("delta_rmsf", lambda values: float(np.mean(np.abs(values)))),
        delta_rmsf_sd=("delta_rmsf", "std"),
        strong_response_fraction=("delta_rmsf_z", lambda values: float(np.mean(np.abs(values) >= 1.0))),
    )
    table = complex_table.merge(outcomes, on=["pdb_id", "pair_id"], how="left", suffixes=("", "_outcome"))
    rows = []
    for feature in GLOBAL_DYNAMIC_FEATURES:
        for outcome in ("mean_abs_delta_rmsf", "delta_rmsf_sd", "strong_response_fraction"):
            value, low, high, count = _bootstrap_spearman_rows(table, feature, outcome)
            rows.append({
                "feature": feature,
                "outcome": outcome,
                "spearman": value,
                "ci_low": low,
                "ci_high": high,
                "pdbs": count,
            })
    return pd.DataFrame(rows)


def _feature_decisions(reliability: pd.DataFrame, associations: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for feature in DYNAMIC_RESIDUE_FEATURES:
        block = reliability[(reliability["feature"] == feature) & (reliability["metric"] == "block_spearman")]
        icc = reliability[(reliability["feature"] == feature) & (reliability["metric"] == "icc")]
        assoc = associations[
            (associations["feature"] == feature)
            & (associations["outcome"] == "delta_rmsf")
            & (associations["metric"] == "partial_spearman")
        ]
        block_value = float(block["macro_mean"].iloc[0]) if len(block) else float("nan")
        icc_value = float(icc["macro_mean"].iloc[0]) if len(icc) else float("nan")
        assoc_value = float(assoc["macro_mean"].iloc[0]) if len(assoc) else float("nan")
        assoc_low = float(assoc["ci_low"].iloc[0]) if len(assoc) else float("nan")
        assoc_high = float(assoc["ci_high"].iloc[0]) if len(assoc) else float("nan")
        reproducible = np.isfinite(block_value) and np.isfinite(icc_value) and block_value >= 0.5 and icc_value >= 0.5
        associated = np.isfinite(assoc_value) and abs(assoc_value) >= 0.08 and not (assoc_low <= 0 <= assoc_high)
        if reproducible and associated:
            decision = "retain_for_followup"
        elif reproducible:
            decision = "stable_but_not_response_specific"
        else:
            decision = "too_noisy_or_definition_sensitive"
        rows.append({
            "feature": feature,
            "block_spearman": block_value,
            "icc": icc_value,
            "cross_half_partial_spearman": assoc_value,
            "association_ci_low": assoc_low,
            "association_ci_high": assoc_high,
            "decision": decision,
        })
    return pd.DataFrame(rows)


def _plot_audit(
    reliability: pd.DataFrame,
    associations: pd.DataFrame,
    loco: pd.DataFrame,
    global_associations: pd.DataFrame,
    output_png: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.0, 9.5))
    rel = reliability[reliability["metric"] == "block_spearman"].set_index("feature")
    assoc = associations[
        (associations["outcome"] == "delta_rmsf")
        & (associations["metric"] == "partial_spearman")
    ].set_index("feature")
    common = [feature for feature in DYNAMIC_RESIDUE_FEATURES if feature in rel.index and feature in assoc.index]
    axes[0, 0].scatter(rel.loc[common, "macro_mean"], assoc.loc[common, "macro_mean"], s=45)
    for feature in common:
        axes[0, 0].annotate(feature.replace("_", " "), (rel.loc[feature, "macro_mean"], assoc.loc[feature, "macro_mean"]), fontsize=7)
    axes[0, 0].axvline(0.5, color="0.5", linestyle="--", linewidth=1)
    axes[0, 0].axhline(0, color="0.5", linewidth=1)
    axes[0, 0].set_xlabel("100-ns block reproducibility (macro Spearman)")
    axes[0, 0].set_ylabel("early-feature vs late-ΔRMSF partial Spearman")
    axes[0, 0].set_title("A. Reproducibility versus response specificity")

    selected_alpha = min(loco["alpha"].unique(), key=lambda value: abs(math.log10(value) - 1.0))
    subset = loco[loco["alpha"] == selected_alpha].copy()
    x = np.arange(len(subset))
    axes[0, 1].bar(x - 0.18, subset["spearman"], width=0.36, label="Spearman")
    axes[0, 1].bar(x + 0.18, -subset["mae"], width=0.36, label="-MAE")
    axes[0, 1].set_xticks(x, subset["feature_set"], rotation=25, ha="right")
    axes[0, 1].axhline(0, color="black", linewidth=0.8)
    axes[0, 1].set_title("B. Leave-one-PDB-out models (ridge α={})".format(selected_alpha))
    axes[0, 1].legend(frameon=False)

    global_sorted = global_associations.sort_values("spearman", key=lambda values: np.abs(values), ascending=False).head(10)
    axes[1, 0].barh(np.arange(len(global_sorted)), global_sorted["spearman"], color=np.where(global_sorted["spearman"] >= 0, "#4477AA", "#CC6677"))
    axes[1, 0].set_yticks(np.arange(len(global_sorted)), ["{} → {}".format(a.replace("_", " "), b.replace("_", " ")) for a, b in zip(global_sorted["feature"], global_sorted["outcome"])], fontsize=7)
    axes[1, 0].invert_yaxis()
    axes[1, 0].axvline(0, color="black", linewidth=0.8)
    axes[1, 0].set_xlabel("Across-complex Spearman")
    axes[1, 0].set_title("C. Global interface dynamics versus response heterogeneity")

    dyn = loco[loco["feature_set"].isin(["static", "static_dynamic"])].pivot(index="alpha", columns="feature_set", values=["spearman", "mae"])
    axes[1, 1].plot(dyn.index, dyn[("spearman", "static_dynamic")] - dyn[("spearman", "static")], marker="o", label="Δ Spearman")
    axes[1, 1].plot(dyn.index, dyn[("mae", "static")] - dyn[("mae", "static_dynamic")], marker="s", label="MAE reduction")
    axes[1, 1].set_xscale("log")
    axes[1, 1].axhline(0, color="black", linewidth=0.8)
    axes[1, 1].set_xlabel("Ridge α")
    axes[1, 1].set_ylabel("Dynamic-feature gain over static model")
    axes[1, 1].set_title("D. Sensitivity to regularization")
    axes[1, 1].legend(frameon=False)
    fig.tight_layout()
    fig.savefig(output_png, dpi=220, bbox_inches="tight")
    fig.savefig(output_png.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pair-manifest", default="artifacts/data/dynarepo_pairs.jsonl")
    parser.add_argument("--root", default=r"D:\data\dynarepo_non_nucleic_full")
    parser.add_argument("--label-dir", default="artifacts/data/dynarepo_labels")
    parser.add_argument("--output-dir", default="artifacts/dynarepo_interface_dynamics")
    parser.add_argument("--frame-count", type=int, default=5000)
    parser.add_argument("--block-frames", type=int, default=1000)
    parser.add_argument("--contact-stride", type=int, default=5)
    parser.add_argument(
        "--contact-cutoff", type=float, default=8.0,
        help="C-alpha contact cutoff in Angstrom; 7/8/9-A persistence is always reported",
    )
    parser.add_argument("--distance-cap", type=float, default=20.0)
    parser.add_argument(
        "--sasa-points", type=int, default=0,
        help="Shrake-Rupley sphere points; 0 skips costly SASA and uses a fast atom-occlusion proxy",
    )
    parser.add_argument("--max-pairs", type=int, default=None)
    parser.add_argument("--only-pdb", nargs="*", default=None)
    parser.add_argument("--skip-analysis", action="store_true")
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    pairs = _read_jsonl(Path(args.pair_manifest))
    if args.only_pdb:
        wanted = {value.upper() for value in args.only_pdb}
        pairs = [pair for pair in pairs if str(pair["metadata"]["pdb_id"]).upper() in wanted]
    if args.max_pairs is not None:
        pairs = pairs[:args.max_pairs]
    residue_tables = []
    block_residue_tables = []
    complex_tables = []
    block_complex_tables = []
    qc_rows = []
    failures = []
    for index, pair in enumerate(pairs, 1):
        pair_id = str(pair["pair_id"])
        print("[{}/{}] {}".format(index, len(pairs), pair_id), flush=True)
        try:
            residue, block_residue, complex_table, block_complex, qc = _extract_pair(
                pair,
                Path(args.root),
                Path(args.label_dir),
                args.frame_count,
                args.block_frames,
                args.contact_stride,
                args.contact_cutoff,
                args.distance_cap,
                args.sasa_points,
            )
            residue_tables.append(residue)
            block_residue_tables.append(block_residue)
            complex_tables.append(complex_table)
            block_complex_tables.append(block_complex)
            qc_rows.append(qc)
        except Exception as exc:
            failures.append({
                "pair_id": pair_id,
                "error": repr(exc),
                "traceback": traceback.format_exc(),
            })
            print("FAILED {}: {}".format(pair_id, exc), flush=True)
    if not residue_tables:
        raise RuntimeError("No pairs were extracted successfully")
    residue = pd.concat(residue_tables, ignore_index=True)
    block_residue = pd.concat(block_residue_tables, ignore_index=True)
    complex_table = pd.concat(complex_tables, ignore_index=True)
    block_complex = pd.concat(block_complex_tables, ignore_index=True)
    qc_table = pd.DataFrame(qc_rows)
    residue.to_csv(output_dir / "residue_features.csv", index=False)
    block_residue.to_csv(output_dir / "block_residue_features.csv", index=False)
    complex_table.to_csv(output_dir / "complex_features.csv", index=False)
    block_complex.to_csv(output_dir / "block_complex_features.csv", index=False)
    qc_table.to_csv(output_dir / "trajectory_qc.csv", index=False)
    with (output_dir / "failures.json").open("w", encoding="utf-8") as handle:
        json.dump(failures, handle, indent=2, ensure_ascii=False)
    if args.skip_analysis or len(complex_table["pdb_id"].unique()) < 3:
        return
    reliability = _feature_reliability(block_residue)
    global_reliability = _global_feature_reliability(block_complex)
    associations = _cross_half_associations(block_residue)
    predictions, loco = _loco_models(residue, alphas=[0.1, 1.0, 10.0, 100.0])
    loco_comparisons = _loco_paired_comparisons(predictions)
    global_associations = _global_associations(complex_table, residue)
    decisions = _feature_decisions(reliability, associations)
    reliability.to_csv(output_dir / "feature_reliability.csv", index=False)
    global_reliability.to_csv(output_dir / "global_feature_reliability.csv", index=False)
    associations.to_csv(output_dir / "cross_half_associations.csv", index=False)
    predictions.to_csv(output_dir / "loco_predictions.csv", index=False)
    loco.to_csv(output_dir / "loco_summary.csv", index=False)
    loco_comparisons.to_csv(output_dir / "loco_paired_comparisons.csv", index=False)
    global_associations.to_csv(output_dir / "global_associations.csv", index=False)
    decisions.to_csv(output_dir / "feature_decisions.csv", index=False)
    _plot_audit(
        reliability,
        associations,
        loco,
        global_associations,
        output_dir / "interface_dynamics_audit.png",
    )
    report = {
        "pairs_requested": len(pairs),
        "pairs_succeeded": len(complex_table),
        "pdbs_succeeded": int(complex_table["pdb_id"].nunique()),
        "residues": len(residue),
        "failures": len(failures),
        "settings": vars(args),
        "feature_decisions": decisions.to_dict(orient="records"),
        "loco_summary": loco.to_dict(orient="records"),
        "loco_paired_comparisons": loco_comparisons.to_dict(orient="records"),
    }
    with (output_dir / "analysis_report.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, ensure_ascii=False)


if __name__ == "__main__":
    main()
