"""Dynamic label files and optional trajectory processing.

RMSF is computed after fitting each chain independently. This prevents rigid
motion of one partner relative to another from being mistaken for internal
residue flexibility. A single long trajectory is divided into fixed windows;
the variance across windows is retained as label uncertainty.
"""

from __future__ import annotations

import math
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from flexppi.constants import CHI_ATOMS, METRIC_NAMES, THREE_TO_ONE
from flexppi.data.records import StateRecord
from flexppi.data.structures import parse_pdb


def load_label_file(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    with np.load(path, allow_pickle=False) as data:
        values = np.asarray(data["values"], dtype=np.float32)
        variances = np.asarray(data["variances"], dtype=np.float32)
        mask = np.asarray(data["mask"], dtype=bool)
        names = [str(value) for value in data["metric_names"].tolist()]
    if values.ndim != 2 or values.shape != variances.shape or values.shape != mask.shape:
        raise ValueError("Malformed label tensors in {}".format(path))
    if names != list(METRIC_NAMES):
        reordered_values = np.full((values.shape[0], len(METRIC_NAMES)), np.nan, np.float32)
        reordered_variances = np.zeros_like(reordered_values)
        reordered_mask = np.zeros_like(reordered_values, dtype=bool)
        for source_index, name in enumerate(names):
            if name in METRIC_NAMES:
                target_index = METRIC_NAMES.index(name)
                reordered_values[:, target_index] = values[:, source_index]
                reordered_variances[:, target_index] = variances[:, source_index]
                reordered_mask[:, target_index] = mask[:, source_index]
        values, variances, mask = reordered_values, reordered_variances, reordered_mask
    return values, variances, mask


def save_label_file(
    path: str,
    values: np.ndarray,
    variances: np.ndarray,
    mask: np.ndarray,
) -> None:
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(
        output,
        values=np.asarray(values, dtype=np.float32),
        variances=np.asarray(variances, dtype=np.float32),
        mask=np.asarray(mask, dtype=bool),
        metric_names=np.asarray(METRIC_NAMES),
    )


def _kabsch_fit(mobile: np.ndarray, reference: np.ndarray) -> np.ndarray:
    mobile_center = mobile.mean(axis=0)
    reference_center = reference.mean(axis=0)
    centered_mobile = mobile - mobile_center
    centered_reference = reference - reference_center
    if mobile.shape[0] < 3:
        return centered_mobile + reference_center
    covariance = centered_mobile.T @ centered_reference
    u, _, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    if np.linalg.det(rotation) < 0:
        u[:, -1] *= -1
        rotation = u @ vt
    return centered_mobile @ rotation + reference_center


def _chainwise_rmsf(coordinates: np.ndarray, chain_index: np.ndarray) -> np.ndarray:
    aligned = np.empty_like(coordinates, dtype=np.float32)
    for chain in np.unique(chain_index):
        indices = np.flatnonzero(chain_index == chain)
        reference = coordinates[0, indices]
        for frame in range(coordinates.shape[0]):
            aligned[frame, indices] = _kabsch_fit(coordinates[frame, indices], reference)
    mean_coordinates = aligned.mean(axis=0)
    return np.sqrt(np.mean(np.sum((aligned - mean_coordinates) ** 2, axis=-1), axis=0))


def _dihedral_degrees(points: np.ndarray) -> np.ndarray:
    """Vectorised dihedral angles for points shaped [Q, 4, 3]."""

    b0 = -(points[:, 1] - points[:, 0])
    b1 = points[:, 2] - points[:, 1]
    b2 = points[:, 3] - points[:, 2]
    b1 = b1 / np.maximum(np.linalg.norm(b1, axis=-1, keepdims=True), 1e-8)
    v = b0 - np.sum(b0 * b1, axis=-1, keepdims=True) * b1
    w = b2 - np.sum(b2 * b1, axis=-1, keepdims=True) * b1
    x = np.sum(v * w, axis=-1)
    y = np.sum(np.cross(b1, v) * w, axis=-1)
    return np.degrees(np.arctan2(y, x))


def _entropy(codes: np.ndarray) -> float:
    _, counts = np.unique(codes, return_counts=True)
    probabilities = counts.astype(np.float64) / max(int(counts.sum()), 1)
    return float(-np.sum(probabilities * np.log(np.maximum(probabilities, 1e-12))))


def _sphere_points(count: int) -> np.ndarray:
    """Approximately uniform unit-sphere points using a Fibonacci spiral."""

    if count < 12:
        raise ValueError("SASA sphere point count must be at least 12")
    indices = np.arange(count, dtype=np.float64) + 0.5
    z = 1.0 - 2.0 * indices / count
    radius = np.sqrt(np.maximum(1.0 - z * z, 0.0))
    azimuth = math.pi * (3.0 - math.sqrt(5.0)) * indices
    return np.stack((radius * np.cos(azimuth), radius * np.sin(azimuth), z), axis=-1)


def _shrake_rupley_residue_sasa(
    positions: np.ndarray,
    expanded_radii: np.ndarray,
    atom_to_residue: np.ndarray,
    n_residues: int,
    sphere: np.ndarray,
) -> np.ndarray:
    """Compute residue SASA in Angstrom^2 for one frame.

    Each atom is surrounded by sampled points on its van-der-Waals radius plus
    a 1.4 Angstrom water probe. A point contributes area only when it is not
    inside a neighbouring atom's expanded sphere. A KD-tree limits checks to
    spatial neighbours instead of comparing every atom with every other atom.
    """

    from scipy.spatial import cKDTree

    coordinates = np.asarray(positions, dtype=np.float64)
    radii = np.asarray(expanded_radii, dtype=np.float64)
    result = np.zeros(n_residues, dtype=np.float64)
    if coordinates.shape[0] == 0:
        return result.astype(np.float32)
    tree = cKDTree(coordinates)
    maximum_radius = float(radii.max())
    for atom_index, (center, atom_radius) in enumerate(zip(coordinates, radii)):
        surface = center[None, :] + atom_radius * sphere
        neighbours = tree.query_ball_point(center, atom_radius + maximum_radius)
        neighbours = [index for index in neighbours if index != atom_index]
        accessible = np.ones(surface.shape[0], dtype=bool)
        if neighbours:
            neighbour_positions = coordinates[neighbours]
            neighbour_radii2 = radii[neighbours] ** 2
            distances2 = np.sum(
                (surface[:, None, :] - neighbour_positions[None, :, :]) ** 2,
                axis=-1,
            )
            accessible = ~np.any(distances2 < neighbour_radii2[None, :], axis=1)
        fraction = float(accessible.mean())
        area = 4.0 * math.pi * atom_radius * atom_radius * fraction
        result[int(atom_to_residue[atom_index])] += area
    return result.astype(np.float32)


def _window_ranges(length: int, window_frames: int) -> List[Tuple[int, int]]:
    if length < 2:
        return []
    if window_frames <= 0 or window_frames >= length:
        return [(0, length)]
    ranges = []
    for start in range(0, length - window_frames + 1, window_frames):
        ranges.append((start, start + window_frames))
    return ranges or [(0, length)]


def _summarise_samples(samples: Sequence[np.ndarray], n_residues: int) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    values = np.full((n_residues, len(METRIC_NAMES)), np.nan, dtype=np.float32)
    variances = np.zeros_like(values)
    mask = np.zeros_like(values, dtype=bool)
    if not samples:
        return values, variances, mask
    stacked = np.stack(samples, axis=0)
    for metric in range(stacked.shape[-1]):
        for residue in range(n_residues):
            finite = stacked[:, residue, metric]
            finite = finite[np.isfinite(finite)]
            if finite.size:
                values[residue, metric] = float(finite.mean())
                variances[residue, metric] = float(finite.var(ddof=1)) if finite.size > 1 else 0.0
                mask[residue, metric] = True
    return values, variances, mask


def compute_dynarepo_labels(
    record: StateRecord,
    output_path: str,
    window_frames: int = 1000,
    frame_stride: int = 1,
    include_sasa: bool = False,
    max_windows: Optional[int] = None,
    sasa_points: int = 96,
) -> Dict[str, int]:
    """Compute RMSF, chi entropy, and optionally SASA from local trajectories.

    MDAnalysis is intentionally an optional dependency because the core model
    and ATLAS pretraining can be tested without it. SASA uses a local
    Shrake-Rupley implementation and is returned directly in Angstrom^2.
    """

    try:
        import MDAnalysis as mda
    except ImportError as exc:
        raise RuntimeError(
            "Trajectory labels require the optional dependency: pip install -e .[md]"
        ) from exc

    structure = parse_pdb(record.structure_path)
    n_residues = structure.num_residues
    samples: List[np.ndarray] = []
    consumed_windows = 0

    for trajectory_path in record.trajectory_paths:
        universe = mda.Universe(record.structure_path, trajectory_path)
        md_residues = [
            residue for residue in universe.residues
            if str(residue.resname).upper() in THREE_TO_ONE
            and any(str(atom.name).upper() == "CA" for atom in residue.atoms)
        ]
        if len(md_residues) != n_residues:
            raise ValueError(
                "PDB/trajectory residue mismatch for {}: {} versus {}".format(
                    record.sample_id, n_residues, len(md_residues)
                )
            )
        ca_indices = []
        chi_quadruples: List[List[int]] = []
        residue_to_chis: Dict[int, List[int]] = {}
        for residue_index, residue in enumerate(md_residues):
            atom_lookup = {str(atom.name).upper(): int(atom.index) for atom in residue.atoms}
            ca_indices.append(atom_lookup["CA"])
            definitions = CHI_ATOMS.get(str(residue.resname).upper(), [])
            local_chis = []
            for definition in definitions:
                if all(name in atom_lookup for name in definition):
                    local_chis.append(len(chi_quadruples))
                    chi_quadruples.append([atom_lookup[name] for name in definition])
            if local_chis and len(local_chis) == len(definitions):
                residue_to_chis[residue_index] = local_chis

        sasa_atom_indices: List[int] = []
        sasa_atom_radii: List[float] = []
        sasa_atom_residues: List[int] = []
        sphere = _sphere_points(sasa_points) if include_sasa else None
        if include_sasa:
            van_der_waals = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "P": 1.80}
            for residue_index, residue in enumerate(md_residues):
                for atom in residue.atoms:
                    name = str(atom.name).strip().upper()
                    element = next((character for character in name if character.isalpha()), "C")
                    if element == "H":
                        continue
                    sasa_atom_indices.append(int(atom.index))
                    sasa_atom_radii.append(van_der_waals.get(element, 1.70) + 1.40)
                    sasa_atom_residues.append(residue_index)

        ca_frames: List[np.ndarray] = []
        chi_state_frames: List[np.ndarray] = []
        sasa_frames: List[np.ndarray] = []
        quad_array = np.asarray(chi_quadruples, dtype=np.int64)
        for timestep in universe.trajectory[::frame_stride]:
            positions = timestep.positions
            ca_frames.append(np.asarray(positions[ca_indices], dtype=np.float32).copy())
            if len(quad_array):
                angles = _dihedral_degrees(positions[quad_array])
                bins = np.floor((angles + 180.0) / 120.0).astype(np.int8)
                chi_state_frames.append(np.clip(bins, 0, 2))
            if include_sasa and sphere is not None:
                sasa_frames.append(
                    _shrake_rupley_residue_sasa(
                        positions[sasa_atom_indices],
                        np.asarray(sasa_atom_radii, dtype=np.float32),
                        np.asarray(sasa_atom_residues, dtype=np.int64),
                        n_residues,
                        sphere,
                    )
                )
            if max_windows is not None:
                remaining_windows = max_windows - consumed_windows
                if len(ca_frames) >= max(remaining_windows, 1) * window_frames:
                    break

        ca_array = np.stack(ca_frames, axis=0)
        chi_array = np.stack(chi_state_frames, axis=0) if chi_state_frames else None
        sasa_array = np.stack(sasa_frames, axis=0) if sasa_frames else None

        for start, stop in _window_ranges(len(ca_array), window_frames):
            sample = np.full((n_residues, len(METRIC_NAMES)), np.nan, dtype=np.float32)
            sample[:, 0] = _chainwise_rmsf(ca_array[start:stop], structure.chain_index)
            if chi_array is not None:
                for residue_index, chi_indices in residue_to_chis.items():
                    codes = np.zeros(stop - start, dtype=np.int32)
                    for power, chi_index in enumerate(chi_indices):
                        codes += chi_array[start:stop, chi_index].astype(np.int32) * (3 ** power)
                    sample[residue_index, 1] = _entropy(codes)
            if sasa_array is not None:
                sample[:, 2] = sasa_array[start:stop].mean(axis=0)
            samples.append(sample)
            consumed_windows += 1
            if max_windows is not None and consumed_windows >= max_windows:
                break
        if max_windows is not None and consumed_windows >= max_windows:
            break

    values, variances, mask = _summarise_samples(samples, n_residues)
    save_label_file(output_path, values, variances, mask)
    return {
        "residues": n_residues,
        "windows": consumed_windows,
        "rmsf_labels": int(mask[:, 0].sum()),
        "chi_labels": int(mask[:, 1].sum()),
        "sasa_labels": int(mask[:, 2].sum()),
    }
