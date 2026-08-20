"""Run the official RMSF-net PDB-only checkpoint on a frozen manifest.

The upstream command requires UCSF Chimera even for ``mode 1`` because it uses
``molmap`` to turn a PDB structure into a 4-Angstrom density grid.  Chimera is
not redistributable in a Docker image, so this adapter reproduces the documented
``molmap`` operation directly: atomic-number-weighted Gaussians, sigma equal to
``resolution / (pi * sqrt(2))``, a 1.5-Angstrom grid, five-sigma truncation,
and 3*resolution edge padding.  The network architecture and PDB-only weights
are unmodified files from the pinned official RMSF-net repository.

RMSF-net's released prediction command writes a within-structure normalized
PDB.  This batch file preserves the pre-normalization positive score.  The
evaluation therefore treats it as a relative-flexibility method and compares
rank, within-system standardized error, and flexible-residue retrieval rather
than assigning an absolute-Angstrom error to it.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path
import time
from typing import Dict, Iterable, List, NamedTuple, Sequence, Tuple

import numpy as np
import torch

from model import RMSF_net_model


SOURCE_COMMIT = "eb353b0d9d909a6c46d2808467a9d0f16e2bbda7"
SPACING_ANGSTROM = 1.5
RESOLUTION_ANGSTROM = 4.0
SIGMA_FACTOR = 1.0 / (math.pi * math.sqrt(2.0))
CUTOFF_SIGMA = 5.0
EDGE_PADDING_ANGSTROM = 3.0 * RESOLUTION_ANGSTROM
BOX_SIZE = 40
CORE_SIZE = 10
BOX_HALO = (BOX_SIZE - CORE_SIZE) // 2

ATOMIC_NUMBERS = {
    "H": 1,
    "C": 6,
    "N": 7,
    "O": 8,
    "P": 15,
    "S": 16,
    "SE": 34,
}

OUTPUT_FIELDS = (
    "cohort",
    "sample_id",
    "residue_index",
    "chain",
    "author_residue_id",
    "residue_number",
    "insertion_code",
    "amino_acid_3",
    "rmsf_net_pdb_score",
    "preprocessing_seconds",
    "inference_seconds",
    "selected_box_count",
    "grid_voxel_count",
    "source_commit",
)


class Atom(NamedTuple):
    coord: np.ndarray
    atomic_number: int
    residue_key: Tuple[str, str, str]


class Residue(NamedTuple):
    chain: str
    number: str
    insertion: str
    name: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--output", required=True)
    parser.add_argument("--model", default="model_cpu/sim_model.pth")
    parser.add_argument("--limit", type=int)
    parser.add_argument("--batch-size", type=int, default=2)
    parser.add_argument("--threads", type=int, default=min(8, os.cpu_count() or 1))
    parser.add_argument("--shard-index", type=int, default=0)
    parser.add_argument("--shard-count", type=int, default=1)
    return parser.parse_args()


def read_manifest(
    path: Path, limit: int | None, shard_index: int, shard_count: int
) -> List[Dict[str, str]]:
    if shard_count < 1 or not 0 <= shard_index < shard_count:
        raise ValueError("shard-index must satisfy 0 <= index < shard-count")
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if limit is not None:
        rows = rows[:limit]
    return [row for index, row in enumerate(rows) if index % shard_count == shard_index]


def completed_samples(path: Path) -> set[str]:
    if not path.is_file() or path.stat().st_size == 0:
        return set()
    with path.open("r", encoding="utf-8", newline="") as handle:
        return {row["sample_id"] for row in csv.DictReader(handle)}


def infer_element(line: str) -> str:
    element = line[76:78].strip().upper() if len(line) >= 78 else ""
    if element:
        return element
    letters = "".join(character for character in line[12:16] if character.isalpha())
    if not letters:
        return "C"
    two = letters[:2].upper()
    return two if two in ATOMIC_NUMBERS else letters[0].upper()


def parse_pdb(path: str) -> Tuple[List[Atom], List[Residue]]:
    atoms: List[Atom] = []
    residue_names: Dict[Tuple[str, str, str], str] = {}
    ca_keys: List[Tuple[str, str, str]] = []
    seen_ca: set[Tuple[str, str, str]] = set()
    seen_atoms = False
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ENDMDL") and seen_atoms:
                break
            if not line.startswith("ATOM  ") or len(line) < 54:
                continue
            seen_atoms = True
            if line[16:17] not in (" ", "A"):
                continue
            atom_name = line[12:16].strip()
            chain = line[21:22].strip() or "_"
            number = line[22:26].strip()
            insertion = line[26:27].strip()
            residue_name = line[17:20].strip().upper()
            key = (chain, number, insertion)
            try:
                coord = np.asarray(
                    [float(line[30:38]), float(line[38:46]), float(line[46:54])],
                    dtype=np.float64,
                )
            except ValueError as error:
                raise ValueError(f"Malformed coordinates in {path}: {line.rstrip()}") from error
            element = infer_element(line)
            atoms.append(Atom(coord, ATOMIC_NUMBERS.get(element, 6), key))
            residue_names.setdefault(key, residue_name)
            if atom_name == "CA" and key not in seen_ca:
                seen_ca.add(key)
                ca_keys.append(key)
    if not atoms or len(ca_keys) < 2:
        raise ValueError(f"{path}: fewer than two protein CA residues")
    residues = [
        Residue(chain, number, insertion, residue_names[(chain, number, insertion)])
        for chain, number, insertion in ca_keys
    ]
    return atoms, residues


def grid_geometry(atoms: Sequence[Atom]) -> Tuple[np.ndarray, np.ndarray]:
    coordinates = np.stack([atom.coord for atom in atoms])
    lower = np.floor(
        (coordinates.min(axis=0) - EDGE_PADDING_ANGSTROM) / SPACING_ANGSTROM
    ) * SPACING_ANGSTROM
    upper = np.ceil(
        (coordinates.max(axis=0) + EDGE_PADDING_ANGSTROM) / SPACING_ANGSTROM
    ) * SPACING_ANGSTROM
    shape_xyz = np.rint((upper - lower) / SPACING_ANGSTROM).astype(int) + 1
    if np.any(shape_xyz < 2):
        raise ValueError(f"Invalid RMSF-net grid shape: {shape_xyz.tolist()}")
    return lower, shape_xyz


def simulated_density(
    atoms: Sequence[Atom], origin_xyz: np.ndarray, shape_xyz: np.ndarray
) -> np.ndarray:
    """Generate the documented Chimera molmap Gaussian density."""

    density = np.zeros(tuple(int(value) for value in shape_xyz[::-1]), dtype=np.float32)
    sigma = SIGMA_FACTOR * RESOLUTION_ANGSTROM
    cutoff = CUTOFF_SIGMA * sigma
    nx, ny, nz = (int(value) for value in shape_xyz)
    for atom in atoms:
        low = np.floor((atom.coord - cutoff - origin_xyz) / SPACING_ANGSTROM).astype(int)
        high = np.ceil((atom.coord + cutoff - origin_xyz) / SPACING_ANGSTROM).astype(int)
        low = np.maximum(low, 0)
        high = np.minimum(high, np.asarray([nx - 1, ny - 1, nz - 1]))
        xs = np.arange(low[0], high[0] + 1)
        ys = np.arange(low[1], high[1] + 1)
        zs = np.arange(low[2], high[2] + 1)
        dx2 = np.square(origin_xyz[0] + xs * SPACING_ANGSTROM - atom.coord[0])
        dy2 = np.square(origin_xyz[1] + ys * SPACING_ANGSTROM - atom.coord[1])
        dz2 = np.square(origin_xyz[2] + zs * SPACING_ANGSTROM - atom.coord[2])
        gaussian = np.exp(
            -(
                dz2[:, None, None]
                + dy2[None, :, None]
                + dx2[None, None, :]
            )
            / (2.0 * sigma * sigma)
        )
        density[np.ix_(zs, ys, xs)] += atom.atomic_number * gaussian.astype(np.float32)
    return density


def atom_indices_zyx(
    atoms: Sequence[Atom], origin_xyz: np.ndarray, shape_zyx: Sequence[int]
) -> np.ndarray:
    coordinates = np.stack([atom.coord for atom in atoms])
    indices_xyz = np.floor((coordinates - origin_xyz) / SPACING_ANGSTROM).astype(int)
    indices_zyx = indices_xyz[:, ::-1]
    maximum = np.asarray(shape_zyx, dtype=int) - 1
    return np.clip(indices_zyx, 0, maximum)


def selected_cells(indices_zyx: np.ndarray) -> List[Tuple[int, int, int]]:
    cells = {tuple(int(value) for value in row // CORE_SIZE) for row in indices_zyx}
    # Upstream increments axis 0 fastest, then axis 1, then axis 2.
    return sorted(cells, key=lambda cell: (cell[2], cell[1], cell[0]))


def extract_patch(padded: np.ndarray, cell: Sequence[int]) -> np.ndarray:
    starts = [int(value) * CORE_SIZE for value in cell]
    slices = tuple(slice(start, start + BOX_SIZE) for start in starts)
    patch = padded[slices]
    if patch.shape != (BOX_SIZE, BOX_SIZE, BOX_SIZE):
        raise ValueError(f"Unexpected RMSF-net patch shape {patch.shape}")
    return patch


def load_model(path: Path) -> RMSF_net_model:
    model = RMSF_net_model(in_channels=1)
    try:
        state = torch.load(path, map_location="cpu", weights_only=True)
    except TypeError:
        state = torch.load(path, map_location="cpu")
    model.load_state_dict(state)
    model.eval()
    return model


def predict_grid(
    model: RMSF_net_model,
    density: np.ndarray,
    cells: Sequence[Tuple[int, int, int]],
    batch_size: int,
) -> Tuple[np.ndarray, float]:
    cell_counts = np.ceil(np.asarray(density.shape) / CORE_SIZE).astype(int)
    covered_shape = cell_counts * CORE_SIZE
    back_padding = covered_shape - np.asarray(density.shape) + BOX_HALO
    padded = np.pad(
        density,
        tuple((BOX_HALO, int(back)) for back in back_padding),
        mode="constant",
    )
    prediction = np.zeros(density.shape, dtype=np.float32)
    inference_seconds = 0.0
    with torch.inference_mode():
        for offset in range(0, len(cells), batch_size):
            batch_cells = cells[offset : offset + batch_size]
            batch = np.stack([extract_patch(padded, cell) for cell in batch_cells])
            tensor = torch.from_numpy(batch[:, None]).float()
            started = time.perf_counter()
            output = model(tensor).squeeze(1).cpu().numpy()
            inference_seconds += time.perf_counter() - started
            if output.shape[1:] != (CORE_SIZE, CORE_SIZE, CORE_SIZE):
                raise ValueError(f"Unexpected RMSF-net output shape {output.shape}")
            for cell, core in zip(batch_cells, output):
                starts = np.asarray(cell) * CORE_SIZE
                stops = np.minimum(starts + CORE_SIZE, np.asarray(density.shape))
                target = tuple(slice(int(left), int(right)) for left, right in zip(starts, stops))
                source = tuple(slice(0, int(right - left)) for left, right in zip(starts, stops))
                prediction[target] = core[source]
    if not np.isfinite(prediction).all():
        raise ValueError("RMSF-net returned non-finite voxel predictions")
    # The released code applies 10**x - 1 before mapping values to atoms.
    return np.expm1(prediction * np.float32(math.log(10.0))), inference_seconds


def residue_predictions(
    atoms: Sequence[Atom],
    residues: Sequence[Residue],
    indices_zyx: np.ndarray,
    prediction: np.ndarray,
) -> List[float]:
    per_residue: Dict[Tuple[str, str, str], List[float]] = {}
    for atom, index in zip(atoms, indices_zyx):
        per_residue.setdefault(atom.residue_key, []).append(float(prediction[tuple(index)]))
    values = []
    for residue in residues:
        key = (residue.chain, residue.number, residue.insertion)
        atom_values = per_residue.get(key, [])
        if not atom_values:
            raise ValueError(f"No RMSF-net atoms mapped to residue {key}")
        values.append(float(np.mean(atom_values)))
    if not np.isfinite(values).all():
        raise ValueError("RMSF-net returned non-finite residue predictions")
    return values


def run_sample(
    model: RMSF_net_model, pdb_path: str, batch_size: int
) -> Tuple[List[Residue], List[float], Dict[str, float | int]]:
    started = time.perf_counter()
    atoms, residues = parse_pdb(pdb_path)
    origin, shape_xyz = grid_geometry(atoms)
    density = simulated_density(atoms, origin, shape_xyz)
    indices = atom_indices_zyx(atoms, origin, density.shape)
    cells = selected_cells(indices)
    preprocessing_seconds = time.perf_counter() - started
    prediction, inference_seconds = predict_grid(model, density, cells, batch_size)
    values = residue_predictions(atoms, residues, indices, prediction)
    metadata: Dict[str, float | int] = {
        "preprocessing_seconds": preprocessing_seconds,
        "inference_seconds": inference_seconds,
        "selected_box_count": len(cells),
        "grid_voxel_count": int(density.size),
    }
    return residues, values, metadata


def main() -> None:
    args = parse_args()
    if args.batch_size < 1 or args.threads < 1:
        raise ValueError("batch-size and threads must be positive")
    torch.set_num_threads(args.threads)
    manifest = read_manifest(
        Path(args.manifest), args.limit, args.shard_index, args.shard_count
    )
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    done = completed_samples(output)
    append = output.is_file() and output.stat().st_size > 0
    model = load_model(Path(args.model))
    with output.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_FIELDS)
        if not append:
            writer.writeheader()
        for position, record in enumerate(manifest, start=1):
            sample_id = record["sample_id"]
            if sample_id in done:
                continue
            pdb_path = record.get("container_pdb_path") or record["pdb_path"]
            residues, values, metadata = run_sample(model, pdb_path, args.batch_size)
            for residue_index, (residue, value) in enumerate(zip(residues, values)):
                writer.writerow(
                    {
                        "cohort": record["cohort"],
                        "sample_id": sample_id,
                        "residue_index": residue_index,
                        "chain": residue.chain,
                        "author_residue_id": residue.number + residue.insertion,
                        "residue_number": residue.number,
                        "insertion_code": residue.insertion,
                        "amino_acid_3": residue.name,
                        "rmsf_net_pdb_score": value,
                        **metadata,
                        "source_commit": SOURCE_COMMIT,
                    }
                )
            handle.flush()
            print(
                f"[{position}/{len(manifest)}] {sample_id}: {len(residues)} residues, "
                f"{metadata['selected_box_count']} boxes, "
                f"preprocess {metadata['preprocessing_seconds']:.3f} s, "
                f"inference {metadata['inference_seconds']:.3f} s",
                flush=True,
            )


if __name__ == "__main__":
    main()
