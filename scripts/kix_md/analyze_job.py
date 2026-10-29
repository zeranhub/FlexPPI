"""Stream one independent KIX MD production repeat into RMSF/H-chi metrics.

Example (run only after all production segments of this repeat are closed):
  python analyze_job.py --protein-pdb prepared/apo/protein.pdb \
      --job-dir runs/apo/rep01

Default inputs: <job-dir>/production/segment_*.dcd, 100000 frames, 1 ps/frame.
Default outputs: <job-dir>/analysis/residue_metrics.csv and summary.json.
Only one independent repeat is analyzed. No water, equilibrium frames, other
repeats, clustering references, inferred experimental values, or predictions
are included. Requires numpy and MDAnalysis. --self-test uses synthetic data.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import re
import sys
import tempfile

import numpy as np


PROTOCOL = "kix_one_repeat_chainB_CA_first_frame_fit_full_repeat_joint120deg_v1"
# Exact frozen CHI_ATOMS from FlexPPI src/flexppi/constants.py, validated by
# --self-test against the repository when available. No symmetry correction or
# alternative rotamer boundaries are introduced in this validation protocol.
CHI_ATOMS = {
    "ARG": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "NE"), ("CG", "CD", "NE", "CZ")],
    "ASN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "ASP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "OD1")],
    "CYS": [("N", "CA", "CB", "SG")],
    "GLN": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")],
    "GLU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "OE1")],
    "HIS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "ND1")],
    "ILE": [("N", "CA", "CB", "CG1"), ("CA", "CB", "CG1", "CD1")],
    "LEU": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "LYS": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD"), ("CB", "CG", "CD", "CE"), ("CG", "CD", "CE", "NZ")],
    "MET": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "SD"), ("CB", "CG", "SD", "CE")],
    "PHE": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "PRO": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD")],
    "SER": [("N", "CA", "CB", "OG")],
    "THR": [("N", "CA", "CB", "OG1")],
    "TRP": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "TYR": [("N", "CA", "CB", "CG"), ("CA", "CB", "CG", "CD1")],
    "VAL": [("N", "CA", "CB", "CG1")],
}
EXPECTED_CHI_HASH = "6b18e55670eb7bbb7179fc7ec701486d52c0a55234d475da4e28cfee400771c6"
ALIASES = {name: "HIS" for name in ("HID", "HIE", "HIP", "HSD", "HSE", "HSP")}
AA1 = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))
RANGES = {"A": (291, 315), "B": (586, 672), "C": (839, 857)}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def chi_hash() -> str:
    return hashlib.sha256(json.dumps(CHI_ATOMS, sort_keys=True, separators=(",", ":")).encode()).hexdigest()


def canonical_name(name: str) -> str:
    return ALIASES.get(name.upper(), name.upper())


def dihedral_degrees(points: np.ndarray) -> np.ndarray:
    """Frozen ITrace formula; leading dimensions may include a frame batch."""
    b0 = -(points[..., 1, :] - points[..., 0, :])
    b1 = points[..., 2, :] - points[..., 1, :]
    b2 = points[..., 3, :] - points[..., 2, :]
    b1 = b1 / np.maximum(np.linalg.norm(b1, axis=-1, keepdims=True), 1e-8)
    v = b0 - np.sum(b0 * b1, axis=-1, keepdims=True) * b1
    w = b2 - np.sum(b2 * b1, axis=-1, keepdims=True) * b1
    x = np.sum(v * w, axis=-1)
    y = np.sum(np.cross(b1, v) * w, axis=-1)
    return np.degrees(np.arctan2(y, x))


def angle_bins(angles: np.ndarray) -> np.ndarray:
    """Match ITrace exactly, including +180 -> bin 2 (no modulo remapping)."""
    return np.clip(np.floor((angles + 180.0) / 120.0).astype(np.int8), 0, 2)


def entropy_from_counts(counts: np.ndarray) -> float:
    counts = counts[counts > 0].astype(np.float64)
    if not len(counts):
        return float("nan")
    probabilities = counts / counts.sum()
    return float(-np.sum(probabilities * np.log(np.maximum(probabilities, 1e-12))))


def align_batch(mobile: np.ndarray, reference: np.ndarray) -> np.ndarray:
    """Proper Kabsch rotations, independently fitting every frame's KIX CA."""
    centered = mobile - mobile.mean(axis=1, keepdims=True)
    ref_mean = reference.mean(axis=0)
    ref = reference - ref_mean
    covariance = np.einsum("bni,nj->bij", centered, ref)
    u, _, vt = np.linalg.svd(covariance)
    rotation = u @ vt
    negative = np.linalg.det(rotation) < 0
    u[negative, :, -1] *= -1
    rotation = u @ vt
    return centered @ rotation + ref_mean


class Moments:
    """Merge block central moments without retaining trajectory coordinates."""

    def __init__(self):
        self.n = 0
        self.mean = None
        self.m2 = None

    def update(self, values: np.ndarray) -> None:
        values = np.asarray(values, dtype=np.float64)
        size = len(values)
        if size == 0:
            return
        mean = values.mean(axis=0)
        m2 = np.sum((values - mean) ** 2, axis=0)
        if not self.n:
            self.n, self.mean, self.m2 = size, mean, m2
        else:
            total = self.n + size
            delta = mean - self.mean
            self.m2 += m2 + delta * delta * self.n * size / total
            self.mean += delta * size / total
            self.n = total

    def rmsf(self) -> np.ndarray:
        if self.n < 2:
            raise ValueError("At least two production frames are required")
        return np.sqrt(np.maximum(np.sum(self.m2 / self.n, axis=-1), 0.0))


def parse_topology(path: Path) -> dict:
    residues, atoms = {}, []
    for raw in path.read_text(encoding="ascii").splitlines():
        line = raw.ljust(80)
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        chain, number, insertion = line[21], int(line[22:26]), line[26].strip()
        raw_name, atom_name = line[17:20], line[12:16].strip()
        name = canonical_name(raw_name)
        if name not in AA1 or chain not in RANGES:
            raise ValueError(f"Noncanonical/nonprotein atom in protein-only topology: {chain}:{number}:{raw_name}")
        if line[16] not in (" ", "A"):
            raise ValueError("Alternative locations are not allowed in prepared topology")
        key = (chain, number, insertion)
        if key not in residues:
            residues[key] = {"chain": chain, "author_residue_number": number, "insertion_code": insertion,
                             "amino_acid_3": name, "topology_residue_name": raw_name,
                             "amino_acid": AA1[name], "atom_lookup": {}}
        residue = residues[key]
        if residue["amino_acid_3"] != name or atom_name in residue["atom_lookup"]:
            raise ValueError(f"Duplicate/conflicting residue atom: {key} {atom_name}")
        residue["atom_lookup"][atom_name] = len(atoms)
        atoms.append([float(line[30:38]), float(line[38:46]), float(line[46:54])])
    if not atoms or not any(key[0] == "B" for key in residues):
        raise ValueError("No KIX chain B in protein topology")
    residue_list = list(residues.values())
    quads, bonds, alias_events, kix_ca, kix_rows = [], set(), [], [], []
    previous = {}
    for row_index, (key, residue) in enumerate(residues.items()):
        lookup, name = residue["atom_lookup"], residue["amino_acid_3"]
        for required in ("N", "CA", "C", "O"):
            if required not in lookup:
                raise ValueError(f"Missing required backbone atom {key} {required}")
        if name == "ILE" and "CD1" not in lookup and "CD" in lookup:
            lookup["CD1"] = lookup["CD"]
            alias_events.append({"chain": key[0], "residue": key[1], "alias": "ILE CD -> CD1"})
        if name != residue["topology_residue_name"]:
            alias_events.append({"chain": key[0], "residue": key[1], "alias": residue["topology_residue_name"] + " -> " + name})
        definitions = CHI_ATOMS.get(name, [])
        local = []
        for definition in definitions:
            missing = [atom for atom in definition if atom not in lookup]
            if missing:
                raise ValueError(f"Incomplete side chain {key}: {missing}; repair before analysis")
            local.append(len(quads))
            indices = [lookup[atom] for atom in definition]
            quads.append(indices)
            for first, second in zip(indices, indices[1:]):
                bonds.add(tuple(sorted((first, second))))
        residue["chi_columns"] = local
        residue["n_chi_angles"] = len(local)
        for first, second in (("N", "CA"), ("CA", "C"), ("C", "O")):
            bonds.add(tuple(sorted((lookup[first], lookup[second]))))
        if key[0] in previous:
            previous_key, previous_c = previous[key[0]]
            if key[1] != previous_key[1] + 1 or key[2] or previous_key[2]:
                raise ValueError(f"Unexpected missing/reordered author residue within chain {key}")
            bonds.add(tuple(sorted((previous_c, lookup["N"]))))
        previous[key[0]] = (key, lookup["C"])
        if key[0] == "B":
            kix_ca.append(lookup["CA"])
            kix_rows.append(row_index)
    for chain in {key[0] for key in residues}:
        low, high = RANGES[chain]
        observed = [(key[1], key[2]) for key in residues if key[0] == chain]
        if observed != [(number, "") for number in range(low, high + 1)]:
            raise ValueError(f"Prepared author-number mapping differs from KIX protocol for chain {chain}")
    return {"residues": residue_list, "atom_count": len(atoms), "reference_positions": np.asarray(atoms),
            "chi_quads": np.asarray(quads, dtype=np.int64), "bond_pairs": np.asarray(sorted(bonds), dtype=np.int64),
            "kix_ca": np.asarray(kix_ca), "kix_rows": np.asarray(kix_rows), "alias_events": alias_events}


def segment_paths(production: Path) -> list[Path]:
    numbered = []
    for path in production.glob("segment_*.dcd"):
        match = re.fullmatch(r"segment_(\d+)\.dcd", path.name)
        if not match:
            raise ValueError(f"Unexpected segment filename: {path.name}")
        numbered.append((int(match.group(1)), path))
    numbered.sort()
    if not numbered:
        raise FileNotFoundError(f"No segment_*.dcd files in {production}")
    numbers = [number for number, _ in numbered]
    if numbers[0] not in (0, 1) or numbers != list(range(numbers[0], numbers[0] + len(numbers))):
        raise ValueError(f"Segments must be uniquely numbered and contiguous from 0 or 1: {numbers}")
    return [path for _, path in numbered]


def geometry_extrema(frames: np.ndarray, pairs: np.ndarray) -> tuple[float, float]:
    lengths = np.linalg.norm(frames[:, pairs[:, 0]] - frames[:, pairs[:, 1]], axis=-1)
    low, high = float(lengths.min()), float(lengths.max())
    if low < 0.5 or high > 3.0:
        raise ValueError(f"Broken/periodic protein geometry: checked covalent distances {low:.4f}..{high:.4f} Angstrom. No wrapping repair is silently applied.")
    return low, high


def analyze(protein_pdb: Path, job_dir: Path, production_dir: Path | None = None,
            output_dir: Path | None = None, expected_frames: int = 100000,
            frame_interval_ps: float = 1.0, batch_frames: int = 256,
            overwrite: bool = False) -> dict:
    import MDAnalysis
    from MDAnalysis.coordinates.DCD import DCDReader

    if expected_frames < 2 or frame_interval_ps <= 0 or not 1 <= batch_frames <= 4096:
        raise ValueError("Require expected_frames>=2, frame_interval_ps>0 and 1<=batch_frames<=4096")
    if chi_hash() != EXPECTED_CHI_HASH:
        raise ValueError("Frozen CHI_ATOMS protocol hash mismatch")
    protein_pdb = protein_pdb.resolve(strict=True)
    job_dir = job_dir.resolve(strict=True)
    production = (production_dir or job_dir / "production").resolve(strict=True)
    output = (output_dir or job_dir / "analysis").resolve()
    if any((output / name).exists() for name in ("summary.json", "residue_metrics.csv")) and not overwrite:
        raise FileExistsError("Analysis outputs already exist; use --overwrite only for intentional re-analysis")
    topology = parse_topology(protein_pdb)
    segments = segment_paths(production)
    protein_hash = sha256_file(protein_pdb)
    provenance, total_header_frames = [], 0
    for path in segments:
        before = path.stat()
        with DCDReader(str(path)) as reader:
            if reader.n_atoms != topology["atom_count"]:
                raise ValueError(f"PDB/DCD atom count mismatch in {path}")
            if len(reader) == 0:
                raise ValueError(f"Empty DCD segment: {path}")
            interval = float(reader.dt)
            if not np.isclose(interval, frame_interval_ps, rtol=5e-5, atol=1e-6):
                raise ValueError(f"DCD frame interval {interval} ps != requested {frame_interval_ps} ps: {path}")
            first_time, last_time = float(reader[0].time), float(reader[-1].time)
            total_header_frames += len(reader)
            provenance.append({"filename": path.name, "num_frames": len(reader), "frame_interval_ps": interval,
                               "dcd_first_time_ps": first_time, "dcd_last_time_ps": last_time,
                               "bytes": before.st_size, "mtime_ns": before.st_mtime_ns,
                               "sha256": sha256_file(path)})
    if total_header_frames != expected_frames:
        raise ValueError(f"Expected exactly {expected_frames} frames in this repeat; found {total_header_frames}. Incomplete or duplicate production is not accepted.")
    residues, quads = topology["residues"], topology["chi_quads"]
    histograms = [np.zeros(3 ** r["n_chi_angles"], dtype=np.int64) if r["n_chi_angles"] else None for r in residues]
    ca_moments, reference = Moments(), None
    minimum_bond, maximum_bond = float("inf"), 0.0
    split_moments = [Moments(), Moments()]
    split_histograms = [[np.zeros_like(h) if h is not None else None for h in histograms] for _ in range(2)]
    split_counts = [0, 0]
    split_boundary = expected_frames // 2
    frame_count = 0
    for path, meta in zip(segments, provenance):
        with DCDReader(str(path)) as reader:
            for start in range(0, len(reader), batch_frames):
                frames = np.stack([ts.positions.copy() for ts in reader[start:start + batch_frames]], axis=0).astype(np.float64)
                if not np.isfinite(frames).all():
                    raise ValueError(f"Nonfinite coordinate in {path} near frame {start}")
                low, high = geometry_extrema(frames, topology["bond_pairs"])
                minimum_bond, maximum_bond = min(low, minimum_bond), max(high, maximum_bond)
                ca = frames[:, topology["kix_ca"]]
                if reference is None:
                    reference = ca[0].copy()
                aligned = align_batch(ca, reference)
                ca_moments.update(aligned)
                # Compute each quadruple with the identical arctan2 and boundary
                # rule used in the training labels; histograms retain exact counts.
                angles = dihedral_degrees(frames[:, quads])
                if not np.isfinite(angles).all():
                    raise ValueError(f"Nonfinite chi angle in {path} near frame {start}")
                bins = angle_bins(angles)
                block_codes = {}
                for residue_index, residue in enumerate(residues):
                    columns = residue["chi_columns"]
                    if columns:
                        codes = np.sum(bins[:, columns].astype(np.int32) * (3 ** np.arange(len(columns))), axis=1)
                        histograms[residue_index] += np.bincount(codes, minlength=len(histograms[residue_index]))
                        block_codes[residue_index] = codes
                global_indices = frame_count + np.arange(len(frames))
                for half in (0, 1):
                    mask = global_indices < split_boundary if half == 0 else global_indices >= split_boundary
                    if not np.any(mask):
                        continue
                    split_moments[half].update(aligned[mask])
                    split_counts[half] += int(mask.sum())
                    for residue_index, codes in block_codes.items():
                        split_histograms[half][residue_index] += np.bincount(codes[mask], minlength=len(histograms[residue_index]))
                frame_count += len(frames)
        after = path.stat()
        if after.st_size != meta["bytes"] or after.st_mtime_ns != meta["mtime_ns"]:
            raise RuntimeError(f"Trajectory changed during analysis: {path}; wait for production to finish")
    if frame_count != expected_frames or ca_moments.n != expected_frames:
        raise AssertionError("Read frame count differs from validated DCD headers")
    if sha256_file(protein_pdb) != protein_hash:
        raise RuntimeError("Protein topology changed during analysis")
    rmsf = ca_moments.rmsf()
    rmsf_by_row = dict(zip(topology["kix_rows"].tolist(), rmsf.tolist()))
    rows = []
    for index, residue in enumerate(residues):
        histogram = histograms[index]
        if histogram is not None and int(histogram.sum()) != expected_frames:
            raise AssertionError("Incomplete side-chain histogram")
        rows.append({
            "state": job_dir.parent.name, "replicate": job_dir.name,
            "chain": residue["chain"], "author_residue_number": residue["author_residue_number"],
            "insertion_code": residue["insertion_code"], "amino_acid_3": residue["amino_acid_3"],
            "amino_acid": residue["amino_acid"], "topology_residue_name": residue["topology_residue_name"],
            "n_frames": frame_count, "frame_interval_ps": frame_interval_ps,
            "rmsf_angstrom": rmsf_by_row.get(index, ""),
            "rmsf_definition": "KIX-chain-B-CA-first-frame-fit-mean-fluctuation" if index in rmsf_by_row else "not_computed_for_partner",
            "n_chi_angles": residue["n_chi_angles"],
            "chi_entropy_nats": entropy_from_counts(histogram) if histogram is not None else "",
            "n_observed_joint_states": int(np.count_nonzero(histogram)) if histogram is not None else "",
            "hchi_available": histogram is not None,
            "hchi_informative_evaluation": histogram is not None and residue["amino_acid_3"] not in ("ALA", "GLY", "PRO"),
        })
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0]))
    writer.writeheader()
    writer.writerows(rows)
    csv_bytes = buffer.getvalue().encode("utf-8")
    kix_entropy = [r["chi_entropy_nats"] for r in rows if r["chain"] == "B" and r["hchi_informative_evaluation"]]
    half_check = {}
    if min(split_counts) >= 2:
        half_rmsf = [moments.rmsf() for moments in split_moments]
        entropy_differences = [
            abs(entropy_from_counts(split_histograms[1][i]) - entropy_from_counts(split_histograms[0][i]))
            for i, residue in enumerate(residues)
            if residue["chain"] == "B" and residue["n_chi_angles"] and residue["amino_acid_3"] != "PRO"
        ]
        half_check = {"half_frame_counts": split_counts,
                      "KIX_CA_RMSF_half1_vs_half2_MAE_angstrom": float(np.mean(np.abs(half_rmsf[1] - half_rmsf[0]))),
                      "KIX_Hchi_excluding_AGP_half1_vs_half2_MAE_nats": float(np.mean(entropy_differences)),
                      "interpretation": "Descriptive split-half sensitivity only; not an independent-replicate uncertainty estimate or convergence certificate."}
    summary = {
        "protocol": PROTOCOL, "state": job_dir.parent.name, "replicate": job_dir.name,
        "protein_pdb": str(protein_pdb), "protein_pdb_sha256": protein_hash,
        "production_directory": str(production), "script_sha256": sha256_file(Path(__file__)),
        "chi_atoms_sha256": chi_hash(), "numpy_version": np.__version__, "MDAnalysis_version": MDAnalysis.__version__,
        "n_frames": frame_count, "frame_interval_ps": frame_interval_ps,
        "nominal_sampled_production_ns": frame_count * frame_interval_ps / 1000.0,
        "first_to_last_sample_span_ns": (frame_count - 1) * frame_interval_ps / 1000.0,
        "protein_atom_count": topology["atom_count"], "protein_residue_count": len(residues),
        "KIX_CA_count": len(rmsf), "KIX_mean_residue_RMSF_angstrom": float(rmsf.mean()),
        "KIX_sqrt_mean_squared_residue_RMSF_angstrom": float(np.sqrt(np.mean(rmsf ** 2))),
        "KIX_mean_Hchi_excluding_AGP_nats": float(np.mean(kix_entropy)),
        "RMSF_definition": "FlexPPI-compatible geometric definition: all KIX chain-B CA fitted to first production frame, then sqrt(mean squared displacement about each aligned CA trajectory mean). Full repeat, not training window averaging. NOT paper dominant-cluster-reference RMSF.",
        "Hchi_definition": "Per independent repeat; each available chi in [-180,-60),[-60,60),[60,180], floor((angle+180)/120) clipped to [0,2]; joint base-3 rotamer state; -sum(p*ln(p)); no pseudocounts, symmetry corrections, or pooling of independent repeats.",
        "PBC_policy": "No rewrapping, unwrapping, or silent repair; check all selected backbone/peptide/chi covalent distances in every frame before fitting KIX independently.",
        "segment_time_policy": "Numeric contiguous segment ordering; each file consumed exactly once. DCD header clocks may restart per segment; they are reported but are not treated as a globally continuous clock. Closed-file production bookkeeping remains the runner's responsibility.",
        "checked_covalent_distance_min_angstrom": minimum_bond,
        "checked_covalent_distance_max_angstrom": maximum_bond,
        "atom_name_aliases": topology["alias_events"], "segments": provenance,
        "split_half_sensitivity": half_check,
        "limitations": ["One MD repeat is not a five-repeat estimate.", "MD observables are simulation references, not experimental truth.", "Frame count does not establish convergence.", "Hchi is not thermodynamic entropy or methyl order parameter S2.", "RMSF is computed for chain B only; partners intentionally have blank RMSF fields."],
        "output_csv_sha256": hashlib.sha256(csv_bytes).hexdigest(),
    }
    output.mkdir(parents=True, exist_ok=True)
    for filename, data in (("residue_metrics.csv", csv_bytes), ("summary.json", (json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n").encode("utf-8"))):
        temporary = output / (filename + ".tmp")
        with temporary.open("wb") as handle:
            handle.write(data)
        os.replace(temporary, output / filename)
    return summary


def self_test() -> dict:
    """Numerical parity + synthetic DCD round trip, never analyze real production."""
    if chi_hash() != EXPECTED_CHI_HASH:
        raise AssertionError("CHI_ATOMS snapshot changed")
    angles = np.array([-180, -60.000001, -60, 59.999999, 60, 180], dtype=np.float64)
    np.testing.assert_array_equal(angle_bins(angles), [0, 0, 1, 1, 2, 2])
    for alias in ALIASES:
        assert canonical_name(alias) == "HIS"
    assert canonical_name("ILE") == "ILE"
    rng = np.random.default_rng(20260914)
    reference = rng.normal(size=(87, 3))
    mobile = reference + rng.normal(scale=.2, size=(65, 87, 3))
    mobile += rng.normal(scale=20, size=(65, 1, 3))
    fitted = align_batch(mobile, mobile[0])
    running = Moments()
    for start in range(0, len(fitted), 7):
        running.update(fitted[start:start + 7])
    expected = np.sqrt(np.mean(np.sum((fitted - fitted.mean(axis=0)) ** 2, axis=-1), axis=0))
    np.testing.assert_allclose(running.rmsf(), expected, rtol=1e-12, atol=1e-12)
    counts = np.zeros(81, dtype=np.int64)
    synthetic_codes = rng.integers(0, 81, size=1234)
    for start in range(0, len(synthetic_codes), 73):
        counts += np.bincount(synthetic_codes[start:start + 73], minlength=81)
    _, direct_counts = np.unique(synthetic_codes, return_counts=True)
    probabilities = direct_counts / direct_counts.sum()
    np.testing.assert_allclose(entropy_from_counts(counts), -np.sum(probabilities * np.log(probabilities)), rtol=1e-14)
    np.testing.assert_allclose(entropy_from_counts(np.ones(9, dtype=int)), np.log(9), rtol=1e-14)
    local_source = "not present in standalone deployment; frozen hash verified"
    root = Path(__file__).resolve().parents[2]
    constants_path = root / "src/flexppi/constants.py"
    if constants_path.is_file():
        spec = importlib.util.spec_from_file_location("kix_selftest_source_constants", constants_path)
        constants = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(constants)
        assert CHI_ATOMS == constants.CHI_ATOMS
        sys.path.insert(0, str(root / "src"))
        sys.path.insert(0, str(root / "scripts"))
        from extract_itrace_hchi import dihedral_degrees as original_dihedral, entropy as original_entropy
        from flexppi.data.labels import _kabsch_fit
        points = rng.normal(size=(33, 17, 4, 3))
        np.testing.assert_allclose(dihedral_degrees(points), np.stack([original_dihedral(frame) for frame in points]), rtol=1e-14, atol=1e-14)
        np.testing.assert_allclose(entropy_from_counts(counts), original_entropy(synthetic_codes), rtol=1e-14)
        np.testing.assert_allclose(fitted, np.stack([_kabsch_fit(frame, mobile[0]) for frame in mobile]), rtol=1e-12, atol=1e-12)
        local_source = "CHI_ATOMS, ITrace dihedral/entropy, and FlexPPI Kabsch numerical parity passed"
    import MDAnalysis as mda
    from MDAnalysis.coordinates.DCD import DCDReader, DCDWriter
    with tempfile.TemporaryDirectory(prefix="kix-analysis-selftest-") as temporary:
        path = Path(temporary) / "synthetic.dcd"
        universe = mda.Universe.empty(8, trajectory=True)
        original = []
        with DCDWriter(str(path), 8, dt=1.0) as writer:
            for index in range(6):
                universe.atoms.positions = rng.normal(size=(8, 3)) + index
                universe.dimensions = [100, 100, 100, 90, 90, 90]
                original.append(universe.atoms.positions.copy())
                writer.write(universe.atoms)
        with DCDReader(str(path)) as reader:
            assert len(reader) == 6 and reader.n_atoms == 8
            np.testing.assert_allclose(reader.dt, 1.0, rtol=1e-6)
            blocks = [np.stack([ts.positions.copy() for ts in reader[start:start + 2]]) for start in range(0, len(reader), 2)]
            np.testing.assert_allclose(np.concatenate(blocks), np.stack(original), rtol=1e-6, atol=1e-6)
    return {"self_test": "passed", "data_type": "synthetic numerical fixtures only; not MD results",
            "bin_boundaries": "passed", "online_RMSF": "passed", "joint_entropy": "passed",
            "histidine_aliases": "passed", "DCD_chunk_read_round_trip": "passed",
            "local_FlexPPI_protocol_parity": local_source, "chi_atoms_sha256": chi_hash()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--protein-pdb", type=Path, help="Prepared protein-only topology in exact DCD atom order")
    parser.add_argument("--job-dir", type=Path, help="One runs/<state>/repXX directory")
    parser.add_argument("--production-dir", type=Path, help="Override <job-dir>/production")
    parser.add_argument("--output-dir", type=Path, help="Override <job-dir>/analysis")
    parser.add_argument("--expected-frames", type=int, default=100000, help="Exact completed frame count, default 100000")
    parser.add_argument("--frame-interval-ps", type=float, default=1.0)
    parser.add_argument("--batch-frames", type=int, default=256, help="Maximum frames in memory; default 256")
    parser.add_argument("--overwrite", action="store_true", help="Intentionally replace this analysis's two output files")
    parser.add_argument("--self-test", action="store_true", help="Test numerical parity and synthetic DCD reading, not real production")
    args = parser.parse_args()
    if args.self_test:
        print(json.dumps(self_test(), indent=2))
        return
    if args.protein_pdb is None or args.job_dir is None:
        parser.error("--protein-pdb and --job-dir are required unless --self-test is used")
    result = analyze(args.protein_pdb, args.job_dir, args.production_dir, args.output_dir,
                     args.expected_frames, args.frame_interval_ps, args.batch_frames, args.overwrite)
    print(json.dumps({"protocol": result["protocol"], "state": result["state"], "replicate": result["replicate"],
                      "n_frames": result["n_frames"], "KIX_CA_count": result["KIX_CA_count"],
                      "analysis_output": str((args.output_dir or args.job_dir / "analysis").resolve())}, indent=2))


if __name__ == "__main__":
    main()
