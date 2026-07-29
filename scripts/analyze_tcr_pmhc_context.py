"""Run FlexPPI context comparisons for the 7L1D and 7RRG TCR-pMHC assemblies.

RCSB legacy biological-assembly PDB files can encode symmetry-generated
components as separate MODEL records.  FlexPPI intentionally reads only the
first MODEL of a PDB (the usual behaviour for NMR ensembles), so this script
first merges the non-overlapping A-E assembly chains into one standard PDB.

For each entry it then compares two fixed-coordinate inputs:

* whole: HLA-A3 heavy chain + beta-2-microglobulin + peptide + TCR alpha/beta
* tcr_peptide: peptide + TCR alpha/beta, with HLA chains removed

FlexPPI predicts residue-level flexibility properties; it does not generate
new atomic coordinates.  PyMOL outputs therefore colour the experimental
coordinates by predicted RMSF or by the whole-minus-subset context response.
"""

from __future__ import annotations

import argparse
import csv
import gc
import hashlib
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, MutableMapping, Sequence, Tuple

import numpy as np
import torch

from flexppi.constants import THREE_TO_ONE
from flexppi.data.dataset import collate_states
from flexppi.data.structures import ProteinStructure, parse_pdb
from flexppi.evaluation import _rankdata
from flexppi.inference import _structure_example, load_model
from flexppi.trainer import move_to_device


ResidueKey = Tuple[str, str]
METRICS = ("rmsf", "chi_entropy", "sasa")
WHOLE_CHAINS = ("A", "B", "C", "D", "E")
SUBSET_CHAINS = ("C", "D", "E")
CHAIN_ROLES = {
    "A": "HLA-A*03:01 heavy chain",
    "B": "beta-2-microglobulin",
    "C": "mutant PIK3CA peptide ALHGGWTTK",
    "D": "TCR alpha chain",
    "E": "TCR beta chain",
}


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _atom_key(line: str) -> Tuple[str, str, str, str, str, str]:
    return (
        line[21:22],
        line[22:26],
        line[26:27],
        line[17:20],
        line[12:16],
        line[16:17],
    )


def _renumber_atom(line: str, serial: int) -> str:
    value = line.rstrip("\r\n").ljust(80)
    return value[:6] + f"{serial:5d}" + value[11:] + "\n"


def _header_before_models(lines: Sequence[str]) -> List[str]:
    header: List[str] = []
    for line in lines:
        record = line[:6].strip()
        if record in {"MODEL", "ATOM", "HETATM"}:
            break
        if record not in {"END", "MASTER", "CONECT"}:
            header.append(line if line.endswith("\n") else line + "\n")
    return header


def merge_assembly_models(source: Path, destination: Path) -> Dict[str, Any]:
    """Merge non-overlapping assembly MODEL records into a single PDB model."""

    lines = source.read_text(encoding="utf-8", errors="replace").splitlines(True)
    model_count = sum(1 for line in lines if line.startswith("MODEL"))
    atoms = [line for line in lines if line.startswith("ATOM  ")]
    if not atoms:
        raise ValueError(f"No ATOM records found in {source}")
    seen: set = set()
    duplicates: List[Tuple[str, str, str, str, str, str]] = []
    for line in atoms:
        key = _atom_key(line)
        if key in seen:
            duplicates.append(key)
        seen.add(key)
    if duplicates:
        raise ValueError(
            "Assembly MODEL records contain overlapping atom identities; refusing to merge: "
            + repr(duplicates[:5])
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    header = _header_before_models(lines)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        handle.writelines(header)
        handle.write(
            "REMARK 900 RCSB BIOLOGICAL ASSEMBLY MODELS MERGED INTO ONE COORDINATE MODEL\n"
        )
        handle.write("REMARK 900 ONLY PROTEIN ATOM RECORDS RETAINED; COORDINATES UNCHANGED\n")
        serial = 1
        previous_chain = None
        for line in atoms:
            chain = line[21:22]
            if previous_chain is not None and chain != previous_chain:
                handle.write(f"TER   {serial:5d}      {line[17:20]} {previous_chain}\n")
                serial += 1
            handle.write(_renumber_atom(line, serial))
            serial += 1
            previous_chain = chain
        handle.write("END\n")

    structure = parse_pdb(str(destination))
    return {
        "source": str(source.resolve()),
        "source_sha256": _sha256(source),
        "destination": str(destination.resolve()),
        "source_model_count": model_count,
        "atom_count": len(atoms),
        "residue_count": structure.num_residues,
        "chains": _chain_summary(structure),
        "sha256": _sha256(destination),
    }


def extract_chains(source: Path, destination: Path, wanted: Sequence[str]) -> Dict[str, Any]:
    wanted_set = set(wanted)
    lines = source.read_text(encoding="utf-8", errors="replace").splitlines(True)
    atoms = [
        line for line in lines
        if line.startswith("ATOM  ") and line[21:22].strip() in wanted_set
    ]
    if not atoms:
        raise ValueError(f"No requested chains {wanted} found in {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        for line in _header_before_models(lines):
            handle.write(line)
        retained = ",".join(str(chain) for chain in wanted)
        handle.write(
            f"REMARK 900 FLEXPPI FIXED-COORDINATE SUBSET: RETAINED CHAINS {retained} ONLY\n"
        )
        serial = 1
        previous_chain = None
        for line in atoms:
            chain = line[21:22]
            if previous_chain is not None and chain != previous_chain:
                handle.write(f"TER   {serial:5d}      {line[17:20]} {previous_chain}\n")
                serial += 1
            handle.write(_renumber_atom(line, serial))
            serial += 1
            previous_chain = chain
        handle.write("END\n")
    structure = parse_pdb(str(destination))
    return {
        "source": str(source.resolve()),
        "destination": str(destination.resolve()),
        "residue_count": structure.num_residues,
        "chains": _chain_summary(structure),
        "sha256": _sha256(destination),
    }


def _chain_summary(structure: ProteinStructure) -> List[Dict[str, Any]]:
    return [
        {
            "chain": chain,
            "role": CHAIN_ROLES.get(chain, "unknown"),
            "residue_count": int(len(structure.indices_for_chain(chain))),
            "sequence": structure.chain_sequence(chain),
        }
        for chain in structure.unique_chains()
    ]


@torch.no_grad()
def predict_state(
    model: torch.nn.Module,
    model_config: Any,
    device: torch.device,
    pdb_path: Path,
    temperature: float,
    window_ns: float,
) -> Dict[str, Any]:
    structure = parse_pdb(str(pdb_path))
    batch = collate_states([
        _structure_example(
            str(pdb_path),
            structure.unique_chains(),
            temperature,
            window_ns,
            0.0,
        )
    ])
    batch = move_to_device(batch, device)
    output = model.forward_state(batch)
    valid = (batch["target_mask"][0] & output["mask"][0]).cpu().numpy()
    mean = output["mean"][0].cpu().numpy().astype(np.float64)
    std = torch.exp(0.5 * output["log_variance"][0]).cpu().numpy().astype(np.float64)
    rows: List[Dict[str, Any]] = []
    for index in np.flatnonzero(valid):
        row: Dict[str, Any] = {
            "residue_index": int(index),
            "chain": structure.chain_ids[index],
            "chain_role": CHAIN_ROLES.get(structure.chain_ids[index], "unknown"),
            "author_residue_id": structure.author_residue_ids[index],
            "amino_acid": structure.sequence[index],
        }
        for metric_index, metric in enumerate(model_config.metrics):
            row[metric] = float(mean[index, metric_index])
            row[metric + "_std"] = float(std[index, metric_index])
        rows.append(row)
    result = {
        "mode": "absolute_state",
        "input_pdb": str(pdb_path.resolve()),
        "residue_count": structure.num_residues,
        "valid_prediction_count": len(rows),
        "metrics": list(model_config.metrics),
        "temperature_kelvin": temperature,
        "window_ns_condition": window_ns,
        "dataset_code": 0.0,
        "chains": _chain_summary(structure),
        "rows": rows,
    }
    del batch, output
    gc.collect()
    return result


def _residue_key(row: Mapping[str, Any]) -> ResidueKey:
    return str(row["chain"]), str(row["author_residue_id"])


def _prediction_map(prediction: Mapping[str, Any]) -> Dict[ResidueKey, Mapping[str, Any]]:
    return {_residue_key(row): row for row in prediction["rows"]}


def _write_json(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)


def _write_csv(path: Path, rows: Sequence[Mapping[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames: List[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _chain_percentiles(rows: Sequence[Mapping[str, Any]], metric: str) -> Dict[ResidueKey, float]:
    result: Dict[ResidueKey, float] = {}
    chains = list(dict.fromkeys(str(row["chain"]) for row in rows))
    for chain in chains:
        selected = [row for row in rows if str(row["chain"]) == chain]
        values = np.asarray([float(row[metric]) for row in selected], dtype=np.float64)
        if len(values) == 1:
            percentiles = np.asarray([50.0])
        else:
            percentiles = 100.0 * (_rankdata(values) - 1.0) / (len(values) - 1.0)
        for row, percentile in zip(selected, percentiles):
            result[_residue_key(row)] = float(percentile)
    return result


def _write_bfactor_pdb(
    source: Path,
    destination: Path,
    values: Mapping[ResidueKey, float],
    remarks: Sequence[str],
    default_value: float | None = None,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True)
    with source.open("r", encoding="utf-8", errors="replace") as input_handle, destination.open(
        "w", encoding="utf-8", newline=""
    ) as output_handle:
        for remark in remarks:
            output_handle.write("REMARK 900 " + remark[:68] + "\n")
        for line in input_handle:
            if not line.startswith("ATOM  ") or len(line) < 66:
                output_handle.write(line)
                continue
            residue_name = line[17:20].strip().upper()
            if residue_name not in THREE_TO_ONE:
                output_handle.write(line)
                continue
            key = (line[21:22].strip() or "_", line[22:26].strip() + line[26:27].strip())
            value = values.get(key, default_value)
            if value is None or not np.isfinite(value):
                output_handle.write(line)
                continue
            padded = line.rstrip("\r\n").ljust(66)
            output_handle.write(padded[:60] + f"{float(value):6.2f}" + padded[66:] + "\n")


def _ca_coordinates(structure: ProteinStructure) -> Dict[ResidueKey, np.ndarray]:
    coordinates: Dict[ResidueKey, np.ndarray] = {}
    for index in range(structure.num_residues):
        if structure.backbone_mask[index, 1]:
            coordinates[(structure.chain_ids[index], structure.author_residue_ids[index])] = (
                structure.backbone[index, 1].astype(np.float64)
            )
    return coordinates


def _interface_flags(structure: ProteinStructure, cutoff: float = 8.0) -> Dict[ResidueKey, Dict[str, bool]]:
    coordinates = _ca_coordinates(structure)
    mhc = np.stack([value for (chain, _), value in coordinates.items() if chain in {"A", "B"}])
    peptide = np.stack([value for (chain, _), value in coordinates.items() if chain == "C"])
    flags: Dict[ResidueKey, Dict[str, bool]] = {}
    for key, coordinate in coordinates.items():
        chain, _ = key
        if chain not in SUBSET_CHAINS:
            continue
        mhc_contact = bool(np.any(np.linalg.norm(mhc - coordinate[None, :], axis=1) <= cutoff))
        peptide_contact = chain == "C" or bool(
            np.any(np.linalg.norm(peptide - coordinate[None, :], axis=1) <= cutoff)
        )
        flags[key] = {
            "ca_contact_with_mhc_at_8A": mhc_contact,
            "ca_contact_with_peptide_at_8A": peptide_contact,
        }
    return flags


def build_context_rows(
    pdb_id: str,
    whole: Mapping[str, Any],
    subset: Mapping[str, Any],
    structure: ProteinStructure,
) -> List[Dict[str, Any]]:
    whole_map = _prediction_map(whole)
    subset_map = _prediction_map(subset)
    flags = _interface_flags(structure)
    rows: List[Dict[str, Any]] = []
    for key in sorted(set(whole_map) & set(subset_map), key=lambda item: (item[0], item[1])):
        if key[0] not in SUBSET_CHAINS:
            continue
        whole_row = whole_map[key]
        subset_row = subset_map[key]
        row: Dict[str, Any] = {
            "pdb_id": pdb_id,
            "chain": key[0],
            "chain_role": CHAIN_ROLES[key[0]],
            "author_residue_id": key[1],
            "amino_acid": whole_row["amino_acid"],
            **flags.get(key, {}),
        }
        for metric in METRICS:
            whole_value = float(whole_row[metric])
            subset_value = float(subset_row[metric])
            combined_std = math.sqrt(
                float(whole_row[metric + "_std"]) ** 2
                + float(subset_row[metric + "_std"]) ** 2
            )
            delta = whole_value - subset_value
            row["whole_" + metric] = whole_value
            row["whole_" + metric + "_std"] = float(whole_row[metric + "_std"])
            row["tcr_peptide_" + metric] = subset_value
            row["tcr_peptide_" + metric + "_std"] = float(subset_row[metric + "_std"])
            row["delta_whole_minus_tcr_peptide_" + metric] = delta
            row["delta_over_combined_std_" + metric] = delta / combined_std if combined_std else math.nan
        rows.append(row)
    return rows


def _safe_correlation(first: np.ndarray, second: np.ndarray, rank: bool = False) -> float:
    valid = np.isfinite(first) & np.isfinite(second)
    if int(valid.sum()) < 3:
        return math.nan
    x = first[valid]
    y = second[valid]
    if rank:
        x = _rankdata(x)
        y = _rankdata(y)
    if float(np.std(x)) == 0.0 or float(np.std(y)) == 0.0:
        return math.nan
    return float(np.corrcoef(x, y)[0, 1])


def _numeric_summary(values: Iterable[float]) -> Dict[str, Any]:
    array = np.asarray(list(values), dtype=np.float64)
    array = array[np.isfinite(array)]
    if not len(array):
        return {"n": 0}
    return {
        "n": int(len(array)),
        "mean": float(array.mean()),
        "median": float(np.median(array)),
        "mean_absolute": float(np.abs(array).mean()),
        "minimum": float(array.min()),
        "maximum": float(array.max()),
    }


def summarize_context(rows: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    groups: Dict[str, List[Mapping[str, Any]]] = {
        "all_tcr_peptide": list(rows),
        "peptide": [row for row in rows if row["chain"] == "C"],
        "tcr": [row for row in rows if row["chain"] in {"D", "E"}],
        "mhc_contacting_tcr_peptide": [
            row for row in rows if bool(row.get("ca_contact_with_mhc_at_8A"))
        ],
        "peptide_contacting_tcr": [
            row for row in rows
            if row["chain"] in {"D", "E"} and bool(row.get("ca_contact_with_peptide_at_8A"))
        ],
    }
    result: Dict[str, Any] = {"groups": {}}
    for group_name, selected in groups.items():
        group_result: Dict[str, Any] = {"residue_count": len(selected)}
        for metric in METRICS:
            whole_values = np.asarray([float(row["whole_" + metric]) for row in selected])
            subset_values = np.asarray([float(row["tcr_peptide_" + metric]) for row in selected])
            delta_name = "delta_whole_minus_tcr_peptide_" + metric
            group_result[metric] = {
                "whole": _numeric_summary(whole_values),
                "tcr_peptide": _numeric_summary(subset_values),
                "delta_whole_minus_tcr_peptide": _numeric_summary(
                    float(row[delta_name]) for row in selected
                ),
                "pearson_whole_vs_subset": _safe_correlation(whole_values, subset_values),
                "spearman_whole_vs_subset": _safe_correlation(
                    whole_values, subset_values, rank=True
                ),
            }
        result["groups"][group_name] = group_result
    return result


def compare_complexes(
    context_by_id: Mapping[str, Sequence[Mapping[str, Any]]],
) -> Dict[str, Any]:
    ids = list(context_by_id)
    if len(ids) != 2:
        raise ValueError("Exactly two complexes are required")
    first_id, second_id = ids
    first = {(str(row["chain"]), str(row["author_residue_id"])): row for row in context_by_id[first_id]}
    second = {(str(row["chain"]), str(row["author_residue_id"])): row for row in context_by_id[second_id]}
    peptide_keys = sorted(
        [key for key in set(first) & set(second) if key[0] == "C"],
        key=lambda key: int("".join(character for character in key[1] if character.isdigit()) or 0),
    )
    peptide_rows: List[Dict[str, Any]] = []
    for key in peptide_keys:
        row: Dict[str, Any] = {
            "chain": key[0],
            "author_residue_id": key[1],
            "amino_acid": first[key]["amino_acid"],
        }
        for metric in METRICS:
            for field_prefix in ("whole_", "tcr_peptide_", "delta_whole_minus_tcr_peptide_"):
                field = field_prefix + metric
                row[first_id + "_" + field] = float(first[key][field])
                row[second_id + "_" + field] = float(second[key][field])
        peptide_rows.append(row)

    metric_summary: Dict[str, Any] = {}
    for metric in METRICS:
        whole_first = np.asarray([row[first_id + "_whole_" + metric] for row in peptide_rows])
        whole_second = np.asarray([row[second_id + "_whole_" + metric] for row in peptide_rows])
        delta_first = np.asarray([
            row[first_id + "_delta_whole_minus_tcr_peptide_" + metric] for row in peptide_rows
        ])
        delta_second = np.asarray([
            row[second_id + "_delta_whole_minus_tcr_peptide_" + metric] for row in peptide_rows
        ])
        metric_summary[metric] = {
            "whole_peptide_pearson": _safe_correlation(whole_first, whole_second),
            "whole_peptide_spearman": _safe_correlation(whole_first, whole_second, rank=True),
            "context_delta_peptide_pearson": _safe_correlation(delta_first, delta_second),
            "context_delta_peptide_spearman": _safe_correlation(
                delta_first, delta_second, rank=True
            ),
            "context_delta_same_sign_fraction": float(
                np.mean(np.sign(delta_first) == np.sign(delta_second))
            ),
        }
    p6 = next((row for row in peptide_rows if row["author_residue_id"] == "6"), None)
    central = [row for row in peptide_rows if row["author_residue_id"] in {"4", "5", "6", "7"}]
    return {
        "complexes": ids,
        "matched_peptide_residue_count": len(peptide_rows),
        "peptide_rows": peptide_rows,
        "metric_summary": metric_summary,
        "pTrp6": p6,
        "central_peptide_positions_4_to_7": central,
    }


def _all_atom_coordinates(path: Path) -> Dict[Tuple[str, str, str], np.ndarray]:
    atoms: Dict[Tuple[str, str, str], np.ndarray] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.startswith("ATOM  ") or line[16:17] not in {" ", "A", "1"}:
            continue
        key = (
            line[21:22].strip() or "_",
            line[22:26].strip() + line[26:27].strip(),
            line[12:16].strip(),
        )
        if key in atoms:
            continue
        atoms[key] = np.asarray(
            [float(line[30:38]), float(line[38:46]), float(line[46:54])],
            dtype=np.float64,
        )
    return atoms


def compare_structural_geometry(reference_path: Path, mobile_path: Path) -> Dict[str, Any]:
    """Compare the two ternary complexes after fitting common HLA chain-A CA atoms."""

    reference = _all_atom_coordinates(reference_path)
    mobile = _all_atom_coordinates(mobile_path)
    fit_keys = sorted(
        set(key for key in reference if key[0] == "A" and key[2] == "CA")
        & set(key for key in mobile if key[0] == "A" and key[2] == "CA")
    )
    target = np.stack([reference[key] for key in fit_keys])
    moving = np.stack([mobile[key] for key in fit_keys])
    target_centroid = target.mean(axis=0)
    moving_centroid = moving.mean(axis=0)
    u_matrix, _, vt_matrix = np.linalg.svd(
        (moving - moving_centroid).T @ (target - target_centroid)
    )
    rotation = u_matrix @ vt_matrix
    if np.linalg.det(rotation) < 0:
        u_matrix[:, -1] *= -1
        rotation = u_matrix @ vt_matrix

    def transform(coordinate: np.ndarray) -> np.ndarray:
        return (coordinate - moving_centroid) @ rotation + target_centroid

    def rmsd(keys: Sequence[Tuple[str, str, str]]) -> float:
        first = np.stack([reference[key] for key in keys])
        second = np.stack([transform(mobile[key]) for key in keys])
        return float(np.sqrt(np.sum((first - second) ** 2) / len(keys)))

    def common_keys(predicate: Any) -> List[Tuple[str, str, str]]:
        return sorted(
            set(key for key in reference if predicate(key))
            & set(key for key in mobile if predicate(key))
        )

    peptide_backbone = common_keys(
        lambda key: key[0] == "C" and key[2] in {"N", "CA", "C", "O"}
    )
    ptrp6_sidechain = common_keys(
        lambda key: key[0] == "C"
        and key[1] == "6"
        and key[2] not in {"N", "CA", "C", "O"}
    )
    pmhc_ca = common_keys(lambda key: key[0] in {"A", "B", "C"} and key[2] == "CA")

    def residue_number(value: str) -> int:
        digits = "".join(character for character in value if character.isdigit())
        return int(digits) if digits else 9999

    reference_peptide_centroid = np.mean(
        [coordinate for (chain, _, atom), coordinate in reference.items()
         if chain == "C" and atom == "CA"],
        axis=0,
    )
    mobile_peptide_centroid = np.mean(
        [transform(coordinate) for (chain, _, atom), coordinate in mobile.items()
         if chain == "C" and atom == "CA"],
        axis=0,
    )
    reference_tcr_centroid = np.mean(
        [coordinate for (chain, residue, atom), coordinate in reference.items()
         if chain in {"D", "E"} and atom == "CA" and residue_number(residue) <= 120],
        axis=0,
    )
    mobile_tcr_centroid = np.mean(
        [transform(coordinate) for (chain, residue, atom), coordinate in mobile.items()
         if chain in {"D", "E"} and atom == "CA" and residue_number(residue) <= 120],
        axis=0,
    )
    reference_vector = reference_tcr_centroid - reference_peptide_centroid
    mobile_vector = mobile_tcr_centroid - mobile_peptide_centroid
    cosine = float(
        np.dot(reference_vector, mobile_vector)
        / (np.linalg.norm(reference_vector) * np.linalg.norm(mobile_vector))
    )
    docking_angle = math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))
    return {
        "reference": str(reference_path.resolve()),
        "mobile": str(mobile_path.resolve()),
        "fit_selection": "common HLA chain A CA atoms",
        "fit_atom_count": len(fit_keys),
        "hla_chain_A_ca_rmsd_angstrom": rmsd(fit_keys),
        "pmhc_ABC_ca_rmsd_angstrom": rmsd(pmhc_ca),
        "peptide_backbone_rmsd_angstrom": rmsd(peptide_backbone),
        "pTrp6_sidechain_rmsd_angstrom": rmsd(ptrp6_sidechain),
        "tcr_variable_domain_centroid_distance_angstrom": float(
            np.linalg.norm(reference_tcr_centroid - mobile_tcr_centroid)
        ),
        "tcr_approach_vector_angle_degrees": docking_angle,
    }


def _metric_maps(prediction: Mapping[str, Any], metric: str) -> Dict[ResidueKey, float]:
    return {_residue_key(row): float(row[metric]) for row in prediction["rows"]}


def _range(values: Iterable[float]) -> Tuple[float, float]:
    array = np.asarray(list(values), dtype=np.float64)
    array = array[np.isfinite(array)]
    return float(array.min()), float(array.max())


def write_context_pymol_script(
    pdb_id: str,
    whole_pdb: Path,
    subset_pdb: Path,
    delta_pdb: Path,
    output_path: Path,
    session_path: Path,
    whole_range: Tuple[float, float],
    subset_range: Tuple[float, float],
    delta_limit: float,
) -> None:
    shared_min = min(whole_range[0], subset_range[0])
    shared_max = max(whole_range[1], subset_range[1])
    prefix = pdb_id.replace("7", "PDB7")
    script = f'''# {pdb_id}: full TCR-pMHC context versus TCR+peptide context
# B-factor is FlexPPI predicted RMSF in Angstrom. Coordinates are experimental.
reinitialize
bg_color white
set ray_opaque_background, off
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
load "{whole_pdb.resolve().as_posix()}", {prefix}_whole
load "{subset_pdb.resolve().as_posix()}", {prefix}_tcr_peptide
load "{delta_pdb.resolve().as_posix()}", {prefix}_context_delta
hide everything, all
show cartoon, {prefix}_whole or {prefix}_tcr_peptide or {prefix}_context_delta
show sticks, ({prefix}_whole or {prefix}_tcr_peptide or {prefix}_context_delta) and chain C
show sticks, byres (({prefix}_whole or {prefix}_tcr_peptide or {prefix}_context_delta) and chain D+E within 5 of chain C)
hide cartoon, {prefix}_context_delta and chain A+B
spectrum b, blue_white_red, {prefix}_whole and chain C+D+E, minimum={shared_min:.4f}, maximum={shared_max:.4f}
spectrum b, blue_white_red, {prefix}_tcr_peptide, minimum={shared_min:.4f}, maximum={shared_max:.4f}
spectrum b, blue_white_red, {prefix}_context_delta and chain C+D+E, minimum={-delta_limit:.4f}, maximum={delta_limit:.4f}
select {prefix}_peptide, chain C
select {prefix}_pTrp6, chain C and resi 6
select {prefix}_TCR, chain D+E
select {prefix}_TCR_peptide_interface, byres ((chain D+E) within 5 of (chain C))
set grid_mode, 1
set grid_slot, 1, {prefix}_whole
set grid_slot, 2, {prefix}_tcr_peptide
set grid_slot, 3, {prefix}_context_delta
orient {prefix}_whole
zoom all, 4
set seq_view, 1
scene {prefix}_context_comparison, store
save "{session_path.resolve().as_posix()}"
'''
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(script, encoding="utf-8")


def write_between_complex_pymol_scripts(
    whole_pdbs: Mapping[str, Path],
    output_dir: Path,
    shared_range: Tuple[float, float],
) -> Dict[str, str]:
    first_id, second_id = list(whole_pdbs)
    first_name = "PDB7L1D_whole"
    second_name = "PDB7RRG_whole"
    flexibility_path = output_dir / "compare_7L1D_vs_7RRG_whole_flexibility.pml"
    flexibility_session = output_dir / "compare_7L1D_vs_7RRG_whole_flexibility.pse"
    flexibility = f'''# Compare FlexPPI whole-complex RMSF for 7L1D (TCR4) and 7RRG (TCR3)
reinitialize
bg_color white
set ray_opaque_background, off
set antialias, 2
load "{whole_pdbs[first_id].resolve().as_posix()}", {first_name}
load "{whole_pdbs[second_id].resolve().as_posix()}", {second_name}
align ({second_name} and chain A and name CA), ({first_name} and chain A and name CA), cycles=5
hide everything, all
show cartoon, {first_name} or {second_name}
show sticks, ({first_name} or {second_name}) and chain C
show sticks, byres (({first_name} or {second_name}) and chain D+E within 5 of chain C)
spectrum b, blue_white_red, {first_name}, minimum={shared_range[0]:.4f}, maximum={shared_range[1]:.4f}
spectrum b, blue_white_red, {second_name}, minimum={shared_range[0]:.4f}, maximum={shared_range[1]:.4f}
set grid_mode, 1
set grid_slot, 1, {first_name}
set grid_slot, 2, {second_name}
orient {first_name} and chain A+C+D+E
zoom all, 4
scene compare_whole_flexibility, store
save "{flexibility_session.resolve().as_posix()}"
'''
    flexibility_path.write_text(flexibility, encoding="utf-8")

    overlay_path = output_dir / "compare_7L1D_vs_7RRG_structural_overlay.pml"
    overlay_session = output_dir / "compare_7L1D_vs_7RRG_structural_overlay.pse"
    overlay = f'''# Structural overlay after HLA-A3 chain-A alignment
reinitialize
bg_color white
set ray_opaque_background, off
set antialias, 2
load "{whole_pdbs[first_id].resolve().as_posix()}", {first_name}
load "{whole_pdbs[second_id].resolve().as_posix()}", {second_name}
align ({second_name} and chain A and name CA), ({first_name} and chain A and name CA), cycles=5
hide everything, all
show cartoon, {first_name} or {second_name}
show sticks, ({first_name} or {second_name}) and chain C
color gray70, ({first_name} or {second_name}) and chain A+B
color cyan, {first_name} and chain D+E
color magenta, {second_name} and chain D+E
color yellow, {first_name} and chain C
color orange, {second_name} and chain C
show spheres, ({first_name} or {second_name}) and chain C and resi 6 and name CA
set sphere_scale, 0.35
orient ({first_name} or {second_name}) and chain A+C+D+E
zoom all, 4
scene structural_overlay, store
save "{overlay_session.resolve().as_posix()}"
'''
    overlay_path.write_text(overlay, encoding="utf-8")
    return {
        "whole_flexibility_pml": str(flexibility_path.resolve()),
        "whole_flexibility_session_after_running": str(flexibility_session.resolve()),
        "structural_overlay_pml": str(overlay_path.resolve()),
        "structural_overlay_session_after_running": str(overlay_session.resolve()),
    }


def _format_number(value: Any, digits: int = 3) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    return "NA" if not math.isfinite(number) else f"{number:.{digits}f}"


def write_report(
    output_path: Path,
    checkpoint: Path,
    preparation: Mapping[str, Any],
    context_summaries: Mapping[str, Any],
    cross_comparison: Mapping[str, Any],
    artifacts: Mapping[str, Any],
) -> None:
    p6 = cross_comparison.get("pTrp6") or {}
    rmsf_cross = cross_comparison["metric_summary"]["rmsf"]
    geometry = cross_comparison["structural_comparison_after_HLA_alignment"]
    lines = [
        "# 7L1D/7RRG FlexPPI TCR-pMHC context analysis",
        "",
        "## What was predicted",
        "",
        "FlexPPI was run with the static-pretraining checkpoint on four fixed-coordinate inputs: "
        "the complete A-E biological assembly and the C+D+E TCR+peptide subset for each PDB. "
        "The model predicts residue-level RMSF, side-chain conformational entropy, and SASA; "
        "it does not generate alternative coordinates or a molecular-dynamics trajectory.",
        "",
        f"Checkpoint: `{checkpoint.resolve()}`  ",
        f"Checkpoint SHA256: `{_sha256(checkpoint)}`",
        "",
        "## Input audit",
        "",
    ]
    for pdb_id, item in preparation.items():
        chain_text = ", ".join(
            f"{chain['chain']}={chain['residue_count']} aa" for chain in item["whole"]["chains"]
        )
        lines.append(
            f"- {pdb_id}: merged {item['whole']['source_model_count']} assembly MODEL records; "
            f"{item['whole']['residue_count']} residues ({chain_text})."
        )
    lines.extend([
        "",
        "## Main RMSF context result",
        "",
        "`delta = whole complex prediction - TCR+peptide-only prediction`. Negative values mean "
        "that adding HLA-A3/beta-2-microglobulin makes the model predict lower flexibility.",
        "",
        "| PDB | peptide mean delta (A) | TCR mean delta (A) | MHC-contacting C+D+E mean delta (A) |",
        "|---|---:|---:|---:|",
    ])
    for pdb_id, summary in context_summaries.items():
        groups = summary["groups"]
        lines.append(
            "| {} | {} | {} | {} |".format(
                pdb_id,
                _format_number(groups["peptide"]["rmsf"]["delta_whole_minus_tcr_peptide"].get("mean")),
                _format_number(groups["tcr"]["rmsf"]["delta_whole_minus_tcr_peptide"].get("mean")),
                _format_number(
                    groups["mhc_contacting_tcr_peptide"]["rmsf"]
                    ["delta_whole_minus_tcr_peptide"].get("mean")
                ),
            )
        )
    lines.extend([
        "",
        "## Cross-TCR peptide consistency",
        "",
        f"Across the nine peptide positions, the whole-complex RMSF profiles have Pearson "
        f"`{_format_number(rmsf_cross['whole_peptide_pearson'])}` and Spearman "
        f"`{_format_number(rmsf_cross['whole_peptide_spearman'])}`. The two context-delta "
        f"profiles have the same sign at `{_format_number(rmsf_cross['context_delta_same_sign_fraction'])}` "
        "of peptide positions.",
        "",
        "For pTrp6:",
        "",
        "| quantity | 7L1D (TCR4) | 7RRG (TCR3) |",
        "|---|---:|---:|",
        "| whole RMSF (A) | {} | {} |".format(
            _format_number(p6.get("7L1D_whole_rmsf")),
            _format_number(p6.get("7RRG_whole_rmsf")),
        ),
        "| TCR+peptide RMSF (A) | {} | {} |".format(
            _format_number(p6.get("7L1D_tcr_peptide_rmsf")),
            _format_number(p6.get("7RRG_tcr_peptide_rmsf")),
        ),
        "| whole-minus-subset delta (A) | {} | {} |".format(
            _format_number(p6.get("7L1D_delta_whole_minus_tcr_peptide_rmsf")),
            _format_number(p6.get("7RRG_delta_whole_minus_tcr_peptide_rmsf")),
        ),
        "",
        "## Structural comparison after HLA alignment",
        "",
        "7RRG was rigidly fitted to 7L1D using common HLA chain-A CA atoms. The fitted "
        "HLA CA RMSD is `{}` A, peptide-backbone RMSD is `{}` A, and pTrp6-side-chain "
        "RMSD is `{}` A. The TCR variable-domain approach vectors differ by `{}` degrees."
        .format(
            _format_number(geometry["hla_chain_A_ca_rmsd_angstrom"]),
            _format_number(geometry["peptide_backbone_rmsd_angstrom"]),
            _format_number(geometry["pTrp6_sidechain_rmsd_angstrom"]),
            _format_number(geometry["tcr_approach_vector_angle_degrees"]),
        ),
        "",
        "## Interpretation boundary",
        "",
        "This experiment can test whether FlexPPI detects an MHC-conditioned flexibility signal "
        "in the bound mutant peptide and whether that signal is consistent across TCR3 and TCR4. "
        "It cannot directly validate mutant-versus-wild-type selectivity because both 7L1D and "
        "7RRG contain the same mutant peptide ALHGGWTTK and there is no experimental WT-TCR "
        "ternary structure. Agreement is therefore supportive only of the narrower MHC-context "
        "part of the mechanism, not proof of the paper's complete dynamic-allostery theory.",
        "",
        "## PyMOL",
        "",
        "Run any `.pml` file in the `pymol` directory. The scripts load absolute RMSF in the "
        "PDB B-factor field, use the same colour range within a comparison, highlight peptide "
        "chain C and pTrp6, and save a `.pse` session after execution.",
        "",
        f"Artifact index: `{artifacts.get('metadata', '')}`",
        "",
    ])
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input-7l1d",
        default=(
            "artifacts/evaluation/7L1D_7RRG_tcr_pmhc_flexppi/"
            "inputs/biological_assembly/7L1D.pdb"
        ),
    )
    parser.add_argument(
        "--input-7rrg",
        default=(
            "artifacts/evaluation/7L1D_7RRG_tcr_pmhc_flexppi/"
            "inputs/biological_assembly/7RRG.pdb"
        ),
    )
    parser.add_argument("--checkpoint", default="artifacts/pretrain/best.pt")
    parser.add_argument(
        "--output-dir",
        default="artifacts/evaluation/7L1D_7RRG_tcr_pmhc_flexppi",
    )
    parser.add_argument("--device", default="auto")
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--window-ns", type=float, default=100.0)
    arguments = parser.parse_args()

    output_dir = Path(arguments.output_dir)
    prepared_dir = output_dir / "inputs" / "prepared"
    prediction_dir = output_dir / "predictions"
    analysis_dir = output_dir / "analysis"
    pymol_dir = output_dir / "pymol"
    for directory in (prepared_dir, prediction_dir, analysis_dir, pymol_dir):
        directory.mkdir(parents=True, exist_ok=True)

    sources = {"7L1D": Path(arguments.input_7l1d), "7RRG": Path(arguments.input_7rrg)}
    checkpoint = Path(arguments.checkpoint)
    preparation: Dict[str, Any] = {}
    prepared_paths: Dict[str, Dict[str, Path]] = {}
    for pdb_id, source in sources.items():
        whole_path = prepared_dir / f"{pdb_id}_biological_assembly_merged.pdb"
        subset_path = prepared_dir / f"{pdb_id}_TCR_peptide_CDE.pdb"
        whole_info = merge_assembly_models(source, whole_path)
        found_chains = {item["chain"] for item in whole_info["chains"]}
        if found_chains != set(WHOLE_CHAINS):
            raise ValueError(f"{pdb_id} expected chains A-E, found {sorted(found_chains)}")
        subset_info = extract_chains(whole_path, subset_path, SUBSET_CHAINS)
        preparation[pdb_id] = {"whole": whole_info, "tcr_peptide": subset_info}
        prepared_paths[pdb_id] = {"whole": whole_path, "tcr_peptide": subset_path}

    model, model_config, device = load_model(str(checkpoint), arguments.device)
    predictions: Dict[str, Dict[str, Dict[str, Any]]] = {}
    for pdb_id in sources:
        predictions[pdb_id] = {}
        for context in ("whole", "tcr_peptide"):
            result = predict_state(
                model,
                model_config,
                device,
                prepared_paths[pdb_id][context],
                arguments.temperature,
                arguments.window_ns,
            )
            predictions[pdb_id][context] = result
            context_dir = prediction_dir / pdb_id / context
            _write_json(context_dir / "prediction.json", result)
            _write_csv(context_dir / "residue_predictions.csv", result["rows"])

    context_rows_by_id: Dict[str, List[Dict[str, Any]]] = {}
    context_summaries: Dict[str, Any] = {}
    pymol_artifacts: Dict[str, Any] = {}
    whole_rmsf_pdbs: Dict[str, Path] = {}
    all_whole_rmsf: List[float] = []
    for pdb_id in sources:
        whole = predictions[pdb_id]["whole"]
        subset = predictions[pdb_id]["tcr_peptide"]
        whole_structure = parse_pdb(str(prepared_paths[pdb_id]["whole"]))
        context_rows = build_context_rows(pdb_id, whole, subset, whole_structure)
        context_rows_by_id[pdb_id] = context_rows
        context_summary = summarize_context(context_rows)
        context_summaries[pdb_id] = context_summary
        _write_csv(analysis_dir / f"{pdb_id}_whole_vs_tcr_peptide_residues.csv", context_rows)
        _write_json(analysis_dir / f"{pdb_id}_whole_vs_tcr_peptide_summary.json", context_summary)

        whole_rmsf = _metric_maps(whole, "rmsf")
        subset_rmsf = _metric_maps(subset, "rmsf")
        delta_rmsf = {
            _residue_key(row): float(row["delta_whole_minus_tcr_peptide_rmsf"])
            for row in context_rows
        }
        whole_percentile = _chain_percentiles(whole["rows"], "rmsf")
        subset_percentile = _chain_percentiles(subset["rows"], "rmsf")
        visual_dir = prediction_dir / pdb_id / "pymol_values"
        whole_rmsf_path = visual_dir / f"{pdb_id}_whole_FlexPPI_RMSF_A.pdb"
        subset_rmsf_path = visual_dir / f"{pdb_id}_TCR_peptide_FlexPPI_RMSF_A.pdb"
        delta_rmsf_path = visual_dir / f"{pdb_id}_whole_minus_TCR_peptide_RMSF_A.pdb"
        whole_percentile_path = visual_dir / f"{pdb_id}_whole_FlexPPI_RMSF_percentile.pdb"
        subset_percentile_path = visual_dir / f"{pdb_id}_TCR_peptide_FlexPPI_RMSF_percentile.pdb"
        _write_bfactor_pdb(
            prepared_paths[pdb_id]["whole"],
            whole_rmsf_path,
            whole_rmsf,
            ("B-factor = FlexPPI predicted RMSF in Angstrom; whole A-E context",),
        )
        _write_bfactor_pdb(
            prepared_paths[pdb_id]["tcr_peptide"],
            subset_rmsf_path,
            subset_rmsf,
            ("B-factor = FlexPPI predicted RMSF in Angstrom; C+D+E context",),
        )
        _write_bfactor_pdb(
            prepared_paths[pdb_id]["whole"],
            delta_rmsf_path,
            delta_rmsf,
            ("B-factor = whole minus C+D+E FlexPPI RMSF in Angstrom",),
            default_value=0.0,
        )
        _write_bfactor_pdb(
            prepared_paths[pdb_id]["whole"],
            whole_percentile_path,
            whole_percentile,
            ("B-factor = within-chain FlexPPI RMSF percentile; whole context",),
        )
        _write_bfactor_pdb(
            prepared_paths[pdb_id]["tcr_peptide"],
            subset_percentile_path,
            subset_percentile,
            ("B-factor = within-chain FlexPPI RMSF percentile; C+D+E context",),
        )
        whole_range = _range(whole_rmsf.values())
        subset_range = _range(subset_rmsf.values())
        delta_limit = max(abs(value) for value in delta_rmsf.values()) or 1.0
        context_pml = pymol_dir / f"compare_{pdb_id}_whole_vs_TCR_peptide.pml"
        context_pse = pymol_dir / f"compare_{pdb_id}_whole_vs_TCR_peptide.pse"
        write_context_pymol_script(
            pdb_id,
            whole_rmsf_path,
            subset_rmsf_path,
            delta_rmsf_path,
            context_pml,
            context_pse,
            whole_range,
            subset_range,
            delta_limit,
        )
        pymol_artifacts[pdb_id] = {
            "context_comparison_pml": str(context_pml.resolve()),
            "session_after_running": str(context_pse.resolve()),
            "whole_rmsf_pdb": str(whole_rmsf_path.resolve()),
            "tcr_peptide_rmsf_pdb": str(subset_rmsf_path.resolve()),
            "context_delta_rmsf_pdb": str(delta_rmsf_path.resolve()),
            "whole_percentile_pdb": str(whole_percentile_path.resolve()),
            "tcr_peptide_percentile_pdb": str(subset_percentile_path.resolve()),
        }
        whole_rmsf_pdbs[pdb_id] = whole_rmsf_path
        all_whole_rmsf.extend(whole_rmsf.values())

    cross = compare_complexes(context_rows_by_id)
    cross["structural_comparison_after_HLA_alignment"] = compare_structural_geometry(
        prepared_paths["7L1D"]["whole"],
        prepared_paths["7RRG"]["whole"],
    )
    _write_json(analysis_dir / "7L1D_vs_7RRG_cross_complex_summary.json", cross)
    _write_csv(analysis_dir / "7L1D_vs_7RRG_peptide_predictions.csv", cross["peptide_rows"])
    between_artifacts = write_between_complex_pymol_scripts(
        whole_rmsf_pdbs,
        pymol_dir,
        _range(all_whole_rmsf),
    )
    pymol_artifacts["between_complexes"] = between_artifacts

    metadata_path = output_dir / "run_metadata.json"
    metadata: Dict[str, Any] = {
        "analysis": "7L1D/7RRG FlexPPI full assembly versus TCR+peptide context",
        "checkpoint": str(checkpoint.resolve()),
        "checkpoint_sha256": _sha256(checkpoint),
        "device": str(device),
        "torch_version": torch.__version__,
        "python_version": sys.version,
        "temperature_kelvin": arguments.temperature,
        "window_ns_condition": arguments.window_ns,
        "preparation": preparation,
        "context_summaries": context_summaries,
        "cross_complex_summary": cross,
        "pymol": pymol_artifacts,
        "warnings": [
            "FlexPPI predicts residue properties on fixed coordinates; it does not generate structures.",
            "7L1D has more residues than the pretraining configuration max_residues=768.",
            "Both PDB entries contain the same mutant peptide; this is not a WT-versus-mutant test.",
            "TCR+peptide without MHC is an in-silico context ablation, not a stable experimental complex.",
        ],
        "metadata": str(metadata_path.resolve()),
    }
    _write_json(metadata_path, metadata)
    write_report(
        output_dir / "REPORT.md",
        checkpoint,
        preparation,
        context_summaries,
        cross,
        metadata,
    )
    print(json.dumps({
        "output_dir": str(output_dir.resolve()),
        "report": str((output_dir / "REPORT.md").resolve()),
        "metadata": str(metadata_path.resolve()),
        "pymol_scripts": pymol_artifacts,
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
