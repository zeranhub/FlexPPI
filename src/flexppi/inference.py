"""Checkpoint loading and state/pair prediction from static PDB files."""

from __future__ import annotations

import math
from importlib.resources import files
from pathlib import Path
from typing import Any, Dict, List, Optional, Sequence

import torch

from flexppi.config import ModelConfig
from flexppi.constants import METRIC_NAMES
from flexppi.data.dataset import collate_pairs, collate_states
from flexppi.data.pairing import match_apo_to_holo
from flexppi.data.structures import ProteinStructure, parse_pdb, target_mask
from flexppi.model.information_ablation import build_information_ablation
from flexppi.trainer import move_to_device, resolve_device


def _structure_example(
    path: str,
    target_chains: Sequence[str],
    temperature: float,
    window_ns: float,
    dataset_code: float,
    structure: Optional[ProteinStructure] = None,
) -> Dict[str, Any]:
    if structure is None:
        structure = parse_pdb(path)
    n_residues = structure.num_residues
    metric_count = len(METRIC_NAMES)
    return {
        "sample_id": path,
        "dataset": "inference",
        "aatype": torch.from_numpy(structure.aatype),
        "residue_index": torch.from_numpy(structure.residue_index),
        "chain_index": torch.from_numpy(structure.chain_index),
        "backbone": torch.from_numpy(structure.backbone),
        "backbone_mask": torch.from_numpy(structure.backbone_mask),
        "geometry_mask": torch.from_numpy(structure.frame_mask),
        "target_mask": torch.from_numpy(target_mask(structure, target_chains)),
        "condition": torch.tensor(
            [
                (temperature - 300.0) / 20.0,
                math.log1p(max(window_ns, 0.0)) / math.log1p(500.0),
                dataset_code,
            ],
            dtype=torch.float32,
        ),
        "label_values": torch.full((n_residues, metric_count), float("nan")),
        "label_variances": torch.zeros((n_residues, metric_count)),
        "label_mask": torch.zeros((n_residues, metric_count), dtype=torch.bool),
        "metric_names": METRIC_NAMES,
    }


def load_model(checkpoint_path: str, device_name: str = "auto") -> tuple:
    device = resolve_device(device_name)
    checkpoint = torch.load(checkpoint_path, map_location="cpu", weights_only=True)
    model_config = ModelConfig(**checkpoint["config"]["model"])
    model = build_information_ablation(model_config, model_config.information_ablation)
    model.load_state_dict(checkpoint["model"])
    model.to(device).eval()
    return model, model_config, device


def _state_rows(
    structure: ProteinStructure,
    indices: Sequence[int],
    mean: torch.Tensor,
    log_variance: torch.Tensor,
    metric_names: Sequence[str],
) -> List[Dict[str, Any]]:
    rows = []
    for output_index, residue_index in enumerate(indices):
        row: Dict[str, Any] = {
            "residue_index": int(residue_index),
            "chain": structure.chain_ids[residue_index],
            "author_residue_id": structure.author_residue_ids[residue_index],
            "amino_acid": structure.sequence[residue_index],
        }
        for metric_index, name in enumerate(metric_names):
            row[name] = float(mean[output_index, metric_index])
            row[name + "_std"] = float(torch.exp(0.5 * log_variance[output_index, metric_index]))
        rows.append(row)
    return rows


@torch.no_grad()
def predict(
    checkpoint_path: str,
    apo_pdb: str,
    holo_pdb: Optional[str] = None,
    target_chains: Optional[Sequence[str]] = None,
    temperature: float = 300.0,
    window_ns: float = 100.0,
    device_name: str = "auto",
) -> Dict[str, Any]:
    model, config, device = load_model(checkpoint_path, device_name)
    apo_structure = parse_pdb(apo_pdb)
    if holo_pdb is None:
        chains = list(target_chains or apo_structure.unique_chains())
        batch = collate_states(
            [_structure_example(apo_pdb, chains, temperature, window_ns, 0.0)]
        )
        batch = move_to_device(batch, device)
        output = model.forward_state(batch)
        indices = torch.nonzero(
            batch["target_mask"][0] & output["mask"][0], as_tuple=False
        ).squeeze(-1)
        return {
            "mode": "state",
            "metrics": config.metrics,
            "rows": _state_rows(
                apo_structure,
                indices.cpu().tolist(),
                output["mean"][0, indices].cpu(),
                output["log_variance"][0, indices].cpu(),
                config.metrics,
            ),
        }

    holo_structure = parse_pdb(holo_pdb)
    apo_indices, holo_indices, matched_chains, identity, coverage = match_apo_to_holo(
        apo_structure, holo_structure
    )
    if target_chains:
        allowed = set(target_chains)
        keep = [
            position for position, holo_index in enumerate(holo_indices)
            if holo_structure.chain_ids[holo_index] in allowed
        ]
        apo_indices = [apo_indices[position] for position in keep]
        holo_indices = [holo_indices[position] for position in keep]
        matched_chains = sorted(allowed)
    pair_example = {
        "pair_id": "inference_pair",
        "apo": _structure_example(
            apo_pdb, apo_structure.unique_chains(), temperature, window_ns, 1.0
        ),
        "holo": _structure_example(
            holo_pdb, matched_chains, temperature, window_ns, 1.0
        ),
        "apo_indices": torch.tensor(apo_indices, dtype=torch.long),
        "holo_indices": torch.tensor(holo_indices, dtype=torch.long),
        "sequence_identity": identity,
        "coverage": coverage,
        "metadata": {
            "apo_pdb": apo_pdb,
            "holo_pdb": holo_pdb,
            "matched_holo_chains": matched_chains,
        },
    }
    batch = move_to_device(collate_pairs([pair_example]), device)
    output = model.forward_pair(batch)
    rows = []
    apo_mean = output["apo"]["mean"][0, batch["apo_indices"][0]].cpu()
    holo_mean = output["holo"]["mean"][0, batch["holo_indices"][0]].cpu()
    delta_mean = output["delta_mean"][0].cpu()
    delta_log_variance = output["delta_log_variance"][0].cpu()
    state_difference = output["state_difference"][0].cpu()
    valid_pair = output["pair_mask"][0].cpu()
    for position, (apo_index, holo_index) in enumerate(zip(apo_indices, holo_indices)):
        if not bool(valid_pair[position]):
            continue
        row: Dict[str, Any] = {
            "apo_residue_index": int(apo_index),
            "holo_residue_index": int(holo_index),
            "chain": holo_structure.chain_ids[holo_index],
            "author_residue_id": holo_structure.author_residue_ids[holo_index],
            "amino_acid": holo_structure.sequence[holo_index],
        }
        for metric_index, name in enumerate(config.metrics):
            row["apo_" + name] = float(apo_mean[position, metric_index])
            row["holo_" + name] = float(holo_mean[position, metric_index])
            row["delta_" + name] = float(delta_mean[position, metric_index])
            row["delta_" + name + "_std"] = float(
                torch.exp(0.5 * delta_log_variance[position, metric_index])
            )
            row["state_subtraction_" + name] = float(
                state_difference[position, metric_index]
            )
        rows.append(row)
    return {
        "mode": "paired_response",
        "metrics": config.metrics,
        "sequence_identity": identity,
        "coverage": coverage,
        "matched_holo_chains": matched_chains,
        "rows": rows,
    }


def _full_model_paths() -> List[Path]:
    try:
        models = Path(str(files("flexppi_weights"))) / "full"
    except ModuleNotFoundError:
        models = Path(__file__).resolve().parents[2] / "models" / "full"
    return [models / f"fold_{fold}.pt" for fold in range(5)]


@torch.inference_mode()
def predict_structure(
    pdb: str,
    checkpoint_path: Optional[str] = None,
    target_chains: Optional[Sequence[str]] = None,
    device_name: str = "auto",
) -> List[Dict[str, Any]]:
    """Predict per-residue RMSF (angstroms) and side-chain entropy (nats).

    The default is the arithmetic mean of the five system-split full models.
    The network sees all chains, independently of the selected output chains.
    """
    structure = parse_pdb(pdb)
    chains = list(target_chains or structure.unique_chains())
    missing = set(chains) - set(structure.chain_ids)
    if missing:
        raise ValueError(f"Unknown PDB chains: {', '.join(sorted(missing))}")
    example = _structure_example(pdb, chains, 300.0, 100.0, 0.0, structure)
    device = resolve_device(device_name)
    batch = move_to_device(collate_states([example]), device)
    selected = batch["target_mask"][0] & batch["geometry_mask"][0]
    if not selected.any():
        raise ValueError("Selected chains have no residues with a complete N, CA, C frame")
    paths = [Path(checkpoint_path)] if checkpoint_path else _full_model_paths()
    member_means = []
    for path in paths:
        model, config, _ = load_model(str(path), device_name)
        if config.metrics != ["rmsf", "chi_entropy"] or config.use_condition_features:
            raise ValueError("Structure prediction requires a two-target, structure-only checkpoint")
        output = model.forward_state(batch)
        member_means.append(output["mean"][0, selected].double().cpu())
    means = torch.stack(member_means).mean(dim=0).tolist()
    indices = torch.nonzero(selected).flatten().tolist()
    rows = []
    for residue, (rmsf, entropy) in zip(indices, means):
        amino_acid = structure.sequence[residue]
        rows.append({
            "chain": structure.chain_ids[residue],
            "author_residue_id": structure.author_residue_ids[residue],
            "amino_acid": amino_acid,
            "rmsf_angstrom": rmsf,
            "chi_entropy_nats": entropy if amino_acid not in {"A", "G", "X"} else None,
        })
    return rows
