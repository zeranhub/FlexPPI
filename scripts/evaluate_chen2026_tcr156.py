"""Chen 2026 (Zenodo 17065125) TCR156 external validation for FlexPPI.

Stages
------
extract : per-replica chain-fitted C-alpha RMSF and side-chain chi entropy from
          the released WT / S32M Amber trajectories (frames at 10 ns; frames
          7..110 are treated as production, matching the archive README).
predict : frozen epoch-40 ATLAS-Hchi v1 five-fold models applied to the
          starting structure of each condition (system_renum.pdb).
compare : per-chain and pooled agreement between predictions and MD labels,
          five-fold variation, retrieval of the most flexible 20 percent, and
          the WT-vs-S32M predicted difference versus the MD difference.

The system is a minimal construct (TCR V-alpha/V-beta, HLA-A*02:01 alpha1/alpha2
platform, PAP22 peptide; no beta-2-microglobulin, no alpha3 domain).  The test
is combination-level: the exact PDB entries 9NMU-9NMY never occur in training,
but the TCR alpha chain shares ~91 percent identity with the A6 TCR alpha chain
(1AO7) that is present in the ITrace training records.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from flexppi.constants import CHI_ATOMS, METRIC_NAMES, THREE_TO_ONE  # noqa: E402
from flexppi.data.dataset import collate_states  # noqa: E402
from flexppi.data.labels import (  # noqa: E402
    _chainwise_rmsf,
    _dihedral_degrees,
    _entropy,
    _kabsch_fit,
    _summarise_samples,
    save_label_file,
)
from flexppi.data.structures import parse_pdb  # noqa: E402
from flexppi.evaluation import (  # noqa: E402
    _average_precision,
    _binary_auroc,
    _correlation,
    _spearman,
)
from flexppi.inference import _structure_example  # noqa: E402
from flexppi.model.network import FlexPPI  # noqa: E402
from flexppi.trainer import move_to_device  # noqa: E402
from flexppi.config import ModelConfig  # noqa: E402

CONDITIONS = ("WT", "S32M")
CHECKPOINT_DIRECTORY = "artifacts/performance_first_v2_atlas_hchi_v1"
PRODUCTION_START_FRAME = 7
PRODUCTION_STOP_FRAME = 111
PRODUCTION_FRAME_NS = 10.0
DEAD_ZONE_RMSF = 0.10


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def git_revision() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=str(ROOT), text=True
        ).strip()
    except Exception:
        return "unknown"


def condition_directory(data_root: Path, condition: str) -> Path:
    return data_root / condition / condition


def load_protein_atoms(universe) -> Tuple[List[int], List[int], List[str], List[str], List[str], List[List[int]], List[List[int]]]:
    """Return CA indices, chain index per residue, chain ids, residue keys, names and chi quadruples."""
    ca_indices: List[int] = []
    residue_chain: List[str] = []
    residue_keys: List[str] = []
    residue_names: List[str] = []
    chi_quadruples: List[List[int]] = []
    residue_chi_indices: List[List[int]] = []

    for residue in universe.residues:
        name = str(residue.resname).upper()
        if name not in THREE_TO_ONE:
            continue
        atoms = {str(atom.name).upper(): int(atom.index) for atom in residue.atoms}
        if "CA" not in atoms:
            continue
        residue_chain.append(str(residue.segid))
        residue_keys.append("{}_{}".format(residue.segid, str(residue.resid)))
        residue_names.append(name)
        ca_indices.append(atoms["CA"])
        local: List[int] = []
        for definition in CHI_ATOMS.get(name, []):
            if all(atom_name in atoms for atom_name in definition):
                local.append(len(chi_quadruples))
                chi_quadruples.append([atoms[atom_name] for atom_name in definition])
            else:
                local = []
                break
        residue_chi_indices.append(local if local else [])

    chain_order = list(dict.fromkeys(residue_chain))
    chain_to_index = {chain: index for index, chain in enumerate(chain_order)}
    chain_index = [chain_to_index[chain] for chain in residue_chain]
    return ca_indices, chain_index, chain_order, residue_keys, residue_names, chi_quadruples, residue_chi_indices


def extract_condition(data_root: Path, condition: str, output_dir: Path) -> Dict[str, object]:
    import MDAnalysis as mda

    directory = condition_directory(data_root, condition)
    pdb_path = directory / "system_renum.pdb"
    prmtop_path = directory / "system.prmtop"
    replica_paths = sorted(directory.glob("{}_reimaged_skip50_r*.dcd".format(condition)),
                           key=lambda p: int(p.stem.split("_r")[-1]))
    if not replica_paths:
        raise FileNotFoundError("No trajectories under {}".format(directory))

    structure = parse_pdb(str(pdb_path))
    n_residues = structure.num_residues
    samples: List[np.ndarray] = []
    replica_rows: List[pd.DataFrame] = []
    replica_rmsf: List[np.ndarray] = []
    replica_manifest: List[Dict[str, object]] = []
    expected_chains = structure.unique_chains()

    for replica_path in replica_paths:
        universe = mda.Universe(str(pdb_path), str(replica_path))
        (ca_indices, chain_index, chain_order, residue_keys, residue_names,
         chi_quadruples, residue_chi_indices) = load_protein_atoms(universe)
        if len(ca_indices) != n_residues:
            raise ValueError("Residue count mismatch in {}: {} vs {}".format(
                replica_path.name, len(ca_indices), n_residues))
        parsed_keys = ["{}_{}".format(c, r) for c, r in
                       zip(structure.chain_ids, structure.author_residue_ids)]
        if residue_keys != parsed_keys:
            raise ValueError("Residue identity mismatch in {} (first difference: {} vs {})".format(
                replica_path.name,
                next((k for k, p in zip(residue_keys, parsed_keys) if k != p), None),
                next((p for k, p in zip(residue_keys, parsed_keys) if k != p), None)))
        if sorted(chain_order) != sorted(expected_chains):
            raise ValueError("Chain mismatch in {}: {} vs {}".format(
                replica_path.name, chain_order, expected_chains))

        n_frames = universe.trajectory.n_frames
        stop = min(PRODUCTION_STOP_FRAME, n_frames)
        if stop - PRODUCTION_START_FRAME < 2:
            raise ValueError("Too few production frames in {}".format(replica_path.name))
        frames = list(range(PRODUCTION_START_FRAME, stop))

        ca_frames: List[np.ndarray] = []
        chi_codes: List[np.ndarray] = []
        quad_array = np.asarray(chi_quadruples, dtype=np.int64)
        for timestep in universe.trajectory[PRODUCTION_START_FRAME:stop]:
            positions = timestep.positions
            ca_frames.append(np.asarray(positions[ca_indices], dtype=np.float32).copy())
            if len(quad_array):
                angles = _dihedral_degrees(positions[quad_array])
                bins = np.clip(np.floor((angles + 180.0) / 120.0).astype(np.int8), 0, 2)
                chi_codes.append(bins)

        ca_array = np.stack(ca_frames, axis=0)
        rmsf = _chainwise_rmsf(ca_array, np.asarray(chain_index))
        sample = np.full((n_residues, len(METRIC_NAMES)), np.nan, dtype=np.float32)
        sample[:, 0] = rmsf
        if chi_codes:
            chi_array = np.stack(chi_codes, axis=0)
            for residue_index, local in enumerate(residue_chi_indices):
                if not local:
                    continue
                codes = np.zeros(len(frames), dtype=np.int32)
                for power, chi_index in enumerate(local):
                    codes += chi_array[:, chi_index].astype(np.int32) * (3 ** power)
                sample[residue_index, 1] = _entropy(codes)
        samples.append(sample)
        replica_rmsf.append(rmsf)
        replica_rows.append(pd.DataFrame({
            "condition": condition, "replica": replica_path.name,
            "chain": [structure.chain_ids[i] for i in range(n_residues)],
            "author_residue_id": structure.author_residue_ids,
            "amino_acid": list(structure.sequence),
            "rmsf": rmsf, "chi_entropy": sample[:, 1]}))
        replica_manifest.append({
            "replica": replica_path.name, "frames_total": int(n_frames),
            "frames_used": len(frames), "production_start": PRODUCTION_START_FRAME,
            "production_stop": int(stop) - 1, "sha256": sha256(replica_path),
            "size_bytes": replica_path.stat().st_size})
        print("{} {}: {} frames used (of {})".format(
            condition, replica_path.name, len(frames), n_frames), flush=True)

    values, variances, mask = _summarise_samples(samples, n_residues)
    output_dir.mkdir(parents=True, exist_ok=True)
    label_path = output_dir / "{}_labels.npz".format(condition)
    save_label_file(str(label_path), values, variances, mask)

    aggregate = pd.DataFrame({
        "condition": condition,
        "chain": [structure.chain_ids[i] for i in range(n_residues)],
        "author_residue_id": structure.author_residue_ids,
        "amino_acid": list(structure.sequence),
        "rmsf": values[:, 0], "rmsf_variance": variances[:, 0], "rmsf_mask": mask[:, 0],
        "chi_entropy": values[:, 1], "chi_entropy_variance": variances[:, 1],
        "chi_entropy_mask": mask[:, 1]})
    aggregate.to_csv(output_dir / "{}_labels.csv".format(condition), index=False)
    pd.concat(replica_rows, ignore_index=True).to_csv(
        output_dir / "{}_per_replica.csv".format(condition), index=False)
    np.savez_compressed(output_dir / "{}_replica_rmsf.npz".format(condition),
                        rmsf=np.stack(replica_rmsf, axis=0),
                        replicas=np.asarray([p.name for p in replica_paths]))

    manifest = {
        "condition": condition, "structure": str(pdb_path),
        "structure_sha256": sha256(pdb_path), "prmtop": str(prmtop_path),
        "prmtop_sha256": sha256(prmtop_path), "residues": int(n_residues),
        "chains": {chain: int(sum(1 for value in structure.chain_ids if value == chain))
                   for chain in structure.unique_chains()},
        "production_start_frame": PRODUCTION_START_FRAME,
        "production_stop_frame_inclusive": PRODUCTION_STOP_FRAME - 1,
        "frame_spacing_ns": PRODUCTION_FRAME_NS,
        "replicas": replica_manifest, "n_replicas": len(replica_manifest),
        "label_definition": "mean and sample variance across replicas; per-replica "
                            "RMSF is chain-fitted C-alpha fluctuation over production frames",
    }
    (output_dir / "{}_label_manifest.json".format(condition)).write_text(
        json.dumps(manifest, indent=2))
    print("{}: wrote labels for {} residues ({} replicas)".format(
        condition, n_residues, len(replica_manifest)), flush=True)
    return manifest


def extract(data_root: Path, output_dir: Path) -> None:
    for condition in CONDITIONS:
        extract_condition(data_root, condition, output_dir)


def load_checkpoint_model(fold: int, device_name: str = "cpu") -> Tuple[FlexPPI, object]:
    import torch

    path = ROOT / CHECKPOINT_DIRECTORY / "fold_{}".format(fold) / "epoch_040.pt"
    checkpoint = torch.load(str(path), map_location="cpu", weights_only=False)
    if checkpoint.get("epoch") != 40:
        raise ValueError("{} is not the frozen epoch-40 checkpoint".format(path))
    config = ModelConfig(**checkpoint["config"]["model"])
    if list(config.metrics) != ["rmsf", "chi_entropy"]:
        raise ValueError("Unexpected metric set in {}: {}".format(path, config.metrics))
    device = torch.device(device_name)
    model = FlexPPI(config)
    model.load_state_dict(checkpoint["model"])
    model.to(device).eval()
    return model, config


def predict(data_root: Path, output_dir: Path, device_name: str = "cpu") -> None:
    import torch

    torch.set_num_threads(6)
    output_dir.mkdir(parents=True, exist_ok=True)
    checkpoint_hashes = {}
    for fold in range(5):
        checkpoint_path = ROOT / CHECKPOINT_DIRECTORY / "fold_{}".format(fold) / "epoch_040.pt"
        checkpoint_hashes[str(fold)] = sha256(checkpoint_path)

    for condition in CONDITIONS:
        directory = condition_directory(data_root, condition)
        pdb_path = directory / "system_renum.pdb"
        structure = parse_pdb(str(pdb_path))
        example = _structure_example(str(pdb_path), structure.unique_chains(), 300.0, 100.0, 0.0)
        batch = collate_states([example])
        rows_by_fold: List[pd.DataFrame] = []
        for fold in range(5):
            model, config = load_checkpoint_model(fold, device_name)
            moved = move_to_device(batch, torch.device(device_name))
            with torch.inference_mode():
                output = model.forward_state(moved)
            mean = output["mean"][0].cpu().numpy()
            log_variance = output["log_variance"][0].cpu().numpy()
            target_mask = moved["target_mask"][0].cpu().numpy().astype(bool)
            geometry_mask = moved["geometry_mask"][0].cpu().numpy().astype(bool)
            output_mask = output["mask"][0].cpu().numpy().astype(bool)
            valid = target_mask & geometry_mask & output_mask
            metric_index = {name: index for index, name in enumerate(config.metrics)}
            frame = pd.DataFrame({
                "condition": condition, "fold": fold,
                "chain": structure.chain_ids,
                "author_residue_id": structure.author_residue_ids,
                "amino_acid": list(structure.sequence),
                "geometry_mask": geometry_mask, "output_mask": output_mask,
                "valid": valid,
                "rmsf": mean[:, metric_index["rmsf"]],
                "rmsf_std": np.exp(0.5 * log_variance[:, metric_index["rmsf"]]),
                "chi_entropy": mean[:, metric_index["chi_entropy"]],
                "chi_entropy_std": np.exp(0.5 * log_variance[:, metric_index["chi_entropy"]]),
            })
            rows_by_fold.append(frame)
            del model
        all_folds = pd.concat(rows_by_fold, ignore_index=True)
        all_folds.to_csv(output_dir / "{}_per_fold.csv".format(condition), index=False)
        summary = (all_folds[all_folds.valid]
                   .groupby(["chain", "author_residue_id", "amino_acid"], as_index=False)
                   .agg(rmsf=("rmsf", "mean"), rmsf_std_fold=("rmsf", "std"),
                        chi_entropy=("chi_entropy", "mean"),
                        chi_entropy_std_fold=("chi_entropy", "std"),
                        n_folds=("rmsf", "size")))
        summary.to_csv(output_dir / "{}_five_fold_mean.csv".format(condition), index=False)
        print("{}: predicted {} residues over 5 folds".format(
            condition, int(all_folds.valid.sum() // 5)), flush=True)

    (output_dir / "checkpoint_hashes.json").write_text(
        json.dumps(checkpoint_hashes, indent=2))


def regression_metrics(prediction: np.ndarray, target: np.ndarray) -> Dict[str, Optional[float]]:
    error = prediction - target
    return {"n": int(len(target)), "mae": float(np.abs(error).mean()),
            "rmse": float(np.sqrt(np.square(error).mean())),
            "bias": float(error.mean()),
            "pearson": _correlation(prediction, target),
            "spearman": _spearman(prediction, target)}


def descriptive(values: np.ndarray) -> Dict[str, float]:
    return {"n": int(len(values)), "mean": float(np.mean(values)),
            "std": float(np.std(values, ddof=1)) if len(values) > 1 else 0.0,
            "min": float(np.min(values)), "max": float(np.max(values))}


def replica_reliability(labels_dir: Path, condition: str, structure_chains: List[str]) -> Dict[str, object]:
    archive = np.load(labels_dir / "{}_replica_rmsf.npz".format(condition))
    matrix = archive["rmsf"]
    keys = pd.read_csv(labels_dir / "{}_labels.csv".format(condition),
                       dtype={"author_residue_id": str})
    output: Dict[str, object] = {"n_replicas": int(matrix.shape[0])}
    for chain in structure_chains:
        indices = np.flatnonzero(keys["chain"].to_numpy(str) == chain)
        values = []
        for first in range(matrix.shape[0]):
            for second in range(first + 1, matrix.shape[0]):
                value = _spearman(matrix[first, indices], matrix[second, indices])
                if value is not None:
                    values.append(value)
        output[chain] = {"mean_pairwise_spearman": float(np.mean(values)) if values else None,
                         "n_pairs": len(values)}
    return output


def compare(data_root: Path, labels_dir: Path, predictions_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    report: Dict[str, object] = {"conditions": {}, "wt_vs_s32m": {}, "provenance": {
        "labels_dir": str(labels_dir), "predictions_dir": str(predictions_dir),
        "checkpoint_directory": CHECKPOINT_DIRECTORY,
        "git_revision": git_revision(),
        "prmtop_sha256": sha256(condition_directory(data_root, "WT") / "system.prmtop"),
        "pdb_sha256": sha256(condition_directory(data_root, "WT") / "system_renum.pdb"),
        "implementation_sha256": sha256(Path(__file__))}}

    rename = {"rmsf": "pred_rmsf", "rmsf_std_fold": "pred_rmsf_fold_sd",
              "chi_entropy": "pred_chi_entropy",
              "chi_entropy_std_fold": "pred_chi_entropy_fold_sd"}
    keys = ["chain", "author_residue_id", "amino_acid"]
    merged_tables = {}
    for condition in CONDITIONS:
        labels = pd.read_csv(labels_dir / "{}_labels.csv".format(condition),
                             dtype={"author_residue_id": str})
        per_fold = pd.read_csv(predictions_dir / "{}_per_fold.csv".format(condition),
                               dtype={"author_residue_id": str})
        mean_prediction = pd.read_csv(
            predictions_dir / "{}_five_fold_mean.csv".format(condition),
            dtype={"author_residue_id": str}).rename(columns=rename)
        merged = labels.merge(mean_prediction, on=keys, how="inner", validate="one_to_one")
        merged_tables[condition] = merged
        merged.to_csv(output_dir / "{}_merged.csv".format(condition), index=False)

        condition_report: Dict[str, object] = {"chains": {}, "overall": {}}
        for metric in ("rmsf", "chi_entropy"):
            label_mask = merged[metric + "_mask"].fillna(False).to_numpy(bool)
            mask = label_mask & merged["n_folds"].to_numpy(float).__eq__(5.0)
            prediction = merged.loc[mask, "pred_" + metric].to_numpy(float)
            target = merged.loc[mask, metric].to_numpy(float)
            baseline = (merged[metric] - merged.groupby("chain")[metric].transform("mean")).abs()
            condition_report["overall"][metric] = {
                "five_fold_mean": regression_metrics(prediction, target),
                "per_chain_mean_baseline_mae": float(baseline[label_mask].mean())}

            fold_rows = []
            for fold in range(5):
                fold_table = per_fold[per_fold.fold == fold].merge(
                    labels, on=keys, suffixes=("_pred", "_label"),
                    how="inner", validate="one_to_one")
                valid = fold_table["valid"].to_numpy(bool) & fold_table[
                    metric + "_mask"].fillna(False).to_numpy(bool)
                fold_rows.append({"fold": fold, **regression_metrics(
                    fold_table.loc[valid, metric + "_pred"].to_numpy(float),
                    fold_table.loc[valid, metric + "_label"].to_numpy(float))})
            fold_frame = pd.DataFrame(fold_rows)
            fold_frame.to_csv(output_dir / "{}_fold_metrics_{}.csv".format(condition, metric),
                              index=False)
            condition_report["overall"][metric]["fold_mean"] = {
                key: float(fold_frame[key].mean()) for key in ("mae", "rmse", "pearson", "spearman")}
            condition_report["overall"][metric]["fold_sd"] = {
                key: float(fold_frame[key].std(ddof=1)) for key in ("mae", "rmse", "pearson", "spearman")}

        for chain, group in merged.groupby("chain"):
            chain_report: Dict[str, object] = {}
            for metric in ("rmsf", "chi_entropy"):
                valid = (group[metric + "_mask"].fillna(False).to_numpy(bool)
                         & group["n_folds"].to_numpy(float).__eq__(5.0))
                if valid.sum() < 5:
                    chain_report[metric] = None
                    continue
                prediction = group.loc[valid, "pred_" + metric].to_numpy(float)
                target = group.loc[valid, metric].to_numpy(float)
                entry = regression_metrics(prediction, target)
                if metric == "rmsf":
                    positive = target >= np.quantile(target, 0.8)
                    entry["top20_auroc"] = _binary_auroc(prediction, positive)
                    entry["top20_average_precision"] = _average_precision(prediction, positive)
                chain_report[metric] = entry
            condition_report["chains"][str(chain)] = chain_report

        descriptives: Dict[str, object] = {}
        for metric in ("rmsf", "chi_entropy"):
            label_mask = merged[metric + "_mask"].fillna(False).to_numpy(bool)
            prediction_mask = merged["n_folds"].to_numpy(float) == 5.0
            descriptives[metric] = {
                "labels": descriptive(merged.loc[label_mask, metric].to_numpy(float)),
                "predictions": descriptive(
                    merged.loc[prediction_mask, "pred_" + metric].to_numpy(float))}
        condition_report["descriptives"] = descriptives
        condition_report["replica_reliability"] = replica_reliability(
            labels_dir, condition, sorted(merged["chain"].unique()))

        corrected: Dict[str, object] = {}
        for metric in ("rmsf", "chi_entropy"):
            label_mask = merged[metric + "_mask"].fillna(False).to_numpy(bool)
            mask = label_mask & (merged["n_folds"].to_numpy(float) == 5.0)
            subset = merged.loc[mask, ["amino_acid", metric, "pred_" + metric]].copy()
            prior = subset.groupby("amino_acid")[metric].transform("mean")
            residual_target = subset[metric].to_numpy(float) - prior.to_numpy(float)
            residual_prediction = subset["pred_" + metric].to_numpy(float) - prior.to_numpy(float)
            within = []
            for _, group in subset.groupby("amino_acid"):
                if len(group) >= 8:
                    value = _spearman(group["pred_" + metric].to_numpy(float),
                                      group[metric].to_numpy(float))
                    if value is not None:
                        within.append(value)
            corrected[metric] = {
                "n": int(mask.sum()),
                "residual_spearman": _spearman(residual_prediction, residual_target),
                "residual_pearson": _correlation(residual_prediction, residual_target),
                "mean_within_amino_acid_spearman": float(np.mean(within)) if within else None,
                "n_amino_acid_groups": len(within)}
        condition_report["amino_acid_corrected"] = corrected
        report["conditions"][condition] = condition_report

    wt = merged_tables["WT"]; s32m = merged_tables["S32M"]
    delta = wt[keys + ["rmsf", "pred_rmsf", "rmsf_mask"]].merge(
        s32m[keys + ["rmsf", "pred_rmsf"]], on=keys, suffixes=("_wt", "_s32m"),
        validate="one_to_one")
    valid = delta["rmsf_mask"].fillna(False).to_numpy(bool)
    md_delta = delta.loc[valid, "rmsf_s32m"].to_numpy(float) - delta.loc[valid, "rmsf_wt"].to_numpy(float)
    pred_delta = delta.loc[valid, "pred_rmsf_s32m"].to_numpy(float) - delta.loc[valid, "pred_rmsf_wt"].to_numpy(float)
    active = np.abs(md_delta) > DEAD_ZONE_RMSF
    wt_vs = {
        "n_residues": int(valid.sum()),
        "md_delta_mean_abs": float(np.abs(md_delta).mean()),
        "mae_pred_delta_vs_md_delta": float(np.abs(pred_delta - md_delta).mean()),
        "mae_zero_baseline": float(np.abs(md_delta).mean()),
        "pearson": _correlation(pred_delta, md_delta),
        "spearman": _spearman(pred_delta, md_delta),
        "active_n": int(active.sum()),
        "active_sign_accuracy": float(
            (np.sign(pred_delta[active]) == np.sign(md_delta[active])).mean())
        if active.any() else None,
        "dead_zone_rmsf": DEAD_ZONE_RMSF,
    }
    delta.assign(md_delta=md_delta, pred_delta=pred_delta, active=active).to_csv(
        output_dir / "wt_vs_s32m_delta_rmsf.csv", index=False)
    report["wt_vs_s32m"]["rmsf"] = wt_vs
    report["caveats"] = [
        "Minimal construct: TCR V domains + HLA-A*02:01 alpha1/alpha2 + PAP22 peptide; "
        "no beta-2-microglobulin and no alpha3 domain.",
        "Trajectory frames are 10 ns apart; only slow modes are represented and "
        "side-chain chi entropy is not converged.",
        "The authors report a water-oxygen Lennard-Jones radius assignment issue; "
        "within-dataset comparisons are the appropriate use.",
        "Model input is the MD starting structure (system_renum.pdb); no crystal "
        "structure and no MD frame was selected using prediction agreement.",
        "Combination-level novelty only: 9NMU-9NMY are absent from training, but the "
        "TCR alpha chain is ~91 percent identical to the A6 TCR alpha chain (1AO7).",
    ]
    (output_dir / "report.json").write_text(json.dumps(report, indent=2))
    write_summary_markdown(report, output_dir / "SUMMARY.md")
    print(json.dumps(report["wt_vs_s32m"], indent=2), flush=True)
    print("wrote", output_dir / "report.json")


def write_summary_markdown(report: Dict[str, object], path: Path) -> None:
    lines = ["# Chen 2026 TCR156 / FlexPPI external validation", ""]
    lines.append("Model: ATLAS-Hchi v1 five-fold epoch-40 checkpoints. "
                 "Reference: chain-fitted C-alpha RMSF from 10 WT / 10 S32M "
                 "replicas, production frames 7-110 at 10 ns spacing.")
    lines.append("")
    for condition in CONDITIONS:
        block = report["conditions"][condition]
        lines.append("## {}".format(condition))
        lines.append("")
        lines.append("| metric | n | pooled Spearman | fold SD | pooled Pearson | MAE | fold SD |")
        lines.append("|---|---|---|---|---|---|---|")
        for metric in ("rmsf", "chi_entropy"):
            overall = block["overall"][metric]
            five = overall["five_fold_mean"]
            lines.append("| {} | {} | {:.3f} | {:.3f} | {:.3f} | {:.3f} | {:.3f} |".format(
                metric, five["n"], five["spearman"], overall["fold_sd"]["spearman"],
                five["pearson"], five["mae"], overall["fold_sd"]["mae"]))
        lines.append("")
        lines.append("| chain | n | Spearman | MAE | top-20% AUROC | top-20% AP |")
        lines.append("|---|---|---|---|---|---|")
        for chain, values in block["chains"].items():
            entry = values["rmsf"]
            if entry is None:
                continue
            lines.append("| {} | {} | {:.3f} | {:.3f} | {} | {} |".format(
                chain, entry["n"], entry["spearman"], entry["mae"],
                "n/a" if entry.get("top20_auroc") is None else "{:.3f}".format(entry["top20_auroc"]),
                "n/a" if entry.get("top20_average_precision") is None else
                "{:.3f}".format(entry["top20_average_precision"])))
        lines.append("")
        lines.append("Top-20% retrieval chance level: AUROC 0.5, AP 0.2.")
        reliability = block["replica_reliability"]
        lines.append("")
        lines.append("Per-replica ranking reliability (mean pairwise Spearman): " + ", ".join(
            "{} {:.3f}".format(chain, value["mean_pairwise_spearman"])
            for chain, value in reliability.items() if isinstance(value, dict)))
        corrected = block.get("amino_acid_corrected", {})
        if corrected:
            lines.append("")
            lines.append("Amino-acid-corrected (type mean removed from both sides): " + "; ".join(
                "{} residual Spearman {} (within-type mean {})".format(
                    metric,
                    "n/a" if values["residual_spearman"] is None else
                    "{:.3f}".format(values["residual_spearman"]),
                    "n/a" if values["mean_within_amino_acid_spearman"] is None else
                    "{:.3f}".format(values["mean_within_amino_acid_spearman"]))
                for metric, values in corrected.items()))
        lines.append("")
    delta = report["wt_vs_s32m"]["rmsf"]
    lines.append("## WT vs S32M predicted difference")
    lines.append("")
    lines.append("- active residues (|MD delta| > {:.2f} A): {} of {}".format(
        delta["dead_zone_rmsf"], delta["active_n"], delta["n_residues"]))
    lines.append("- MD mean absolute delta: {:.3f} A; model delta MAE vs zero baseline: "
                 "{:.3f} vs {:.3f} A".format(
                     delta["md_delta_mean_abs"], delta["mae_pred_delta_vs_md_delta"],
                     delta["mae_zero_baseline"]))
    lines.append("- Spearman/Pearson of predicted vs MD delta: {} / {}".format(
        "n/a" if delta["spearman"] is None else "{:.3f}".format(delta["spearman"]),
        "n/a" if delta["pearson"] is None else "{:.3f}".format(delta["pearson"])))
    lines.append("- active sign accuracy: {}".format(
        "n/a" if delta["active_sign_accuracy"] is None else
        "{:.3f}".format(delta["active_sign_accuracy"])))
    lines.append("")
    lines.append("## Caveats")
    lines.append("")
    for caveat in report["caveats"]:
        lines.append("- " + caveat)
    lines.append("")
    lines.append("Absolute amplitudes are not comparable: the released frames are 10 ns "
                 "apart, so the MD reference contains only slow modes (per-chain RMSF "
                 "spread of ~0.05 A), whereas the model was trained on finer-sampled "
                 "100 ns windows with a wider absolute range. Within-chain ranking and "
                 "top-20% retrieval are the interpretable readouts.")
    path.write_text("\n".join(lines) + "\n")


def assign_regions(chain: str, sequence: str) -> List[str]:
    """Operational region labels for the minimal construct.

    The released construct has the two conserved V-domain cysteines removed
    (author numbers 22-23 are absent), so the anchors are the conserved W of the
    W-YRQ motif, the following N-G pair, the S-A-[T/L]-Y framework motif and the
    J-region F-G-x-G motif:
      CDR1 : W - 11 through W - 3 (alpha: TYSDRGSQS, author 24-32, contains S32).
      CDR2 : N-G - 4 through N-G + 2.
      CDR3 : two residues after the Y of SA[T/L]Y through the J-Phe minus one.
    Chain C is the PAP22 peptide; chain A is the HLA-A*02:01 alpha1/alpha2 platform.
    """
    if chain not in ("D", "E"):
        return ["peptide" if chain == "C" else "HLA"] * len(sequence)
    labels = ["framework"] * len(sequence)
    w_motif = re.search(r"WY[RQ]", sequence)
    ng = sequence.find("NG", w_motif.start() + 1)
    j_motif = re.search(r"FG.G", sequence)
    sa_motif = re.search(r"SA[TL]Y", sequence)
    if w_motif is None or ng < 0 or j_motif is None or sa_motif is None:
        raise ValueError("Region anchors not found for chain {}".format(chain))
    tryptophan = w_motif.start()
    j_phe = j_motif.start()
    cdr3_start = sa_motif.start() + 5
    spans = (("CDR1", tryptophan - 11, tryptophan - 3),
             ("CDR2", ng - 4, ng + 2),
             ("CDR3", cdr3_start, j_phe - 1))
    for name, start, stop in spans:
        for index in range(max(start, 0), min(stop, len(sequence) - 1) + 1):
            labels[index] = name
    return labels


def start_structure_difference(data_root: Path, output_dir: Path) -> Dict[str, object]:
    wt = parse_pdb(str(condition_directory(data_root, "WT") / "system_renum.pdb"))
    mutant = parse_pdb(str(condition_directory(data_root, "S32M") / "system_renum.pdb"))
    wt_ca = {("{}_{}".format(c, r)): i for i, (c, r) in enumerate(
        zip(wt.chain_ids, wt.author_residue_ids)) if wt.backbone_mask[i, 1]}
    mutant_ca = {("{}_{}".format(c, r)): i for i, (c, r) in enumerate(
        zip(mutant.chain_ids, mutant.author_residue_ids)) if mutant.backbone_mask[i, 1]}
    keys = sorted(set(wt_ca) & set(mutant_ca),
                  key=lambda key: (wt.chain_index[wt_ca[key]], wt_ca[key]))
    wt_xyz = np.asarray([wt.backbone[wt_ca[key], 1] for key in keys], dtype=np.float64)
    mutant_xyz = np.asarray([mutant.backbone[mutant_ca[key], 1] for key in keys], dtype=np.float64)
    fitted = _kabsch_fit(mutant_xyz, wt_xyz)
    displacement = np.linalg.norm(fitted - wt_xyz, axis=1)
    rows = pd.DataFrame({"key": keys, "displacement_angstrom": displacement.astype(float)})
    rows["chain"] = [key.split("_")[0] for key in keys]
    rows["author_residue_id"] = [key.split("_", 1)[1] for key in keys]
    rows.to_csv(output_dir / "start_structure_difference.csv", index=False)
    report = {
        "n_common_ca": len(keys),
        "mean_displacement": float(displacement.mean()),
        "max_displacement": float(displacement.max()),
        "n_above_0.5A": int((displacement > 0.5).sum()),
        "residues_above_0.5A": rows.loc[displacement > 0.5, ["chain", "author_residue_id"]]
        .astype(str).to_dict(orient="records"),
        "note": "WT and S32M starting coordinates differ only by the engineered "
                "single point mutation and local relaxation; no MD frame was fitted.",
    }
    (output_dir / "start_structure_difference.json").write_text(json.dumps(report, indent=2))
    return report


def residues(data_root: Path, labels_dir: Path, predictions_dir: Path, output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    definitions: Dict[str, object] = {}
    merged_conditions: Dict[str, pd.DataFrame] = {}
    region_rows = []
    for condition in CONDITIONS:
        labels = pd.read_csv(labels_dir / "{}_labels.csv".format(condition),
                             dtype={"author_residue_id": str})
        per_replica = pd.read_csv(labels_dir / "{}_per_replica.csv".format(condition),
                                  dtype={"author_residue_id": str})
        prediction = pd.read_csv(predictions_dir / "{}_five_fold_mean.csv".format(condition),
                                 dtype={"author_residue_id": str})
        per_fold = pd.read_csv(predictions_dir / "{}_per_fold.csv".format(condition),
                               dtype={"author_residue_id": str})

        replica_stats = (per_replica.groupby(["chain", "author_residue_id"], as_index=False)
                         .agg(md_rmsf=("rmsf", "mean"), md_rmsf_sd=("rmsf", "std"),
                              md_rmsf_sem=("rmsf", "sem"), n_replicas=("rmsf", "size"),
                              md_chi=("chi_entropy", "mean"),
                              md_chi_sd=("chi_entropy", "std")))
        fold_sd = (per_fold.groupby(["chain", "author_residue_id"], as_index=False)
                   .agg(pred_rmsf_fold_sd=("rmsf", "std"),
                        pred_chi_fold_sd=("chi_entropy", "std")))
        merged = (labels.merge(replica_stats, on=["chain", "author_residue_id"], how="left",
                               validate="one_to_one")
                  .merge(prediction.rename(columns={"rmsf": "pred_rmsf",
                                                    "chi_entropy": "pred_chi_entropy"})
                         [["chain", "author_residue_id", "pred_rmsf", "pred_chi_entropy"]],
                         on=["chain", "author_residue_id"], how="left", validate="one_to_one")
                  .merge(fold_sd, on=["chain", "author_residue_id"], how="left",
                         validate="one_to_one"))
        merged["region"] = ""
        for chain, group in merged.groupby("chain", sort=False):
            sequence = "".join(group["amino_acid"].astype(str))
            merged.loc[group.index, "region"] = assign_regions(chain, sequence)
        merged["rmsf_error"] = merged["pred_rmsf"] - merged["rmsf"]
        merged["cond"] = condition
        merged_conditions[condition] = merged
        merged.to_csv(output_dir / "per_residue_comparison_{}.csv".format(condition), index=False)

        for (chain, region), group in merged.groupby(["chain", "region"], sort=False):
            for metric, prediction_column in (("rmsf", "pred_rmsf"), ("chi_entropy", "pred_chi_entropy")):
                target_column = metric
                mask = group[metric + "_mask"].fillna(False).to_numpy(bool)
                mask &= group[prediction_column].notna().to_numpy()
                if metric == "chi_entropy":
                    mask &= group["md_chi"].notna().to_numpy()
                    target = group.loc[mask, "md_chi"].to_numpy(float)
                else:
                    target = group.loc[mask, target_column].to_numpy(float)
                if mask.sum() == 0:
                    continue
                prediction_values = group.loc[mask, prediction_column].to_numpy(float)
                entry = {"condition": condition, "chain": chain, "region": region,
                         "metric": metric, "n": int(mask.sum()),
                         "md_mean": float(target.mean()), "md_sd": float(target.std(ddof=1))
                         if mask.sum() > 1 else 0.0,
                         "pred_mean": float(prediction_values.mean()),
                         "pred_sd": float(prediction_values.std(ddof=1))
                         if mask.sum() > 1 else 0.0,
                         "bias": float((prediction_values - target).mean()),
                         "mae": float(np.abs(prediction_values - target).mean()),
                         "pearson": _correlation(prediction_values, target) if mask.sum() >= 6 else None,
                         "spearman": _spearman(prediction_values, target) if mask.sum() >= 6 else None}
                region_rows.append(entry)
        definitions[condition] = {chain: sorted(set(group["region"]))
                                  for chain, group in merged.groupby("chain")}
    pd.DataFrame(region_rows).to_csv(output_dir / "region_metrics.csv", index=False)

    wt = merged_conditions["WT"]; mutant = merged_conditions["S32M"]
    keys = ["chain", "author_residue_id"]
    delta = wt[keys + ["region", "amino_acid", "rmsf", "md_rmsf", "md_rmsf_sd",
                       "pred_rmsf", "pred_rmsf_fold_sd"]].merge(
        mutant[keys + ["rmsf", "md_rmsf", "md_rmsf_sd", "pred_rmsf"]], on=keys,
        suffixes=("_wt", "_s32m"), validate="one_to_one")
    delta["md_delta"] = delta["md_rmsf_s32m"] - delta["md_rmsf_wt"]
    delta["md_delta_sd"] = np.sqrt(delta["md_rmsf_sd_wt"] ** 2 + delta["md_rmsf_sd_s32m"] ** 2)
    delta["pred_delta"] = delta["pred_rmsf_s32m"] - delta["pred_rmsf_wt"]
    delta.to_csv(output_dir / "wt_vs_s32m_per_residue.csv", index=False)

    region_delta = []
    for (chain, region), group in delta.groupby(["chain", "region"], sort=False):
        active = group["md_delta"].abs() > DEAD_ZONE_RMSF
        region_delta.append({
            "chain": chain, "region": region, "n": len(group),
            "md_delta_mean": float(group["md_delta"].mean()),
            "md_delta_abs_mean": float(group["md_delta"].abs().mean()),
            "pred_delta_mean": float(group["pred_delta"].mean()),
            "n_active": int(active.sum()),
            "pearson": _correlation(group["pred_delta"].to_numpy(float),
                                    group["md_delta"].to_numpy(float)) if len(group) >= 6 else None})
    pd.DataFrame(region_delta).to_csv(output_dir / "wt_vs_s32m_region_delta.csv", index=False)

    start_report = start_structure_difference(data_root, output_dir)
    summary = {"region_definitions": definitions,
               "start_structure_difference": start_report}
    (output_dir / "regions_summary.json").write_text(json.dumps(summary, indent=2))
    write_regions_markdown(output_dir, region_rows, region_delta, start_report)
    print("per-residue tables and region metrics written", flush=True)


def write_regions_markdown(output_dir: Path, region_rows: List[dict],
                           region_delta: List[dict], start_report: Dict[str, object]) -> None:
    lines = ["# WT / S32M relationship and per-region analysis", ""]
    lines.append("## What WT and S32M are")
    lines.append("")
    lines.append("S32M is the engineered catch-bond variant of TCR156: a single point "
                 "mutation of the TCR alpha-chain CDR1 serine 32 to methionine, selected "
                 "in the source study for strongly increased TCR-pMHC bond lifetime and "
                 "T-cell potency without a change in 3D affinity. The two released "
                 "datasets are independent MD replica sets of otherwise identical "
                 "constructs (10 replicas each), so all observed differences should "
                 "trace back to that single site, plus construction/relaxation noise.")
    lines.append("")
    lines.append("Starting-structure comparison (all 401 common C-alpha, WT superposed on "
                 "S32M): mean displacement {:.3f} A, max {:.3f} A, {} residues above 0.5 A."
                 .format(start_report["mean_displacement"], start_report["max_displacement"],
                         start_report["n_above_0.5A"]))
    lines.append("The starting coordinates therefore differ beyond alpha S32M itself "
                 "(mostly in the HLA platform); the model's predicted WT-vs-S32M delta "
                 "reflects the supplied structures, not mutation knowledge.")
    lines.append("")
    lines.append("## Region definitions (operational, author numbering)")
    lines.append("")
    lines.append("CDR1 = W-11..W-3 of the W-Y[RQ] motif; CDR2 = N-G-4..N-G+2; CDR3 = two "
                 "residues after the Y of SA[T/L]Y through the J-Phe minus one. Chain A "
                 "= HLA platform, chain C = PAP22 peptide. The construct lacks the two "
                 "conserved V-domain cysteines, so classical Cys anchors are unavailable.")
    lines.append("")
    for condition in CONDITIONS:
        lines.append("## {} region metrics (RMSF)".format(condition))
        lines.append("")
        lines.append("| chain | region | n | MD mean (A) | MD SD | model mean (A) | model SD | MAE | Spearman |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        for row in region_rows:
            if row["condition"] != condition or row["metric"] != "rmsf":
                continue
            lines.append("| {} | {} | {} | {:.3f} | {:.3f} | {:.3f} | {:.3f} | {:.3f} | {} |".format(
                row["chain"], row["region"], row["n"], row["md_mean"], row["md_sd"],
                row["pred_mean"], row["pred_sd"], row["mae"],
                "n/a" if row["spearman"] is None else "{:.3f}".format(row["spearman"])))
        lines.append("")
    lines.append("## WT versus S32M by region (RMSF)")
    lines.append("")
    lines.append("| chain | region | n | MD delta mean (A) | |MD delta| mean | model delta mean (A) | active |")
    lines.append("|---|---|---|---|---|---|---|")
    for row in region_delta:
        lines.append("| {} | {} | {} | {:+.3f} | {:.3f} | {:+.3f} | {} |".format(
            row["chain"], row["region"], row["n"], row["md_delta_mean"],
            row["md_delta_abs_mean"], row["pred_delta_mean"], row["n_active"]))
    lines.append("")
    lines.append("Region-level Spearman with n=7-11 residues is indicative only; the WT-vs-S32M "
                 "MD differences are dominated by replica noise at this 10 ns frame spacing.")
    (output_dir / "REGIONS.md").write_text("\n".join(lines) + "\n")


def figures(labels_dir: Path, predictions_dir: Path, output_dir: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=True)
    region_color = {"CDR1": "#f6c85f", "CDR2": "#9bd4a0", "CDR3": "#ef8a8a"}
    chain_titles = {"A": "HLA-A*02:01 alpha1/alpha2 (chain A)",
                    "C": "PAP22 peptide (chain C)",
                    "D": "TCR156 V-alpha (chain D)",
                    "E": "TCR156 V-beta (chain E)"}

    for condition in CONDITIONS:
        table = pd.read_csv(output_dir / "per_residue_comparison_{}.csv".format(condition),
                            dtype={"author_residue_id": str})
        fig, axes = plt.subplots(4, 1, figsize=(13, 14), sharex=False)
        for axis, chain in zip(axes, ("A", "C", "D", "E")):
            group = table[table.chain == chain].reset_index(drop=True)
            x = np.arange(len(group))
            axis.fill_between(x, group.md_rmsf - group.md_rmsf_sd,
                              group.md_rmsf + group.md_rmsf_sd, color="#4c72b0", alpha=0.25,
                              label="MD mean +/- replica SD")
            axis.plot(x, group.md_rmsf, color="#4c72b0", lw=1.6, label="MD chain-fitted RMSF")
            axis.plot(x, group.pred_rmsf, color="#c44e52", lw=1.4, label="FlexPPI 5-fold mean")
            axis.fill_between(x, group.pred_rmsf - group.pred_rmsf_fold_sd,
                              group.pred_rmsf + group.pred_rmsf_fold_sd, color="#c44e52",
                              alpha=0.15, label="model fold SD")
            for region, color in region_color.items():
                mask = (group.region == region).to_numpy()
                if mask.any():
                    axis.axvspan(x[mask].min() - 0.5, x[mask].max() + 0.5, color=color, alpha=0.35)
                    axis.text(x[mask].mean(), axis.get_ylim()[1] * 0.98, region,
                              ha="center", va="top", fontsize=8)
            axis.set_title(chain_titles[chain], fontsize=10)
            axis.set_ylabel("RMSF (A)")
            axis.set_xticks(x[::5])
            axis.set_xticklabels(group.author_residue_id[::5], fontsize=7)
            axis.legend(fontsize=7, loc="upper right")
        axes[-1].set_xlabel("author residue number")
        fig.suptitle("Chen 2026 TCR156 {}: per-residue RMSF, MD versus FlexPPI".format(condition))
        fig.tight_layout(rect=(0, 0, 1, 0.985))
        fig.savefig(figure_dir / "profiles_{}.png".format(condition), dpi=200)
        fig.savefig(figure_dir / "profiles_{}.pdf".format(condition))
        plt.close(fig)

        fig, axes = plt.subplots(4, 1, figsize=(13, 14), sharex=False)
        for axis, chain in zip(axes, ("A", "C", "D", "E")):
            group = table[table.chain == chain].reset_index(drop=True)
            x = np.arange(len(group))
            def zscore(values):
                values = np.asarray(values, dtype=float)
                return (values - values.mean()) / (values.std(ddof=1) or 1.0)
            axis.plot(x, zscore(group.md_rmsf), color="#4c72b0", lw=1.6, label="MD (z)")
            axis.plot(x, zscore(group.pred_rmsf), color="#c44e52", lw=1.4, label="FlexPPI (z)")
            for region, color in region_color.items():
                mask = (group.region == region).to_numpy()
                if mask.any():
                    axis.axvspan(x[mask].min() - 0.5, x[mask].max() + 0.5, color=color, alpha=0.35)
                    axis.text(x[mask].mean(), axis.get_ylim()[1] * 0.98, region,
                              ha="center", va="top", fontsize=8)
            value = _spearman(group.pred_rmsf.to_numpy(float), group.md_rmsf.to_numpy(float))
            axis.set_title("{}  (chain Spearman {})".format(
                chain_titles[chain], "n/a" if value is None else "{:.2f}".format(value)), fontsize=10)
            axis.set_ylabel("RMSF z-score")
            axis.set_xticks(x[::5])
            axis.set_xticklabels(group.author_residue_id[::5], fontsize=7)
            axis.legend(fontsize=7, loc="upper right")
        axes[-1].set_xlabel("author residue number")
        fig.suptitle("Chen 2026 TCR156 {}: within-chain standardized RMSF".format(condition))
        fig.tight_layout(rect=(0, 0, 1, 0.985))
        fig.savefig(figure_dir / "profiles_z_{}.png".format(condition), dpi=200)
        fig.savefig(figure_dir / "profiles_z_{}.pdf".format(condition))
        plt.close(fig)

    delta = pd.read_csv(output_dir / "wt_vs_s32m_per_residue.csv",
                        dtype={"author_residue_id": str})
    fig, axes = plt.subplots(2, 2, figsize=(13, 8))
    for column, chain in enumerate(("D", "E")):
        group = delta[delta.chain == chain].reset_index(drop=True)
        x = np.arange(len(group))
        axis = axes[0, column]
        axis.axhline(0.0, color="grey", lw=0.8)
        axis.fill_between(x, group.md_delta - group.md_delta_sd,
                          group.md_delta + group.md_delta_sd, color="#4c72b0", alpha=0.25)
        axis.plot(x, group.md_delta, color="#4c72b0", lw=1.6, label="MD delta (S32M - WT)")
        axis.plot(x, group.pred_delta, color="#c44e52", lw=1.4, label="FlexPPI delta")
        for region, color in region_color.items():
            mask = (group.region == region).to_numpy()
            if mask.any():
                axis.axvspan(x[mask].min() - 0.5, x[mask].max() + 0.5, color=color, alpha=0.35)
        axis.set_title("chain {} delta RMSF".format(chain), fontsize=10)
        axis.set_ylabel("delta RMSF (A)")
        axis.set_xticks(x[::5])
        axis.set_xticklabels(group.author_residue_id[::5], fontsize=7)
        axis.legend(fontsize=7)
        axis = axes[1, column]
        axis.axhline(0.0, color="grey", lw=0.8)
        axis.axvline(0.0, color="grey", lw=0.8)
        axis.scatter(group.md_delta, group.pred_delta, s=12, alpha=0.7)
        axis.set_xlabel("MD delta RMSF (A)")
        axis.set_ylabel("FlexPPI delta RMSF (A)")
        axis.set_title("chain {}: model delta vs MD delta".format(chain), fontsize=10)
    fig.suptitle("WT versus S32M (TCR156 alpha S32M catch-bond mutant)")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    fig.savefig(figure_dir / "wt_vs_s32m_delta.png", dpi=200)
    fig.savefig(figure_dir / "wt_vs_s32m_delta.pdf")
    plt.close(fig)
    print("figures written to", figure_dir, flush=True)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=("extract", "predict", "compare", "residues",
                                          "figures", "run"))
    parser.add_argument("--data-root", default=r"E:\data\CatchBond")
    parser.add_argument("--labels-dir", default=str(ROOT / "artifacts/chen2026_tcr156/labels"))
    parser.add_argument("--predictions-dir", default=str(ROOT / "artifacts/chen2026_tcr156/predictions"))
    parser.add_argument("--output-dir", default=str(ROOT / "artifacts/chen2026_tcr156"))
    parser.add_argument("--device", default="cpu")
    parser.add_argument("--provenance-only", action="store_true")
    args = parser.parse_args()
    data_root = Path(args.data_root)
    labels_dir = Path(args.labels_dir)
    predictions_dir = Path(args.predictions_dir)
    output_dir = Path(args.output_dir)

    implementation = sha256(Path(__file__))
    if args.stage in ("extract", "run") and not args.provenance_only:
        extract(data_root, labels_dir)
    if args.stage in ("predict", "run") and not args.provenance_only:
        predict(data_root, predictions_dir, args.device)
    if args.stage in ("compare", "run") and not args.provenance_only:
        compare(data_root, labels_dir, predictions_dir, output_dir)
    if args.stage in ("residues", "run") and not args.provenance_only:
        residues(data_root, labels_dir, predictions_dir, output_dir)
    if args.stage in ("figures", "run") and not args.provenance_only:
        figures(labels_dir, predictions_dir, output_dir)
    if args.provenance_only:
        print(json.dumps({"implementation_sha256": implementation,
                          "started": time.strftime("%Y-%m-%dT%H:%M:%S")}, indent=2))


if __name__ == "__main__":
    main()
