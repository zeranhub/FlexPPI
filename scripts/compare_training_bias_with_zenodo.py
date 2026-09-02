"""Compare FlexPPI paired-training response bias before and after Zenodo 6638504.

The script is deliberately model-free.  It puts the current DynaRepo paired
fine-tuning set and the 16 strict Zenodo chain sides on the same definitions:
100-ns RMSF windows, chain-wise rigid-body fitting, target/partner C-alpha
distance, and PDB-level statistics.  It also reports the actual pair-level
sampling weight implied by the current batch-size-one trainer.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from flexppi.constants import THREE_TO_ONE
from flexppi.data.labels import _chainwise_rmsf, load_label_file
from flexppi.data.pairing import global_sequence_alignment
from flexppi.data.structures import parse_pdb


DISTANCE_LABELS = ["<=7A", "7-9A", "9-12A", ">12A"]


def read_jsonl(path: Path) -> List[dict[str, Any]]:
    rows = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                rows.append(json.loads(line))
    return rows


def safe_spearman(first: Sequence[float], second: Sequence[float]) -> float:
    x = np.asarray(first, dtype=np.float64)
    y = np.asarray(second, dtype=np.float64)
    valid = np.isfinite(x) & np.isfinite(y)
    if valid.sum() < 6 or np.std(x[valid]) < 1e-12 or np.std(y[valid]) < 1e-12:
        return float("nan")
    return float(spearmanr(x[valid], y[valid]).statistic)


def bootstrap_mean(
    values: Sequence[float],
    repetitions: int = 5000,
    seed: int = 29,
) -> Tuple[float, float, float]:
    array = np.asarray(values, dtype=np.float64)
    array = array[np.isfinite(array)]
    if not len(array):
        return float("nan"), float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    samples = rng.choice(array, size=(repetitions, len(array)), replace=True).mean(axis=1)
    return float(array.mean()), float(np.quantile(samples, 0.025)), float(np.quantile(samples, 0.975))


def distance_bins(distance: pd.Series) -> pd.Series:
    return pd.cut(
        distance,
        bins=[-np.inf, 7.0, 9.0, 12.0, np.inf],
        labels=DISTANCE_LABELS,
        right=True,
    ).astype(str)


def load_dynarepo(
    residue_path: Path,
    block_path: Path,
) -> pd.DataFrame:
    residue = pd.read_csv(residue_path)
    block = pd.read_csv(block_path)
    keys = ["pdb_id", "pair_id", "target_local_index"]
    pivot = block.pivot_table(index=keys, columns="block", values="delta_rmsf_block", aggfunc="first")
    values = pivot.to_numpy(dtype=np.float64)
    mean = np.nanmean(values, axis=1)
    sd = np.nanstd(values, axis=1, ddof=1)
    count = np.sum(np.isfinite(values), axis=1)
    # t_0.975,df=4 because every current DynaRepo pair has five 100-ns blocks.
    half_width = 2.7764451051977987 * sd / np.sqrt(np.maximum(count, 1))
    stability = pivot.reset_index()[keys].copy()
    stability["ci95_low"] = mean - half_width
    stability["ci95_high"] = mean + half_width
    stability["block_count"] = count
    stability["block_sd_delta_rmsf"] = sd

    selected = residue[
        keys
        + [
            "target_global_index",
            "chain_id",
            "author_residue_id",
            "residue_name",
            "static_distance",
            "delta_rmsf",
            "delta_rmsf_variance",
        ]
    ].merge(stability, on=keys, how="left", validate="one_to_one")
    selected["dataset"] = "DynaRepo-current"
    selected["side"] = selected["chain_id"]
    selected["window_count"] = selected["block_count"]
    selected["label_qc_abs_error"] = np.nan
    return annotate_response(selected)


def standard_residues(universe: Any) -> Tuple[List[Any], List[str]]:
    residues = []
    chains = []
    for residue in universe.residues:
        name = str(residue.resname).upper()
        atoms = {str(atom.name).upper(): atom for atom in residue.atoms}
        if name not in THREE_TO_ONE or "CA" not in atoms:
            continue
        residues.append(residue)
        chain = str(getattr(residue, "chainID", "")).strip()
        if not chain:
            chain = str(getattr(residue, "segid", "")).strip()
        chains.append(chain or "_")
    return residues, chains


def state_block_rmsf(record: Mapping[str, Any], window_frames: int) -> np.ndarray:
    import MDAnalysis as mda

    structure_path = str(record["structure_path"])
    trajectories = [str(path) for path in record["trajectory_paths"]]
    parsed = parse_pdb(structure_path)
    blocks: List[np.ndarray] = []
    for trajectory_path in trajectories:
        universe = mda.Universe(structure_path, trajectory_path)
        residues, _ = standard_residues(universe)
        if len(residues) != parsed.num_residues:
            raise ValueError(
                f"PDB/XTC residue mismatch for {record['sample_id']}: {parsed.num_residues} vs {len(residues)}"
            )
        ca_indices = []
        for residue in residues:
            ca_indices.append(next(int(atom.index) for atom in residue.atoms if str(atom.name).upper() == "CA"))
        frames = np.stack(
            [np.asarray(ts.positions[ca_indices], dtype=np.float32).copy() for ts in universe.trajectory],
            axis=0,
        )
        for start in range(0, len(frames) - window_frames + 1, window_frames):
            blocks.append(_chainwise_rmsf(frames[start : start + window_frames], parsed.chain_index))
    if not blocks:
        raise ValueError(f"No complete RMSF windows for {record['sample_id']}")
    return np.stack(blocks)


def zenodo_residues(
    manifest_path: Path,
    window_frames: int,
) -> Tuple[pd.DataFrame, Dict[str, str], pd.DataFrame]:
    pairs = read_jsonl(manifest_path)
    state_cache: Dict[str, np.ndarray] = {}
    sequence_by_pair: Dict[str, str] = {}
    rows: List[pd.DataFrame] = []
    qc_rows: List[dict[str, object]] = []
    for pair_index, pair in enumerate(pairs, start=1):
        pair_id = str(pair["pair_id"])
        pdb_id = str(pair["metadata"]["pdb_id"])
        print(f"[{pair_index}/{len(pairs)}] {pair_id}", flush=True)
        apo = pair["apo"]
        holo = pair["holo"]
        for state in (apo, holo):
            sample_id = str(state["sample_id"])
            if sample_id not in state_cache:
                state_cache[sample_id] = state_block_rmsf(state, window_frames)

        apo_indices = np.asarray(pair["apo_indices"], dtype=np.int64)
        holo_indices = np.asarray(pair["holo_indices"], dtype=np.int64)
        apo_blocks = state_cache[str(apo["sample_id"])][:, apo_indices]
        holo_blocks = state_cache[str(holo["sample_id"])][:, holo_indices]
        block_count = min(len(apo_blocks), len(holo_blocks))
        delta_blocks = holo_blocks[:block_count] - apo_blocks[:block_count]
        block_mean = delta_blocks.mean(axis=0)
        block_sd = delta_blocks.std(axis=0, ddof=1)
        if block_count == 3:
            t_critical = 4.302652729749464
        elif block_count == 5:
            t_critical = 2.7764451051977987
        else:
            raise ValueError(f"Unexpected Zenodo block count {block_count} for {pair_id}")
        half_width = t_critical * block_sd / np.sqrt(block_count)

        parsed = parse_pdb(str(holo["structure_path"]))
        ca = parsed.backbone[:, 1]
        ca_mask = parsed.backbone_mask[:, 1]
        target_set = set(int(value) for value in holo_indices)
        partner_indices = np.asarray(
            [index for index in range(parsed.num_residues) if index not in target_set and ca_mask[index]],
            dtype=np.int64,
        )
        target_ca = ca[holo_indices]
        partner_ca = ca[partner_indices]
        if not len(partner_ca):
            raise ValueError(f"No partner residues found for {pair_id}")
        static_distance = np.sqrt(
            np.sum((target_ca[:, None, :] - partner_ca[None, :, :]) ** 2, axis=-1)
        ).min(axis=1)

        apo_values, apo_variances, _ = load_label_file(str(apo["label_path"]))
        holo_values, holo_variances, _ = load_label_file(str(holo["label_path"]))
        label_delta = holo_values[holo_indices, 0] - apo_values[apo_indices, 0]
        label_variance = holo_variances[holo_indices, 0] + apo_variances[apo_indices, 0]
        sequence = "".join(parsed.sequence[int(index)] for index in holo_indices)
        sequence_by_pair[pair_id] = sequence
        table = pd.DataFrame(
            {
                "dataset": "Zenodo-6638504",
                "pdb_id": pdb_id,
                "pair_id": pair_id,
                "side": str(pair["metadata"]["side"]),
                "target_local_index": np.arange(len(holo_indices), dtype=np.int64),
                "target_global_index": holo_indices,
                "chain_id": [parsed.chain_ids[int(index)] for index in holo_indices],
                "author_residue_id": [parsed.author_residue_ids[int(index)] for index in holo_indices],
                "residue_name": [parsed.residue_names[int(index)] for index in holo_indices],
                "static_distance": static_distance,
                "delta_rmsf": label_delta,
                "delta_rmsf_variance": label_variance,
                "ci95_low": block_mean - half_width,
                "ci95_high": block_mean + half_width,
                "block_count": block_count,
                "block_sd_delta_rmsf": block_sd,
                "window_count": block_count,
                "label_qc_abs_error": np.abs(label_delta - block_mean),
            }
        )
        rows.append(annotate_response(table))
        qc_rows.append(
            {
                "pdb_id": pdb_id,
                "pair_id": pair_id,
                "residues": len(table),
                "windows": block_count,
                "label_vs_recomputed_spearman": safe_spearman(label_delta, block_mean),
                "label_vs_recomputed_mae": float(np.mean(np.abs(label_delta - block_mean))),
            }
        )
    return pd.concat(rows, ignore_index=True), sequence_by_pair, pd.DataFrame(qc_rows)


def annotate_response(table: pd.DataFrame) -> pd.DataFrame:
    table = table.copy()
    table["distance_bin"] = distance_bins(table["static_distance"])
    table["interface_8a"] = table["static_distance"] <= 8.0
    table["stiffening"] = table["delta_rmsf"] < 0
    table["softening"] = table["delta_rmsf"] > 0
    table["robust_stiffening"] = table["ci95_high"] < 0
    table["robust_softening"] = table["ci95_low"] > 0
    return table


def per_unit_metrics(table: pd.DataFrame, unit_columns: Sequence[str]) -> pd.DataFrame:
    rows = []
    for key, group in table.groupby(list(unit_columns), sort=False):
        if not isinstance(key, tuple):
            key = (key,)
        base = dict(zip(unit_columns, key))
        for region, selected in (
            ("all", group),
            ("interface_8a", group[group["interface_8a"]]),
            ("non_interface", group[~group["interface_8a"]]),
        ):
            if selected.empty:
                continue
            row = dict(base)
            row.update(
                {
                    "region": region,
                    "residues": len(selected),
                    "mean_delta_rmsf": float(selected["delta_rmsf"].mean()),
                    "mean_abs_delta_rmsf": float(np.abs(selected["delta_rmsf"]).mean()),
                    "fraction_stiffening": float(selected["stiffening"].mean()),
                    "fraction_softening": float(selected["softening"].mean()),
                    "fraction_robust_stiffening": float(selected["robust_stiffening"].mean()),
                    "fraction_robust_softening": float(selected["robust_softening"].mean()),
                    "distance_response_spearman": safe_spearman(selected["static_distance"], selected["delta_rmsf"]),
                    "median_label_variance": float(selected["delta_rmsf_variance"].median()),
                    "has_both_robust_directions": bool(
                        selected["robust_stiffening"].any() and selected["robust_softening"].any()
                    ),
                    "mean_interface_positive": bool(
                        region == "interface_8a" and selected["delta_rmsf"].mean() > 0
                    ),
                }
            )
            rows.append(row)
    return pd.DataFrame(rows)


def summarise_datasets(combined: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    dataset_tables = {
        "DynaRepo-current": combined[combined["dataset"] == "DynaRepo-current"],
        "Zenodo-6638504": combined[combined["dataset"] == "Zenodo-6638504"],
        "Combined": combined,
    }
    summary_rows = []
    distance_rows = []
    per_pdb_all = []
    metrics = [
        "mean_delta_rmsf",
        "mean_abs_delta_rmsf",
        "fraction_stiffening",
        "fraction_softening",
        "fraction_robust_stiffening",
        "fraction_robust_softening",
        "distance_response_spearman",
        "median_label_variance",
    ]
    for scenario, table in dataset_tables.items():
        per_pdb = per_unit_metrics(table, ["dataset", "pdb_id"])
        per_pdb["scenario"] = scenario
        per_pdb_all.append(per_pdb)
        for region, region_table in per_pdb.groupby("region", sort=False):
            row: Dict[str, object] = {
                "scenario": scenario,
                "region": region,
                "pdbs": int(region_table[["dataset", "pdb_id"]].drop_duplicates().shape[0]),
                "pairs": int(table["pair_id"].nunique()),
                "residues": int(table["interface_8a"].sum()) if region == "interface_8a" else len(table),
            }
            for metric in metrics:
                mean, low, high = bootstrap_mean(region_table[metric])
                row[metric] = mean
                row[metric + "_ci_low"] = low
                row[metric + "_ci_high"] = high
            row["pdbs_with_mean_interface_softening"] = int(region_table["mean_interface_positive"].sum())
            row["pdbs_with_both_robust_directions"] = int(region_table["has_both_robust_directions"].sum())
            summary_rows.append(row)

        for label in DISTANCE_LABELS:
            subset = table[table["distance_bin"] == label]
            pdb_rows = []
            for (dataset, pdb_id), group in subset.groupby(["dataset", "pdb_id"], sort=False):
                pdb_rows.append(
                    {
                        "dataset": dataset,
                        "pdb_id": pdb_id,
                        "mean_delta_rmsf": float(group["delta_rmsf"].mean()),
                        "mean_abs_delta_rmsf": float(np.abs(group["delta_rmsf"]).mean()),
                        "fraction_stiffening": float(group["stiffening"].mean()),
                        "fraction_robust_stiffening": float(group["robust_stiffening"].mean()),
                        "fraction_robust_softening": float(group["robust_softening"].mean()),
                    }
                )
            pdb_table = pd.DataFrame(pdb_rows)
            row = {
                "scenario": scenario,
                "distance_bin": label,
                "pdbs": len(pdb_table),
                "residues": len(subset),
            }
            for metric in [
                "mean_delta_rmsf",
                "mean_abs_delta_rmsf",
                "fraction_stiffening",
                "fraction_robust_stiffening",
                "fraction_robust_softening",
            ]:
                mean, low, high = bootstrap_mean(pdb_table[metric])
                row[metric] = mean
                row[metric + "_ci_low"] = low
                row[metric + "_ci_high"] = high
            distance_rows.append(row)
    return pd.DataFrame(summary_rows), pd.DataFrame(distance_rows), pd.concat(per_pdb_all, ignore_index=True)


def training_weight_scenarios(combined: pd.DataFrame) -> pd.DataFrame:
    rows = []
    total_pairs = combined[["dataset", "pair_id"]].drop_duplicates().shape[0]
    total_pdbs = combined[["dataset", "pdb_id"]].drop_duplicates().shape[0]
    total_residues = len(combined)
    for dataset, group in combined.groupby("dataset", sort=False):
        pair_count = group["pair_id"].nunique()
        pdb_count = group["pdb_id"].nunique()
        rows.append(
            {
                "dataset": dataset,
                "pairs": pair_count,
                "pdbs": pdb_count,
                "residues": len(group),
                "interface_8a_residues": int(group["interface_8a"].sum()),
                "share_of_optimizer_steps_batch1_pair_sampling": pair_count / total_pairs,
                "share_if_pdb_equal_weighted": pdb_count / total_pdbs,
                "share_if_residue_pooled": len(group) / total_residues,
                "mean_pairs_per_pdb": pair_count / pdb_count,
            }
        )
    return pd.DataFrame(rows)


def sampling_scheme_summary(combined: pd.DataFrame) -> pd.DataFrame:
    """Show how the same labels look under residue, pair, or PDB weighting."""

    scenarios = {
        "DynaRepo-current": combined[combined["dataset"] == "DynaRepo-current"],
        "Zenodo-6638504": combined[combined["dataset"] == "Zenodo-6638504"],
        "Combined": combined,
    }
    rows = []
    for scenario, table in scenarios.items():
        for region, selected in (
            ("all", table),
            ("interface_8a", table[table["interface_8a"]]),
        ):
            for weighting, units in (
                ("residue_pooled", None),
                ("pair_equal", ["dataset", "pair_id"]),
                ("pdb_equal", ["dataset", "pdb_id"]),
            ):
                if units is None:
                    metric_table = pd.DataFrame(
                        [
                            {
                                "mean_delta_rmsf": float(selected["delta_rmsf"].mean()),
                                "mean_abs_delta_rmsf": float(np.abs(selected["delta_rmsf"]).mean()),
                                "fraction_stiffening": float(selected["stiffening"].mean()),
                                "fraction_robust_stiffening": float(selected["robust_stiffening"].mean()),
                                "fraction_robust_softening": float(selected["robust_softening"].mean()),
                                "median_label_variance": float(selected["delta_rmsf_variance"].median()),
                            }
                        ]
                    )
                else:
                    unit_rows = []
                    for _, group in selected.groupby(units, sort=False):
                        unit_rows.append(
                            {
                                "mean_delta_rmsf": float(group["delta_rmsf"].mean()),
                                "mean_abs_delta_rmsf": float(np.abs(group["delta_rmsf"]).mean()),
                                "fraction_stiffening": float(group["stiffening"].mean()),
                                "fraction_robust_stiffening": float(group["robust_stiffening"].mean()),
                                "fraction_robust_softening": float(group["robust_softening"].mean()),
                                "median_label_variance": float(group["delta_rmsf_variance"].median()),
                            }
                        )
                    metric_table = pd.DataFrame(unit_rows)
                row = {
                    "scenario": scenario,
                    "region": region,
                    "weighting": weighting,
                    "units": len(metric_table),
                }
                for metric in metric_table.columns:
                    row[metric] = float(metric_table[metric].mean())
                rows.append(row)
    return pd.DataFrame(rows)


def distance_prior_fit(combined: pd.DataFrame) -> pd.DataFrame:
    """Fit a transparent pair-equal one-feature prior for scale diagnostics."""

    scenarios = {
        "DynaRepo-current": combined[combined["dataset"] == "DynaRepo-current"],
        "Zenodo-6638504": combined[combined["dataset"] == "Zenodo-6638504"],
        "Combined": combined,
    }
    probe_distances = np.asarray([0.0, 5.0, 8.0, 12.0, 20.0], dtype=np.float64)
    rows = []
    for scenario, table in scenarios.items():
        x = -np.exp(-table["static_distance"].to_numpy(dtype=np.float64) / 5.0)
        y = table["delta_rmsf"].to_numpy(dtype=np.float64)
        pair_sizes = table.groupby(["dataset", "pair_id"])["pair_id"].transform("size").to_numpy(dtype=np.float64)
        weights = 1.0 / pair_sizes
        design = np.column_stack([np.ones(len(x)), x])
        coefficients = np.linalg.solve(
            design.T @ (weights[:, None] * design),
            design.T @ (weights * y),
        )
        predictions = coefficients[0] + coefficients[1] * (-np.exp(-probe_distances / 5.0))
        row: Dict[str, object] = {
            "scenario": scenario,
            "intercept": float(coefficients[0]),
            "distance_exp_coefficient": float(coefficients[1]),
        }
        for distance, prediction in zip(probe_distances, predictions):
            row[f"predicted_delta_rmsf_at_{int(distance)}a"] = float(prediction)
        rows.append(row)
    return pd.DataFrame(rows)


def dynarepo_sequences(dynarepo: pd.DataFrame) -> Dict[str, str]:
    sequences = {}
    for pair_id, group in dynarepo.groupby("pair_id", sort=False):
        ordered = group.sort_values("target_local_index")
        sequences[str(pair_id)] = "".join(
            THREE_TO_ONE.get(str(name).upper(), "X") for name in ordered["residue_name"]
        )
    return sequences


def sequence_novelty(
    zenodo_sequences: Mapping[str, str],
    dynarepo_sequence_map: Mapping[str, str],
    dynarepo: pd.DataFrame,
) -> pd.DataFrame:
    pair_to_pdb = dynarepo.groupby("pair_id", sort=False)["pdb_id"].first().to_dict()
    rows = []
    for zenodo_pair, query in zenodo_sequences.items():
        candidates = []
        for dynarepo_pair, target in dynarepo_sequence_map.items():
            alignment = global_sequence_alignment(query, target)
            candidates.append(
                (
                    alignment.identity,
                    alignment.coverage,
                    alignment.score,
                    dynarepo_pair,
                    str(pair_to_pdb[dynarepo_pair]),
                )
            )
        best = max(candidates, key=lambda value: (value[0] * value[1], value[1], value[0], value[2]))
        rows.append(
            {
                "zenodo_pair_id": zenodo_pair,
                "zenodo_pdb_id": zenodo_pair.split("__", 1)[0],
                "length": len(query),
                "best_dynarepo_pair_id": best[3],
                "best_dynarepo_pdb_id": best[4],
                "global_alignment_identity": best[0],
                "query_exact_match_coverage": best[1],
                "identity_times_coverage": best[0] * best[1],
            }
        )
    return pd.DataFrame(rows)


def plot_comparison(
    summary: pd.DataFrame,
    distance_summary: pd.DataFrame,
    weights: pd.DataFrame,
    per_pdb: pd.DataFrame,
    output: Path,
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(13.2, 9.4))
    scenarios = ["DynaRepo-current", "Zenodo-6638504", "Combined"]
    colors = ["#4C78A8", "#F58518", "#54A24B"]
    interface = summary[summary["region"] == "interface_8a"].set_index("scenario")
    means = interface.loc[scenarios, "mean_delta_rmsf"].to_numpy()
    lower = means - interface.loc[scenarios, "mean_delta_rmsf_ci_low"].to_numpy()
    upper = interface.loc[scenarios, "mean_delta_rmsf_ci_high"].to_numpy() - means
    axes[0, 0].bar(np.arange(3), means, color=colors)
    axes[0, 0].errorbar(np.arange(3), means, yerr=[lower, upper], fmt="none", color="black", capsize=4)
    axes[0, 0].axhline(0, color="0.25", linewidth=0.8)
    axes[0, 0].set_xticks(np.arange(3), ["Current", "Zenodo", "Combined"])
    axes[0, 0].set_ylabel(r"PDB-macro interface $\Delta$RMSF ($\AA$)")
    axes[0, 0].set_title("A. Interface bias remains negative")

    x = np.arange(3)
    stiff = interface.loc[scenarios, "fraction_robust_stiffening"].to_numpy()
    soft = interface.loc[scenarios, "fraction_robust_softening"].to_numpy()
    axes[0, 1].bar(x - 0.18, stiff, width=0.36, color="#E45756", label="robust stiffening")
    axes[0, 1].bar(x + 0.18, soft, width=0.36, color="#72B7B2", label="robust softening")
    axes[0, 1].set_xticks(x, ["Current", "Zenodo", "Combined"])
    axes[0, 1].set_ylabel("PDB-macro residue fraction")
    axes[0, 1].set_title("B. Zenodo does not add many softening interfaces")
    axes[0, 1].legend(frameon=False, fontsize=9)

    positions = np.arange(len(DISTANCE_LABELS))
    for scenario, color in zip(scenarios, colors):
        group = distance_summary[distance_summary["scenario"] == scenario].set_index("distance_bin").loc[DISTANCE_LABELS]
        axes[1, 0].plot(positions, group["mean_delta_rmsf"], marker="o", color=color, label=scenario)
    axes[1, 0].axhline(0, color="0.25", linewidth=0.8)
    axes[1, 0].set_xticks(positions, DISTANCE_LABELS)
    axes[1, 0].set_ylabel(r"PDB-macro $\Delta$RMSF ($\AA$)")
    axes[1, 0].set_xlabel("Initial target-partner C-alpha distance")
    axes[1, 0].set_title("C. The same distance gradient is reinforced")
    axes[1, 0].legend(frameon=False, fontsize=8)

    weight_table = weights.set_index("dataset").loc[["DynaRepo-current", "Zenodo-6638504"]]
    labels = ["pair-step", "PDB-equal", "residue-pooled"]
    dyna_values = [
        weight_table.loc["DynaRepo-current", "share_of_optimizer_steps_batch1_pair_sampling"],
        weight_table.loc["DynaRepo-current", "share_if_pdb_equal_weighted"],
        weight_table.loc["DynaRepo-current", "share_if_residue_pooled"],
    ]
    zenodo_values = [1.0 - value for value in dyna_values]
    axes[1, 1].bar(labels, dyna_values, color=colors[0], label="current")
    axes[1, 1].bar(labels, zenodo_values, bottom=dyna_values, color=colors[1], label="Zenodo")
    axes[1, 1].set_ylim(0, 1)
    axes[1, 1].set_ylabel("Share after concatenation")
    axes[1, 1].set_title("D. Two chain sides over-weight each Zenodo PDB")
    axes[1, 1].legend(frameon=False, fontsize=9)

    fig.suptitle("Paired-training rigidification bias: current DynaRepo vs Zenodo 6638504", fontsize=15)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=220, bbox_inches="tight")
    fig.savefig(output.with_suffix(".pdf"), bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--dynarepo-dir", type=Path, default=Path("artifacts/dynarepo_interface_dynamics")
    )
    parser.add_argument(
        "--zenodo-pair-manifest",
        type=Path,
        default=Path("artifacts/zenodo_6638504/data/zenodo_6638504_pairs.jsonl"),
    )
    parser.add_argument("--zenodo-window-frames", type=int, default=2500)
    parser.add_argument(
        "--output-dir", type=Path, default=Path("artifacts/training_bias_with_zenodo_6638504")
    )
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    dynarepo = load_dynarepo(
        args.dynarepo_dir / "residue_features.csv",
        args.dynarepo_dir / "block_residue_features.csv",
    )
    zenodo, zenodo_sequences, qc = zenodo_residues(
        args.zenodo_pair_manifest, args.zenodo_window_frames
    )
    combined = pd.concat([dynarepo, zenodo], ignore_index=True, sort=False)
    summary, distance_summary, per_pdb = summarise_datasets(combined)
    weights = training_weight_scenarios(combined)
    sampling_summary = sampling_scheme_summary(combined)
    prior_fit = distance_prior_fit(combined)
    novelty = sequence_novelty(zenodo_sequences, dynarepo_sequences(dynarepo), dynarepo)

    zenodo.to_csv(args.output_dir / "zenodo_residue_response_calpha8.csv", index=False)
    combined.to_csv(args.output_dir / "combined_residue_response.csv", index=False)
    qc.to_csv(args.output_dir / "zenodo_label_qc.csv", index=False)
    summary.to_csv(args.output_dir / "dataset_bias_summary.csv", index=False)
    distance_summary.to_csv(args.output_dir / "distance_bin_bias.csv", index=False)
    per_pdb.to_csv(args.output_dir / "per_pdb_bias.csv", index=False)
    weights.to_csv(args.output_dir / "training_weight_scenarios.csv", index=False)
    sampling_summary.to_csv(args.output_dir / "sampling_scheme_bias.csv", index=False)
    prior_fit.to_csv(args.output_dir / "distance_prior_fit.csv", index=False)
    novelty.to_csv(args.output_dir / "zenodo_sequence_novelty.csv", index=False)
    plot_comparison(
        summary,
        distance_summary,
        weights,
        per_pdb,
        args.output_dir / "combined_training_bias.png",
    )


if __name__ == "__main__":
    main()
