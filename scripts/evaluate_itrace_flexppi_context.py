"""Strict zero-shot FlexPPI context ablation on the ITrace dataset.

This script reuses the full-complex FlexPPI predictions produced by
``evaluate_itrace_benchmark.py`` and runs the *same frozen checkpoint* on each
chain extracted from the same bound-state topology.  The comparison answers
whether multichain context improves bound-state RMSF prediction.  It does not
estimate an apo-to-holo flexibility change because ITrace has no matched apo
structures for these complexes.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import time
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon
import torch

from flexppi.comparison import _bootstrap_mean_ci, _system_metrics
from flexppi.data.dataset import collate_states
from flexppi.data.structures import ProteinStructure, parse_pdb
from flexppi.inference import _structure_example, load_model
from flexppi.trainer import move_to_device


METRICS = (
    "spearman",
    "pearson",
    "z_mae",
    "z_rmse",
    "top20_auroc",
    "top20_average_precision",
    "top20_precision_at_k",
    "absolute_mae_angstrom",
    "absolute_rmse_angstrom",
)
ADVANTAGE_DIRECTION = {
    "spearman": 1.0,
    "pearson": 1.0,
    "z_mae": -1.0,
    "z_rmse": -1.0,
    "top20_auroc": 1.0,
    "top20_average_precision": 1.0,
    "top20_precision_at_k": 1.0,
    "absolute_mae_angstrom": -1.0,
    "absolute_rmse_angstrom": -1.0,
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--benchmark-dir", default="output/itrace_benchmark")
    parser.add_argument("--checkpoint", default="artifacts/pretrain/best.pt")
    parser.add_argument("--output-dir", default="output/itrace_zero_shot_context")
    parser.add_argument("--bootstrap-repetitions", type=int, default=2000)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--reuse", action="store_true")
    parser.add_argument("--seed", type=int, default=7)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def json_dump(path: Path, value: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False),
        encoding="utf-8",
    )


def _json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): _json_safe(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    return value


def _renumber_atom(line: str, serial: int) -> str:
    value = line.rstrip("\r\n").ljust(80)
    return value[:6] + f"{serial:5d}" + value[11:] + "\n"


def _header_before_coordinates(lines: Sequence[str]) -> List[str]:
    header: List[str] = []
    for line in lines:
        record = line[:6].strip()
        if record in {"MODEL", "ATOM", "HETATM"}:
            break
        if record not in {"END", "MASTER", "CONECT"}:
            header.append(line if line.endswith("\n") else line + "\n")
    return header


def extract_chain(source: Path, destination: Path, chain: str) -> ProteinStructure:
    """Write one chain without changing any retained atom coordinates."""

    lines = source.read_text(encoding="utf-8", errors="replace").splitlines(True)
    atoms = [
        line
        for line in lines
        if line.startswith("ATOM  ") and line[21:22].strip() == chain
    ]
    if not atoms:
        raise ValueError(f"Chain {chain!r} has no ATOM records in {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as handle:
        handle.writelines(_header_before_coordinates(lines))
        handle.write(
            "REMARK 900 FLEXPPI ZERO-SHOT CONTEXT ABLATION: "
            f"BOUND-STATE CHAIN {chain} RETAINED; COORDINATES UNCHANGED\n"
        )
        for serial, line in enumerate(atoms, start=1):
            handle.write(_renumber_atom(line, serial))
        handle.write("END\n")
    structure = parse_pdb(str(destination))
    if structure.unique_chains() != [chain]:
        raise ValueError(
            f"Expected extracted chain {[chain]!r}, got {structure.unique_chains()!r}"
        )
    return structure


def author_residue_parts(author_id: str) -> tuple[str, str]:
    match = re.match(r"^(-?\d+)(.*)$", str(author_id).strip())
    if match is None:
        return str(author_id).strip(), ""
    return match.group(1), match.group(2)


@torch.no_grad()
def predict_isolated_chains(
    samples: Sequence[Mapping[str, str]],
    checkpoint: Path,
    prepared_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    model, model_config, device = load_model(str(checkpoint), "auto")
    model.eval()
    rmsf_index = model_config.metrics.index("rmsf")
    prediction_rows: List[Dict[str, Any]] = []
    timing_rows: List[Dict[str, Any]] = []

    for index, sample in enumerate(samples, start=1):
        pdb_id = str(sample["pdb_id"]).lower()
        chain = str(sample["chain"])
        source = Path(str(sample["topology_path"]))
        destination = prepared_dir / f"{pdb_id}_{chain}.pdb"
        started = time.perf_counter()
        try:
            structure = extract_chain(source, destination, chain)
            example = _structure_example(
                str(destination),
                [chain],
                310.0,
                200.0,
                1.0,
            )
            batch = move_to_device(collate_states([example]), device)
            output = model.forward_state(batch)
            prediction = output["mean"][0, :, rmsf_index].cpu().numpy()
            uncertainty = torch.exp(
                0.5 * output["log_variance"][0, :, rmsf_index]
            ).cpu().numpy()
            for residue_index in range(structure.num_residues):
                residue_number, insertion_code = author_residue_parts(
                    structure.author_residue_ids[residue_index]
                )
                prediction_rows.append(
                    {
                        "pdb_id": pdb_id,
                        "chain": chain,
                        "residue_number": residue_number,
                        "insertion_code": insertion_code,
                        "amino_acid_3": structure.residue_names[residue_index],
                        "isolated_flexppi_rmsf_angstrom": float(
                            prediction[residue_index]
                        ),
                        "isolated_flexppi_std_angstrom": float(
                            uncertainty[residue_index]
                        ),
                    }
                )
            status, error = "ok", ""
            residue_count = structure.num_residues
        except Exception as exception:
            status, error = "error", repr(exception)
            residue_count = np.nan
        timing_rows.append(
            {
                "pdb_id": pdb_id,
                "chain": chain,
                "residue_count": residue_count,
                "end_to_end_seconds": time.perf_counter() - started,
                "status": status,
                "error": error,
            }
        )
        if index % 25 == 0 or index == len(samples):
            failures = sum(row["status"] != "ok" for row in timing_rows)
            print(
                f"isolated-chain FlexPPI {index}/{len(samples)} "
                f"(failures={failures})",
                flush=True,
            )

    return pd.DataFrame(prediction_rows), pd.DataFrame(timing_rows)


def normalize_keys(frame: pd.DataFrame) -> pd.DataFrame:
    output = frame.copy()
    output["pdb_id"] = output["pdb_id"].astype(str).str.lower()
    output["chain"] = output["chain"].fillna("").astype(str)
    output["residue_number"] = output["residue_number"].astype(str)
    output["insertion_code"] = output["insertion_code"].fillna("").astype(str)
    return output


def chain_metric_table(frame: pd.DataFrame) -> pd.DataFrame:
    rows: List[Dict[str, Any]] = []
    for (pdb_id, chain), group in frame.groupby(["pdb_id", "chain"], sort=True):
        target = group["target_chainfit_rmsf_angstrom"].to_numpy(dtype=float)
        full = _system_metrics(
            group["full_flexppi_rmsf_angstrom"].to_numpy(dtype=float),
            target,
            absolute_scale=True,
        )
        isolated = _system_metrics(
            group["isolated_flexppi_rmsf_angstrom"].to_numpy(dtype=float),
            target,
            absolute_scale=True,
        )
        row: Dict[str, Any] = {
            "pdb_id": pdb_id,
            "chain": chain,
            "chain_role": str(group["chain_role"].iloc[0]),
            "residue_count": int(len(group)),
        }
        for metric in METRICS:
            row[f"full_{metric}"] = full.get(metric)
            row[f"isolated_{metric}"] = isolated.get(metric)
            if full.get(metric) is not None and isolated.get(metric) is not None:
                row[f"full_advantage_{metric}"] = ADVANTAGE_DIRECTION[metric] * (
                    float(full[metric]) - float(isolated[metric])
                )
        rows.append(row)
    return pd.DataFrame(rows)


def complex_metric_table(chain_metrics: pd.DataFrame) -> pd.DataFrame:
    metric_columns = [
        column
        for column in chain_metrics.columns
        if column.startswith(("full_", "isolated_"))
    ]
    return (
        chain_metrics.groupby("pdb_id", sort=True)[metric_columns]
        .mean(numeric_only=True)
        .reset_index()
        .merge(
            chain_metrics.groupby("pdb_id", sort=True)
            .agg(chain_count=("chain", "nunique"), residue_count=("residue_count", "sum"))
            .reset_index(),
            on="pdb_id",
            how="left",
            validate="one_to_one",
        )
    )


def summarize_method(
    frame: pd.DataFrame,
    prefix: str,
    bootstrap_repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    output: Dict[str, Any] = {}
    for metric in METRICS:
        values = frame[f"{prefix}_{metric}"].dropna().to_numpy(dtype=float)
        if not len(values):
            continue
        output[metric] = {
            "macro_mean": float(np.mean(values)),
            "macro_median": float(np.median(values)),
            "complex_count": int(len(values)),
            "complex_bootstrap_95_ci": _bootstrap_mean_ci(
                values.tolist(), bootstrap_repetitions, seed
            ),
        }
    return output


def summarize_paired(
    frame: pd.DataFrame,
    bootstrap_repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    output: Dict[str, Any] = {}
    for metric in METRICS:
        pair = frame[[f"full_{metric}", f"isolated_{metric}"]].dropna()
        if pair.empty:
            continue
        differences = ADVANTAGE_DIRECTION[metric] * (
            pair[f"full_{metric}"].to_numpy(dtype=float)
            - pair[f"isolated_{metric}"].to_numpy(dtype=float)
        )
        nonzero = differences[np.abs(differences) > 1e-12]
        p_value: Optional[float] = None
        if len(nonzero):
            try:
                p_value = float(wilcoxon(nonzero, alternative="two-sided").pvalue)
            except ValueError:
                pass
        output[metric] = {
            "advantage_definition": (
                "positive means full-complex input is better than "
                "bound-state isolated-chain input"
            ),
            "mean_advantage": float(np.mean(differences)),
            "median_advantage": float(np.median(differences)),
            "complex_bootstrap_95_ci": _bootstrap_mean_ci(
                differences.tolist(), bootstrap_repetitions, seed
            ),
            "win_rate": float(np.mean(differences > 1e-12)),
            "loss_rate": float(np.mean(differences < -1e-12)),
            "tie_rate": float(np.mean(np.abs(differences) <= 1e-12)),
            "paired_complex_count": int(len(differences)),
            "wilcoxon_two_sided_p": p_value,
        }
    return output


def context_response_summary(
    frame: pd.DataFrame,
    bootstrap_repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    masks = {
        "overall": pd.Series(True, index=frame.index),
        "interface_8A": frame["is_tcr_pmhc_interface_8A"].astype(bool),
        "non_interface": ~frame["is_tcr_pmhc_interface_8A"].astype(bool),
    }
    for role in sorted(frame["chain_role"].dropna().astype(str).unique()):
        masks[f"role:{role}"] = frame["chain_role"].astype(str).eq(role)

    output: Dict[str, Any] = {}
    for label, mask in masks.items():
        subset = frame.loc[mask]
        if subset.empty:
            continue
        per_complex = (
            subset.groupby("pdb_id")["context_response_angstrom"]
            .agg(["mean", lambda values: float(np.mean(np.abs(values)))])
            .reset_index()
        )
        per_complex.columns = ["pdb_id", "signed_mean", "absolute_mean"]
        signed = per_complex["signed_mean"].to_numpy(dtype=float)
        absolute = per_complex["absolute_mean"].to_numpy(dtype=float)
        output[label] = {
            "residue_count": int(len(subset)),
            "complex_count": int(per_complex["pdb_id"].nunique()),
            "complex_macro_signed_mean_angstrom": float(np.mean(signed)),
            "complex_macro_signed_mean_95_ci": _bootstrap_mean_ci(
                signed.tolist(), bootstrap_repetitions, seed
            ),
            "complex_macro_absolute_mean_angstrom": float(np.mean(absolute)),
            "complex_macro_absolute_mean_95_ci": _bootstrap_mean_ci(
                absolute.tolist(), bootstrap_repetitions, seed
            ),
        }
    return output


def condition_results(
    complex_metrics: pd.DataFrame,
    bootstrap_repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    return {
        "full_complex": summarize_method(
            complex_metrics, "full", bootstrap_repetitions, seed
        ),
        "isolated_chain": summarize_method(
            complex_metrics, "isolated", bootstrap_repetitions, seed
        ),
        "paired_full_vs_isolated": summarize_paired(
            complex_metrics, bootstrap_repetitions, seed
        ),
    }


def chain_role_results(
    chain_metrics: pd.DataFrame,
    bootstrap_repetitions: int,
    seed: int,
) -> Dict[str, Any]:
    output: Dict[str, Any] = {}
    for role, group in chain_metrics.groupby("chain_role", sort=True):
        output[str(role)] = condition_results(
            group, bootstrap_repetitions, seed
        )
    return output


def markdown_summary(report: Mapping[str, Any]) -> str:
    full = report["external_clean_results"]["full_complex"]
    isolated = report["external_clean_results"]["isolated_chain"]
    paired = report["external_clean_results"]["paired_full_vs_isolated"]
    non_peptide = report["sensitivity_analysis"]["exclude_peptide_chain"]
    non_peptide_full = non_peptide["full_complex"]
    non_peptide_isolated = non_peptide["isolated_chain"]
    non_peptide_paired = non_peptide["paired_full_vs_isolated"]
    role_results = report["sensitivity_analysis"]["by_chain_role"]

    def mean(source: Mapping[str, Any], metric: str) -> float:
        return float(source[metric]["macro_mean"])

    role_advantages = []
    for role in ("mhc", "b2m", "peptide", "tcr_alpha", "tcr_beta"):
        advantage = role_results[role]["paired_full_vs_isolated"]["spearman"][
            "mean_advantage"
        ]
        role_advantages.append(f"{role} {advantage:+.4f}")

    lines = [
        "# ITrace 严格零样本：FlexPPI 复合物上下文消融",
        "",
        "本分析冻结同一个预训练检查点，不进行微调、再训练或阈值选择。完整输入为 ITrace "
        "结合态五链复合物，消融输入为从同一结合态拓扑中原坐标抽出的单链；因此比较的是"
        "“是否提供其他链的结构上下文”，不是有真实 apo 对照的结合前后变化。",
        "",
        "## 外部干净集主结果",
        "",
        f"- PDB 复合物数：{report['cohort']['external_clean_complexes']}（精确训练重叠 "
        f"{', '.join(report['cohort']['excluded_exact_training_overlap'])} 已排除）。",
        f"- Spearman：完整复合物 {mean(full, 'spearman'):.4f}，孤立链 "
        f"{mean(isolated, 'spearman'):.4f}，完整输入优势 "
        f"{paired['spearman']['mean_advantage']:+.4f}。",
        f"- Pearson：完整复合物 {mean(full, 'pearson'):.4f}，孤立链 "
        f"{mean(isolated, 'pearson'):.4f}，完整输入优势 "
        f"{paired['pearson']['mean_advantage']:+.4f}。",
        f"- 绝对 MAE：完整复合物 {mean(full, 'absolute_mae_angstrom'):.4f} Å，"
        f"孤立链 {mean(isolated, 'absolute_mae_angstrom'):.4f} Å；按“越低越好”换算的"
        f"完整输入优势 {paired['absolute_mae_angstrom']['mean_advantage']:+.4f} Å。",
        f"- top-20% AUROC：完整复合物 {mean(full, 'top20_auroc'):.4f}，孤立链 "
        f"{mean(isolated, 'top20_auroc'):.4f}，完整输入优势 "
        f"{paired['top20_auroc']['mean_advantage']:+.4f}。",
        "",
        "## 排除短肽后的敏感性分析",
        "",
        f"- 排除 peptide 链后，Spearman 仍由孤立链 "
        f"{mean(non_peptide_isolated, 'spearman'):.4f} 提高至完整复合物 "
        f"{mean(non_peptide_full, 'spearman'):.4f}，优势 "
        f"{non_peptide_paired['spearman']['mean_advantage']:+.4f}。",
        f"- 但绝对 MAE 由孤立链 "
        f"{mean(non_peptide_isolated, 'absolute_mae_angstrom'):.4f} Å 变为完整复合物 "
        f"{mean(non_peptide_full, 'absolute_mae_angstrom'):.4f} Å；按越低越好换算，完整"
        f"输入优势为 {non_peptide_paired['absolute_mae_angstrom']['mean_advantage']:+.4f} Å，"
        "即绝对标定略差。",
        f"- 各链角色的 Spearman 完整输入优势：{'；'.join(role_advantages)}。这表明"
        "复合物上下文对不同链并非一致有利。",
        "- 孤立短肽产生明显的长度外推失稳，因此包含 peptide 的总体绝对 MAE 改善不能"
        "单独作为跨链机制有效的证据；相关性和检索指标、排除肽链分析更有解释价值。",
        "",
        "## 解释边界",
        "",
        "这里的 `完整预测 − 孤立链预测` 只能称为上下文响应或上下文敏感性。两种输入的"
        "原子坐标完全相同，但完整输入同时改变了跨链接触、输入总长度以及全局注意力归一化；"
        "ITrace 没有匹配 apo 结构，所以它不能被解释成真实的结合诱导 ΔRMSF。",
        "",
        "详细数值见 `itrace_context_report.json`、`metrics/per_complex_metrics.csv` 和 "
        "`raw/common_context_predictions.csv`。",
        "",
    ]
    return "\n".join(lines)


def main() -> None:
    args = parse_arguments()
    benchmark_dir = Path(args.benchmark_dir).resolve()
    checkpoint = Path(args.checkpoint).resolve()
    output_dir = Path(args.output_dir).resolve()
    raw_dir = output_dir / "raw"
    metrics_dir = output_dir / "metrics"
    prepared_dir = output_dir / "prepared_chains"
    for directory in (raw_dir, metrics_dir, prepared_dir):
        directory.mkdir(parents=True, exist_ok=True)

    benchmark_report_path = benchmark_dir / "itrace_benchmark_report.json"
    benchmark_report = json.loads(benchmark_report_path.read_text(encoding="utf-8"))
    expected_sha = benchmark_report["models"]["flexppi"]["checkpoint_sha256"]
    actual_sha = sha256(checkpoint)
    if actual_sha != expected_sha:
        raise ValueError(
            "Checkpoint mismatch: existing full-complex predictions used "
            f"{expected_sha}, but requested checkpoint is {actual_sha}"
        )

    manifest = pd.read_csv(benchmark_dir / "raw" / "itrace_complex_manifest.csv")
    common = pd.read_csv(
        benchmark_dir / "raw" / "common_residue_predictions.csv",
        dtype={"residue_number": str, "insertion_code": str},
    )
    manifest["pdb_id"] = manifest["pdb_id"].astype(str).str.lower()
    common = normalize_keys(common)
    selected_ids = sorted(common["pdb_id"].unique())
    if args.limit is not None:
        selected_ids = selected_ids[: args.limit]
    manifest = manifest[manifest["pdb_id"].isin(selected_ids)].copy()
    common = common[common["pdb_id"].isin(selected_ids)].copy()

    chain_table = (
        common[["pdb_id", "chain"]]
        .drop_duplicates()
        .merge(
            manifest[["pdb_id", "topology_path"]],
            on="pdb_id",
            how="left",
            validate="many_to_one",
        )
        .sort_values(["pdb_id", "chain"])
    )
    if chain_table["topology_path"].isna().any():
        missing = chain_table.loc[chain_table["topology_path"].isna(), "pdb_id"].unique()
        raise ValueError(f"Missing topology paths for: {missing.tolist()}")

    chain_prediction_path = raw_dir / "isolated_chain_predictions.csv"
    timing_path = raw_dir / "isolated_chain_timings.csv"
    if args.reuse and chain_prediction_path.exists() and timing_path.exists():
        chain_predictions = pd.read_csv(
            chain_prediction_path,
            dtype={"residue_number": str, "insertion_code": str},
        )
        timings = pd.read_csv(timing_path)
    else:
        chain_predictions, timings = predict_isolated_chains(
            chain_table.to_dict("records"),
            checkpoint,
            prepared_dir,
        )
        chain_predictions.to_csv(chain_prediction_path, index=False)
        timings.to_csv(timing_path, index=False)

    failures = timings[timings["status"] != "ok"]
    if not failures.empty:
        raise RuntimeError(
            "Isolated-chain inference failures:\n"
            + failures[["pdb_id", "chain", "error"]].to_string(index=False)
        )
    chain_predictions = normalize_keys(chain_predictions)
    keys = ["pdb_id", "chain", "residue_number", "insertion_code"]
    common = common.rename(
        columns={
            "flexppi_rmsf_angstrom": "full_flexppi_rmsf_angstrom",
            "flexppi_std_angstrom": "full_flexppi_std_angstrom",
        }
    )
    merged = common.merge(
        chain_predictions,
        on=keys,
        how="left",
        validate="one_to_one",
        suffixes=("", "_isolated"),
        indicator=True,
    )
    unmatched = merged[merged["_merge"] != "both"]
    if not unmatched.empty:
        raise ValueError(
            f"{len(unmatched)} benchmark residues lack an isolated-chain prediction"
        )
    merged = merged.drop(columns=["_merge"])
    merged["context_response_angstrom"] = (
        merged["full_flexppi_rmsf_angstrom"]
        - merged["isolated_flexppi_rmsf_angstrom"]
    )
    if merged["is_tcr_pmhc_interface_8A"].dtype != bool:
        merged["is_tcr_pmhc_interface_8A"] = (
            merged["is_tcr_pmhc_interface_8A"]
            .astype(str)
            .str.lower()
            .map({"true": True, "false": False})
            .fillna(False)
        )
    merged.to_csv(raw_dir / "common_context_predictions.csv", index=False)

    chain_metrics = chain_metric_table(merged)
    complex_metrics = complex_metric_table(chain_metrics)
    chain_metrics.to_csv(metrics_dir / "per_chain_metrics.csv", index=False)
    complex_metrics.to_csv(metrics_dir / "per_complex_metrics.csv", index=False)

    overlap = {
        str(value).lower()
        for value in benchmark_report.get("training_overlap_exact_pdb", [])
    }
    external_clean = complex_metrics[~complex_metrics["pdb_id"].isin(overlap)].copy()
    external_clean_chains = chain_metrics[
        ~chain_metrics["pdb_id"].isin(overlap)
    ].copy()
    external_clean_residues = merged[~merged["pdb_id"].isin(overlap)].copy()
    non_peptide_complex_metrics = complex_metric_table(
        external_clean_chains[
            external_clean_chains["chain_role"].astype(str) != "peptide"
        ]
    )
    main_results = condition_results(
        external_clean, args.bootstrap_repetitions, args.seed
    )
    main_results["context_response"] = context_response_summary(
        external_clean_residues,
        args.bootstrap_repetitions,
        args.seed,
    )
    report = {
        "analysis": {
            "name": "ITrace FlexPPI strict zero-shot context ablation",
            "question": (
                "Does the same frozen FlexPPI checkpoint predict bound-state "
                "chain-fitted RMSF better with the full complex than with each "
                "bound-state chain isolated?"
            ),
            "not_tested": (
                "True apo-to-holo or binding-induced delta RMSF; ITrace has no "
                "matched apo structures in this benchmark."
            ),
            "parameter_updates": 0,
            "target_used_during_inference": False,
        },
        "model": {
            "checkpoint": str(checkpoint),
            "checkpoint_sha256": actual_sha,
            "same_checkpoint_as_full_complex_benchmark": True,
            "temperature_K": 310.0,
            "window_ns": 200.0,
            "dataset_code": 1.0,
        },
        "cohort": {
            "requested_complexes": int(len(selected_ids)),
            "chains": int(len(chain_table)),
            "common_residues": int(len(merged)),
            "external_clean_complexes": int(len(external_clean)),
            "excluded_exact_training_overlap": sorted(overlap.intersection(selected_ids)),
        },
        "runtime": {
            "isolated_chain_total_seconds": float(
                timings["end_to_end_seconds"].sum()
            ),
            "isolated_chain_median_seconds": float(
                timings["end_to_end_seconds"].median()
            ),
            "failures": int(len(failures)),
        },
        "external_clean_results": main_results,
        "sensitivity_analysis": {
            "exclude_peptide_chain": condition_results(
                non_peptide_complex_metrics,
                args.bootstrap_repetitions,
                args.seed,
            ),
            "by_chain_role": chain_role_results(
                external_clean_chains,
                args.bootstrap_repetitions,
                args.seed,
            ),
        },
        "interpretation_limits": [
            (
                "Full-minus-isolated is an input-context response, not a causal "
                "binding-induced flexibility change."
            ),
            (
                "The ablation changes partner availability, total sequence length, "
                "and global attention normalization at the same time."
            ),
            (
                "Most full complexes exceed the training maximum length, while "
                "isolated chains are generally length-in-distribution."
            ),
        ],
        "artifacts": {
            "isolated_chain_predictions": str(chain_prediction_path),
            "common_context_predictions": str(
                raw_dir / "common_context_predictions.csv"
            ),
            "per_chain_metrics": str(metrics_dir / "per_chain_metrics.csv"),
            "per_complex_metrics": str(metrics_dir / "per_complex_metrics.csv"),
        },
    }
    safe_report = _json_safe(report)
    json_dump(output_dir / "itrace_context_report.json", safe_report)
    (output_dir / "SUMMARY_ZH.md").write_text(
        markdown_summary(safe_report), encoding="utf-8"
    )
    print(json.dumps(safe_report["cohort"], ensure_ascii=False), flush=True)
    print(f"Wrote {output_dir}", flush=True)


if __name__ == "__main__":
    main()
