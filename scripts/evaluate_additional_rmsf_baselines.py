"""Evaluate BackFlip, ProDy, and RMSF-net with frozen FlexPPI RMSF protocols.

Static benchmarks are macro-averaged over proteins; multi-chain DynaRepo and
ITrace first average chain-level metrics within each complex and then give each
complex equal weight.  This reproduces the estimands used by the existing
paper tables rather than pooling residues from differently sized proteins.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any, Dict, Iterable, List, Mapping, Sequence

import numpy as np
import pandas as pd

from flexppi.comparison import (
    _bootstrap_mean_ci,
    _method_summary,
    _paired_comparison,
    _system_metrics,
)


ROOT = Path(__file__).resolve().parents[1]
BASE = ROOT / "artifacts" / "evaluation" / "additional_baselines"
RAW = BASE / "raw"
OUT = BASE / "results"
STATIC_COMMON = ROOT / "artifacts" / "evaluation" / "fastprotflex_comparison"
ITRACE_COMMON = ROOT / "output" / "itrace_benchmark" / "raw" / "common_residue_predictions.csv"

METHOD_COLUMNS = {
    "fastprotflex": "fastprotflex_score",
    "flexppi": "flexppi_pretrain_rmsf_angstrom",
    "backflip_structure": "backflip_structure_rmsf_angstrom",
    "backflip_seq": "backflip_seq_rmsf_angstrom",
    "prody_anm": "anm_sq_fluct",
    "prody_gnm": "gnm_sq_fluct",
    "rmsf_net_pdb": "rmsf_net_pdb_score",
}
ABSOLUTE_SCALE = {
    "fastprotflex": False,
    "flexppi": True,
    "backflip_structure": True,
    "backflip_seq": True,
    "prody_anm": False,
    "prody_gnm": False,
    "rmsf_net_pdb": False,
}
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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bootstrap-repetitions", type=int, default=2_000)
    parser.add_argument("--seed", type=int, default=2026)
    return parser.parse_args()


def normalize_identity(frame: pd.DataFrame, sample_column: str) -> pd.DataFrame:
    frame = frame.copy()
    frame[sample_column] = frame[sample_column].astype(str).str.lower()
    frame["chain"] = frame["chain"].fillna("_").astype(str).replace("", "_")
    frame["residue_number"] = frame["residue_number"].astype(str).str.replace(
        r"\.0$", "", regex=True
    )
    frame["insertion_code"] = frame["insertion_code"].fillna("").astype(str)
    frame["amino_acid_3"] = frame["amino_acid_3"].astype(str).str.upper()
    return frame


def load_backflip(path: Path, output_column: str) -> pd.DataFrame:
    frame = pd.read_csv(
        path,
        dtype={
            "sample_id": str,
            "author_residue_id": str,
            "residue_number": str,
            "insertion_code": str,
        },
        low_memory=False,
    )
    return frame.rename(columns={"backflip_rmsf_angstrom": output_column})


def load_new_predictions(kind: str) -> pd.DataFrame:
    suffix = "static" if kind == "static" else "itrace"
    structure = load_backflip(
        RAW / f"backflip_structure_{suffix}.csv",
        "backflip_structure_rmsf_angstrom",
    )
    sequence = load_backflip(
        RAW / f"backflip_seq_{suffix}.csv", "backflip_seq_rmsf_angstrom"
    )
    if kind == "static":
        structure_chainwise = load_backflip(
            RAW / "backflip_structure_static_complex_chainwise.csv",
            "backflip_structure_rmsf_angstrom",
        )
        sequence_chainwise = load_backflip(
            RAW / "backflip_seq_static_complex_chainwise.csv",
            "backflip_seq_rmsf_angstrom",
        )
        for name, full, replacement in (
            ("BackFlip-structure", structure, structure_chainwise),
            ("BackFlip-seq", sequence, sequence_chainwise),
        ):
            replacement_ids = set(replacement["sample_id"])
            full.drop(full.index[full["sample_id"].isin(replacement_ids)], inplace=True)
            full.reset_index(drop=True, inplace=True)
            replacement.reset_index(drop=True, inplace=True)
            if not replacement_ids:
                raise ValueError(f"{name} chain-wise replacement is empty")
            if set(replacement["cohort"]) != {"dynarepo_complex"}:
                raise ValueError(f"{name} chain-wise replacement has unexpected cohorts")
            combined = pd.concat([full, replacement], ignore_index=True)
            if name == "BackFlip-structure":
                structure = combined
            else:
                sequence = combined
    prody = pd.read_csv(
        RAW / f"prody_{suffix}.csv",
        dtype={
            "sample_id": str,
            "author_residue_id": str,
            "residue_number": str,
            "insertion_code": str,
        },
        low_memory=False,
    )
    rmsf_net = pd.read_csv(
        RAW / f"rmsf_net_{suffix}.csv",
        dtype={
            "sample_id": str,
            "author_residue_id": str,
            "residue_number": str,
            "insertion_code": str,
        },
        low_memory=False,
    )
    for frame in (structure, sequence, prody, rmsf_net):
        frame["sample_id"] = frame["sample_id"].astype(str)
        frame["chain"] = frame["chain"].fillna("_").astype(str).replace("", "_")
        frame["author_residue_id"] = (
            frame["author_residue_id"]
            .fillna("")
            .astype(str)
            .str.replace(r"\.0$", "", regex=True)
        )
        frame["residue_number"] = (
            frame["residue_number"]
            .fillna("")
            .astype(str)
            .str.replace(r"\.0$", "", regex=True)
        )
        frame["insertion_code"] = frame["insertion_code"].fillna("").astype(str)
        frame["amino_acid_3"] = frame["amino_acid_3"].astype(str).str.upper()
    identity = [
        "sample_id",
        "residue_index",
        "chain",
        "author_residue_id",
        "residue_number",
        "insertion_code",
        "amino_acid_3",
    ]
    stable_key = ["sample_id", "chain", "author_residue_id"]
    for name, frame in (
        ("BackFlip-structure", structure),
        ("RMSF-net-PDB", rmsf_net),
        ("BackFlip-seq", sequence),
        ("ProDy", prody),
    ):
        if frame.duplicated(stable_key).any():
            raise ValueError(f"{name} contains duplicate stable residue identities")
    merged = prody[identity + ["anm_sq_fluct", "gnm_sq_fluct"]].merge(
        structure[stable_key + ["backflip_structure_rmsf_angstrom"]],
        on=stable_key,
        how="inner",
        validate="one_to_one",
    )
    merged = merged.merge(
        sequence[stable_key + ["backflip_seq_rmsf_angstrom"]],
        on=stable_key,
        how="inner",
        validate="one_to_one",
    )
    merged = merged.merge(
        rmsf_net[stable_key + ["rmsf_net_pdb_score"]],
        on=stable_key,
        how="inner",
        validate="one_to_one",
    )
    if merged.empty:
        raise ValueError("The additional baselines have no common residue identities")
    return merged


def metric_bundle(group: pd.DataFrame, target_column: str) -> Dict[str, Dict[str, Any]]:
    required = [target_column, *METHOD_COLUMNS.values()]
    valid = group[required].apply(np.isfinite).all(axis=1)
    selected = group.loc[valid]
    if len(selected) < 2:
        return {}
    target = selected[target_column].to_numpy(float)
    output: Dict[str, Dict[str, Any]] = {}
    for method, column in METHOD_COLUMNS.items():
        metrics = _system_metrics(
            selected[column].to_numpy(float),
            target,
            absolute_scale=ABSOLUTE_SCALE[method],
        )
        if not metrics:
            return {}
        output[method] = metrics
    return output


def average_chain_rows(rows: Sequence[Mapping[str, Any]]) -> List[Dict[str, Any]]:
    grouped: Dict[str, List[Mapping[str, Any]]] = {}
    for row in rows:
        grouped.setdefault(str(row["sample_id"]), []).append(row)
    output: List[Dict[str, Any]] = []
    for sample_id, chain_rows in sorted(grouped.items()):
        result: Dict[str, Any] = {
            "sample_id": sample_id,
            "chain_count": len(chain_rows),
            "residue_count": sum(int(row["residue_count"]) for row in chain_rows),
        }
        for method in METHOD_COLUMNS:
            result[method] = {}
            for metric in METRICS:
                values = [
                    row[method].get(metric)
                    for row in chain_rows
                    if row.get(method, {}).get(metric) is not None
                ]
                if values:
                    result[method][metric] = float(np.mean(values))
        output.append(result)
    return output


def flatten_rows(rows: Sequence[Mapping[str, Any]]) -> pd.DataFrame:
    flat = []
    for row in rows:
        item = {key: value for key, value in row.items() if not isinstance(value, dict)}
        for method in METHOD_COLUMNS:
            for metric, value in row.get(method, {}).items():
                item[f"{method}__{metric}"] = value
        flat.append(item)
    return pd.DataFrame(flat)


def summary_rows(
    cohort: str,
    rows: Sequence[Mapping[str, Any]],
    repetitions: int,
    seed: int,
) -> List[Dict[str, Any]]:
    summaries = {
        method: _method_summary(rows, method, repetitions, seed)
        for method in METHOD_COLUMNS
    }
    output: List[Dict[str, Any]] = []
    for method, metrics in summaries.items():
        for metric, values in metrics.items():
            output.append(
                {
                    "cohort": cohort,
                    "method": method,
                    "metric": metric,
                    **values,
                }
            )
    return output


def paired_rows(
    cohort: str,
    rows: Sequence[Mapping[str, Any]],
    repetitions: int,
    seed: int,
) -> List[Dict[str, Any]]:
    output = []
    for contender in METHOD_COLUMNS:
        if contender == "flexppi":
            continue
        comparisons = _paired_comparison(rows, "flexppi", contender, repetitions, seed)
        for metric, values in comparisons.items():
            output.append(
                {
                    "cohort": cohort,
                    "contender": "flexppi",
                    "reference": contender,
                    "metric": metric,
                    **values,
                }
            )
    return output


def merge_static_predictions() -> pd.DataFrame:
    atlas = pd.read_csv(
        STATIC_COMMON / "atlas_common_residue_predictions.csv",
        dtype={"sample_id": str, "author_residue_id": str},
    )
    atlas["cohort"] = "atlas_single_chain"
    dynarepo = pd.read_csv(
        STATIC_COMMON / "dynarepo" / "dynarepo_common_residue_predictions.csv",
        dtype={"sample_id": str, "author_residue_id": str},
    )
    complex_ids = set(
        pd.read_csv(
            STATIC_COMMON / "chain_context" / "complex_per_system_metrics.csv",
            dtype={"sample_id": str},
        )["sample_id"]
    )
    independent = pd.read_csv(
        STATIC_COMMON
        / "chain_context"
        / "fastprotflex_independent_chain_predictions.csv",
        dtype={
            "sample_id": str,
            "residue_number": str,
            "insertion_code": str,
        },
    )
    independent["original_sample_id"] = independent["sample_id"].str.split(
        "__chain__", n=1
    ).str[0]
    independent["chain"] = independent["chain"].fillna("_").replace("", "_")
    independent["author_residue_id"] = (
        independent["residue_number"].astype(str)
        + independent["insertion_code"].fillna("").astype(str)
    )
    independent = independent.rename(
        columns={"fastprotflex_score": "fastprotflex_independent_chain_score"}
    )
    complex_frame = dynarepo.loc[dynarepo["sample_id"].isin(complex_ids)].merge(
        independent[
            [
                "original_sample_id",
                "chain",
                "author_residue_id",
                "fastprotflex_independent_chain_score",
            ]
        ],
        left_on=["sample_id", "chain", "author_residue_id"],
        right_on=["original_sample_id", "chain", "author_residue_id"],
        how="left",
        validate="one_to_one",
    )
    if complex_frame["fastprotflex_independent_chain_score"].isna().any():
        raise ValueError("Independent-chain FastProtFlex does not cover DynaRepo complexes")
    dynarepo = dynarepo.loc[~dynarepo["sample_id"].isin(complex_ids)].copy()
    complex_frame["fastprotflex_score"] = complex_frame[
        "fastprotflex_independent_chain_score"
    ]
    dynarepo = pd.concat(
        [dynarepo, complex_frame[dynarepo.columns]], ignore_index=True
    )
    dynarepo["cohort"] = np.where(
        dynarepo["sample_id"].isin(complex_ids),
        "dynarepo_complex",
        "dynarepo_single_chain",
    )
    truth = pd.concat([atlas, dynarepo], ignore_index=True)
    predictions = load_new_predictions("static")
    for frame in (truth, predictions):
        frame["sample_id"] = frame["sample_id"].astype(str)
        frame["chain"] = frame["chain"].fillna("_").astype(str).replace("", "_")
        frame["author_residue_id"] = frame["author_residue_id"].astype(str)
    stable_key = ["sample_id", "chain", "author_residue_id"]
    if truth.duplicated(stable_key).any() or predictions.duplicated(stable_key).any():
        raise ValueError("Static data contain duplicate stable residue identities")
    value_columns = [*stable_key, *list(METHOD_COLUMNS.values())[2:]]
    common = truth.merge(
        predictions[value_columns],
        on=stable_key,
        how="left",
        validate="one_to_one",
    )
    missing = common[list(METHOD_COLUMNS.values())[2:]].isna().any(axis=1)
    if missing.any():
        preview = common.loc[missing, stable_key].head(10).to_dict("records")
        raise ValueError(
            f"New baselines miss {int(missing.sum())} frozen static residues: {preview}"
        )
    return common


def evaluate_static(common: pd.DataFrame, repetitions: int, seed: int) -> None:
    all_rows: List[Dict[str, Any]] = []
    summary: List[Dict[str, Any]] = []
    paired: List[Dict[str, Any]] = []
    for cohort in ("atlas_single_chain", "dynarepo_single_chain", "dynarepo_complex"):
        selected = common.loc[common["cohort"] == cohort]
        chain_rows: List[Dict[str, Any]] = []
        grouping = ["sample_id", "chain"] if cohort == "dynarepo_complex" else ["sample_id"]
        for keys, group in selected.groupby(grouping, sort=True):
            metrics = metric_bundle(group, "target_rmsf_angstrom")
            if not metrics:
                continue
            if not isinstance(keys, tuple):
                keys = (keys,)
            chain_rows.append(
                {
                    "sample_id": keys[0],
                    "chain": keys[1] if len(keys) > 1 else "__whole_system__",
                    "cohort": cohort,
                    "residue_count": len(group),
                    **metrics,
                }
            )
        system_rows = (
            average_chain_rows(chain_rows)
            if cohort == "dynarepo_complex"
            else chain_rows
        )
        for row in system_rows:
            row["cohort"] = cohort
        all_rows.extend(system_rows)
        summary.extend(summary_rows(cohort, system_rows, repetitions, seed))
        paired.extend(paired_rows(cohort, system_rows, repetitions, seed))
    common.to_csv(OUT / "static_common_residue_predictions.csv", index=False)
    flatten_rows(all_rows).to_csv(OUT / "static_per_system_metrics.csv", index=False)
    pd.DataFrame(summary).to_csv(OUT / "static_summary_metrics.csv", index=False)
    pd.DataFrame(paired).to_csv(OUT / "static_paired_tests.csv", index=False)


def itrace_subset(frame: pd.DataFrame, subset: str) -> np.ndarray:
    if subset == "overall":
        return np.ones(len(frame), dtype=bool)
    if subset == "interface":
        return frame["is_tcr_pmhc_interface_8A"].fillna(False).to_numpy(bool)
    if subset == "non_interface":
        return ~frame["is_tcr_pmhc_interface_8A"].fillna(False).to_numpy(bool)
    if subset == "tcr":
        return frame["chain_role"].isin(["tcr_alpha", "tcr_beta"]).to_numpy()
    if subset == "pmhc":
        return frame["chain_role"].isin(["mhc", "b2m", "peptide"]).to_numpy()
    if subset in {"peptide", "tcr_alpha", "tcr_beta"}:
        return (frame["chain_role"] == subset).to_numpy()
    if subset == "cdr_all":
        return frame["region_group"].isin(["cdr1", "cdr2", "cdr3"]).to_numpy()
    if subset in {"cdr1", "cdr2", "cdr3", "tcr_framework"}:
        return (frame["region_group"] == subset).to_numpy()
    raise ValueError(subset)


def merge_itrace_predictions() -> pd.DataFrame:
    truth = pd.read_csv(
        ITRACE_COMMON,
        dtype={"pdb_id": str, "residue_number": str, "insertion_code": str},
    )
    truth = normalize_identity(truth, "pdb_id")
    truth = truth.loc[truth["pdb_id"] != "3kxf"].copy()
    truth = truth.rename(
        columns={"flexppi_rmsf_angstrom": "flexppi_pretrain_rmsf_angstrom"}
    )
    predictions = normalize_identity(load_new_predictions("itrace"), "sample_id")
    predictions = predictions.rename(columns={"sample_id": "pdb_id"})
    key = ["pdb_id", "chain", "residue_number", "insertion_code", "amino_acid_3"]
    value_columns = [*list(METHOD_COLUMNS.values())[2:]]
    if predictions.duplicated(key).any():
        raise ValueError("ITrace baseline predictions contain duplicate residue identities")
    common = truth.merge(
        predictions[key + value_columns],
        on=key,
        how="inner",
        validate="one_to_one",
    )
    if len(common) != len(truth):
        missing = len(truth) - len(common)
        raise ValueError(
            f"ITrace new baselines do not cover the frozen common set: "
            f"truth={len(truth)}, common={len(common)}, missing={missing}"
        )
    return common


def evaluate_itrace(common: pd.DataFrame, repetitions: int, seed: int) -> None:
    subsets = (
        "overall",
        "interface",
        "non_interface",
        "tcr",
        "pmhc",
        "peptide",
        "tcr_alpha",
        "tcr_beta",
        "cdr_all",
        "cdr1",
        "cdr2",
        "cdr3",
        "tcr_framework",
    )
    all_rows: List[Dict[str, Any]] = []
    summary: List[Dict[str, Any]] = []
    paired: List[Dict[str, Any]] = []
    for subset in subsets:
        selected = common.loc[itrace_subset(common, subset)]
        chain_rows: List[Dict[str, Any]] = []
        for (pdb_id, chain), group in selected.groupby(["pdb_id", "chain"], sort=True):
            if len(group) < 4:
                continue
            metrics = metric_bundle(group, "target_chainfit_rmsf_angstrom")
            if not metrics:
                continue
            chain_rows.append(
                {
                    "sample_id": pdb_id,
                    "chain": chain,
                    "subset": subset,
                    "residue_count": len(group),
                    **metrics,
                }
            )
        complex_rows = average_chain_rows(chain_rows)
        for row in complex_rows:
            row["subset"] = subset
        all_rows.extend(complex_rows)
        summary.extend(summary_rows(subset, complex_rows, repetitions, seed))
        paired.extend(paired_rows(subset, complex_rows, repetitions, seed))
    common.to_csv(OUT / "itrace_common_residue_predictions.csv", index=False)
    flatten_rows(all_rows).to_csv(OUT / "itrace_per_complex_metrics.csv", index=False)
    pd.DataFrame(summary).rename(columns={"cohort": "subset"}).to_csv(
        OUT / "itrace_summary_metrics.csv", index=False
    )
    pd.DataFrame(paired).rename(columns={"cohort": "subset"}).to_csv(
        OUT / "itrace_paired_tests.csv", index=False
    )


def timing_report() -> Dict[str, Any]:
    report: Dict[str, Any] = {}
    for cohort in ("static", "itrace"):
        report[cohort] = {}
        for method, path, column in (
            (
                "backflip_structure",
                RAW / f"backflip_structure_{cohort}.csv",
                "inference_seconds",
            ),
            ("backflip_seq", RAW / f"backflip_seq_{cohort}.csv", "inference_seconds"),
            ("prody_anm", RAW / f"prody_{cohort}.csv", "anm_seconds"),
            (
                "rmsf_net_pdb",
                RAW / f"rmsf_net_{cohort}.csv",
                "__rmsf_net_total__",
            ),
            ("prody_gnm", RAW / f"prody_{cohort}.csv", "gnm_seconds"),
        ):
            frame = pd.read_csv(path, low_memory=False)
            if column == "__rmsf_net_total__":
                first = frame.groupby("sample_id")[
                    ["preprocessing_seconds", "inference_seconds"]
                ].first()
                per_sample = first.sum(axis=1).to_numpy(float)
            else:
                per_sample = frame.groupby("sample_id")[column].first().to_numpy(float)
            report[cohort][method] = {
                "sample_count": int(len(per_sample)),
                "mean_seconds": float(np.mean(per_sample)),
                "median_seconds": float(np.median(per_sample)),
                "p95_seconds": float(np.quantile(per_sample, 0.95)),
                "total_inference_seconds": float(np.sum(per_sample)),
            }
    return report


def main() -> None:
    args = parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    static = merge_static_predictions()
    evaluate_static(static, args.bootstrap_repetitions, args.seed)
    itrace = merge_itrace_predictions()
    evaluate_itrace(itrace, args.bootstrap_repetitions, args.seed)
    with (OUT / "runtime_summary.json").open("w", encoding="utf-8") as handle:
        json.dump(timing_report(), handle, ensure_ascii=False, indent=2)
    print(
        f"Evaluation complete: static residues={len(static)}, "
        f"ITrace residues={len(itrace)}; outputs={OUT}"
    )


if __name__ == "__main__":
    main()
