"""Combine the ThermoPCD pilot with existing DynaRepo/Zenodo bias audits.

All combinations are descriptive PDB-equal arithmetic, not model retraining.
When the same PDB occurs in two sources (2XWT), the two protocol estimates are
averaged first so that the structure is counted once.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


METRICS = [
    "mean_delta_rmsf",
    "fraction_stiffening",
    "fraction_robust_stiffening",
    "fraction_robust_softening",
]


def load_sources(prior_path: Path, thermo_path: Path) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    prior = pd.read_csv(prior_path)
    prior = prior.loc[
        (prior["scenario"].isin(["DynaRepo-current", "Zenodo-6638504"]))
        & (prior["region"] == "interface_8a")
    ].copy()
    prior = prior[
        ["dataset", "pdb_id", "mean_delta_rmsf", "fraction_stiffening", "fraction_robust_stiffening", "fraction_robust_softening"]
    ].drop_duplicates()
    thermo = pd.read_csv(thermo_path)
    thermo = thermo.loc[thermo["region"] == "interface_calpha8"].copy()
    thermo = thermo.rename(columns={"mean_delta_rmsf_angstrom": "mean_delta_rmsf"})
    thermo["dataset"] = "ThermoPCD-refit-pilot"
    thermo = thermo[
        ["dataset", "pdb_id", "mean_delta_rmsf", "fraction_stiffening", "fraction_robust_stiffening", "fraction_robust_softening"]
    ]
    dyna = prior.loc[prior["dataset"] == "DynaRepo-current"].copy()
    zenodo = prior.loc[prior["dataset"] == "Zenodo-6638504"].copy()
    return dyna, zenodo, thermo


def summarize_scenario(name: str, frames: list[pd.DataFrame]) -> dict[str, float | int | str]:
    stacked = pd.concat(frames, ignore_index=True)
    by_pdb = stacked.groupby("pdb_id", as_index=False)[METRICS].mean()
    row: dict[str, float | int | str] = {
        "scenario": name,
        "unique_pdbs": int(len(by_pdb)),
        "source_rows_before_pdb_deduplication": int(len(stacked)),
    }
    for metric in METRICS:
        row[metric] = float(by_pdb[metric].mean())
    return row


def make_figure(scenarios: pd.DataFrame, thermo: pd.DataFrame, output: Path) -> None:
    plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(1, 3, figsize=(13, 4.1), constrained_layout=True)

    order = thermo.sort_values("mean_delta_rmsf")["pdb_id"].tolist()
    selected = thermo.set_index("pdb_id").loc[order]
    colors = ["#D95F02" if value < 0 else "#1B9E77" for value in selected["mean_delta_rmsf"]]
    axes[0].barh(order, selected["mean_delta_rmsf"], color=colors)
    axes[0].axvline(0, color="black", lw=0.8)
    axes[0].set_xlabel("Interface mean ΔRMSF (Å)")
    axes[0].set_title("ThermoPCD target-chain refit")

    x = np.arange(len(order))
    axes[1].bar(x - 0.18, selected["fraction_stiffening"], 0.36, label="raw stiffening", color="#E45756")
    axes[1].bar(
        x + 0.18,
        selected["fraction_robust_softening"],
        0.36,
        label="robust softening",
        color="#4C78A8",
    )
    axes[1].set_xticks(x, order, rotation=45, ha="right")
    axes[1].set_ylim(0, 1.05)
    axes[1].set_ylabel("Interface residue fraction")
    axes[1].set_title("Sign balance within each system")
    axes[1].legend(frameon=False, fontsize=8)

    labels = scenarios["scenario"].tolist()
    sx = np.arange(len(labels))
    axes[2].bar(sx - 0.2, scenarios["fraction_stiffening"], 0.4, color="#E45756", label="raw stiffening")
    axes[2].bar(
        sx + 0.2,
        scenarios["fraction_robust_softening"],
        0.4,
        color="#4C78A8",
        label="robust softening",
    )
    axes[2].set_xticks(sx, labels, rotation=35, ha="right")
    axes[2].set_ylim(0, 1.05)
    axes[2].set_title("PDB-equal descriptive mixtures")
    axes[2].legend(frameon=False, fontsize=8)
    fig.suptitle("ThermoPCD does not by itself rebalance contact-associated rigidification", fontsize=12)
    fig.savefig(output.with_suffix(".png"), dpi=240)
    fig.savefig(output.with_suffix(".pdf"))
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--prior-per-pdb",
        type=Path,
        default=Path("artifacts/training_bias_with_zenodo_6638504/per_pdb_bias.csv"),
    )
    parser.add_argument(
        "--thermo-summary",
        type=Path,
        default=Path("artifacts/thermopcd_audit/refit_pair_region_summary.csv"),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/thermopcd_audit"))
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    dyna, zenodo, thermo = load_sources(args.prior_per_pdb, args.thermo_summary)
    rows = [
        summarize_scenario("DynaRepo current", [dyna]),
        summarize_scenario("DynaRepo + Zenodo", [dyna, zenodo]),
        summarize_scenario("ThermoPCD refit pilot", [thermo]),
        summarize_scenario("DynaRepo + ThermoPCD", [dyna, thermo]),
        summarize_scenario("DynaRepo + Zenodo + ThermoPCD", [dyna, zenodo, thermo]),
    ]
    scenarios = pd.DataFrame(rows)
    scenarios.to_csv(args.output_dir / "training_bias_scenarios.csv", index=False)

    overlap_ids = sorted(set(dyna["pdb_id"]) & set(thermo["pdb_id"]))
    overlap = pd.concat(
        [dyna.loc[dyna["pdb_id"].isin(overlap_ids)], thermo.loc[thermo["pdb_id"].isin(overlap_ids)]],
        ignore_index=True,
    )
    overlap.to_csv(args.output_dir / "cross_protocol_exact_pdb_overlap.csv", index=False)
    make_figure(scenarios, thermo, args.output_dir / "thermopcd_bias_audit")
    report = {
        "thermopcd_pilot_pdbs": sorted(thermo["pdb_id"].tolist()),
        "exact_pdb_overlap_with_dynarepo": overlap_ids,
        "new_pdbs_beyond_dynarepo_and_zenodo": sorted(
            set(thermo["pdb_id"]) - set(dyna["pdb_id"]) - set(zenodo["pdb_id"])
        ),
        "combination_rule": "average source estimates within duplicate PDB, then average PDBs",
        "scenarios": rows,
    }
    (args.output_dir / "training_bias_scenarios.json").write_text(
        json.dumps(report, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    print(json.dumps(report, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
