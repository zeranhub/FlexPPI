"""Expose the residue-level biology hidden in the 7L1D/7RRG CSV outputs.

Both structures contain the same ALHGGWTTK peptide.  The comparison is between
the fixed-coordinate pMHC subset and the corresponding complete TCR-pMHC input.
This script adds heavy-atom TCR contacts, explicitly masks undefined H_chi for
alanine/glycine in interpreted fields, and writes a compact reproducibility
table for manuscript figures and text.
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Any, Dict, List

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from analyze_itrace_residue_biology import (
    CHI_UNDEFINED_AA,
    CONTACT_CUTOFF_ANGSTROM,
    contact_annotations_for_complex,
)


METRICS = ("rmsf", "chi_entropy", "sasa")
PDB_IDS = ("7L1D", "7RRG")
CHAIN_ROLES = {
    "A": "mhc",
    "B": "b2m",
    "C": "peptide",
    "D": "tcr_alpha",
    "E": "tcr_beta",
}


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-dir",
        default=(
            "output/all_experiment_results/01_raw_evaluation_outputs/"
            "7L1D_7RRG_tcr_pmhc_flexppi"
        ),
    )
    parser.add_argument(
        "--output-dir", default="paper/data/itrace_residue_biology"
    )
    parser.add_argument(
        "--itrace-residue-table",
        default=(
            "output/itrace_residue_biology/raw/"
            "matched_residue_context_responses.csv"
        ),
    )
    return parser.parse_args()


def json_safe(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(key): json_safe(item) for key, item in value.items()}
    if isinstance(value, list):
        return [json_safe(item) for item in value]
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return None if not np.isfinite(value) else float(value)
    if isinstance(value, (np.bool_,)):
        return bool(value)
    return value


def main() -> None:
    args = parse_arguments()
    source_dir = Path(args.source_dir)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    itrace_table = Path(args.itrace_residue_table)
    if itrace_table.exists():
        all_residues = pd.read_csv(
            itrace_table, dtype={"author_residue_id": str}, low_memory=False
        )
        all_residues["pdb_id"] = all_residues["pdb_id"].str.upper()
        residues = all_residues[
            all_residues["pdb_id"].isin(PDB_IDS)
            & (all_residues["chain_role"] == "peptide")
        ].copy()
        residues["amino_acid"] = residues["amino_acid_1"]
        residues["peptide_position"] = residues["peptide_position_n"].astype(int)
        residues["direct_tcr_contact_5A"] = (
            residues["min_distance_to_tcr_angstrom"] <= CONTACT_CUTOFF_ANGSTROM
        )
        role_map = {
            (str(row.pdb_id), str(row.chain)): str(row.chain_role)
            for row in all_residues[["pdb_id", "chain", "chain_role"]]
            .drop_duplicates()
            .itertuples(index=False)
        }
        residues["nearest_tcr_role"] = [
            role_map.get((str(pdb_id), str(chain)), "unknown")
            for pdb_id, chain in zip(
                residues["pdb_id"], residues["nearest_tcr_chain"]
            )
        ]
        for metric in METRICS:
            residues["pmhc_before_" + metric] = residues["reference_" + metric]
            residues["tcr_pmhc_after_" + metric] = residues["full_" + metric]
            residues["delta_after_minus_before_" + metric] = residues[
                "binding_context_delta_" + metric
            ]
            residues["delta_over_combined_std_" + metric] = residues[
                "binding_context_z_" + metric
            ]
        residues["interpretable_delta_chi_entropy"] = residues[
            "delta_after_minus_before_chi_entropy"
        ].where(residues["chi_defined"])
        frames = [residues]
    else:
        frames = []
    for pdb_id in (() if frames else PDB_IDS):
        comparison_path = (
            source_dir
            / "analysis"
            / f"{pdb_id}_pMHC_before_vs_TCR_pMHC_after_residues.csv"
        )
        topology_path = (
            source_dir
            / "inputs"
            / "prepared"
            / f"{pdb_id}_biological_assembly_merged.pdb"
        )
        comparison = pd.read_csv(
            comparison_path, dtype={"author_residue_id": str}
        )
        comparison = comparison[comparison["chain"] == "C"].copy()
        comparison["pdb_id"] = pdb_id
        comparison["peptide_position"] = np.arange(1, len(comparison) + 1)
        comparison["chi_defined"] = ~comparison["amino_acid"].isin(
            CHI_UNDEFINED_AA
        )
        contacts = contact_annotations_for_complex(
            pdb_id.lower(), topology_path, CHAIN_ROLES
        )
        contacts = contacts[contacts["chain"] == "C"].copy()
        contacts["pdb_id"] = pdb_id
        comparison = comparison.merge(
            contacts,
            on=["pdb_id", "chain", "author_residue_id"],
            how="left",
            validate="one_to_one",
        )
        comparison["direct_tcr_contact_5A"] = (
            comparison["min_distance_to_tcr_angstrom"]
            <= CONTACT_CUTOFF_ANGSTROM
        )
        comparison["nearest_tcr_role"] = comparison["nearest_tcr_chain"].map(
            CHAIN_ROLES
        )
        comparison["interpretable_delta_chi_entropy"] = comparison[
            "delta_after_minus_before_chi_entropy"
        ].where(comparison["chi_defined"])
        frames.append(comparison)
    residues = pd.concat(frames, ignore_index=True)
    residues.to_csv(
        output_dir / "7L1D_7RRG_peptide_residue_context.csv", index=False
    )

    paired_rows: List[Dict[str, Any]] = []
    for position, group in residues.groupby("peptide_position", sort=True):
        row: Dict[str, Any] = {
            "peptide_position": int(position),
            "amino_acid": str(group["amino_acid"].iloc[0]),
            "chi_defined": bool(group["chi_defined"].iloc[0]),
            "mean_min_distance_to_tcr_angstrom": float(
                group["min_distance_to_tcr_angstrom"].mean()
            ),
            "contacted_complex_count_5A": int(group["direct_tcr_contact_5A"].sum()),
        }
        for metric in METRICS:
            values = group["delta_after_minus_before_" + metric].to_numpy(dtype=float)
            if metric == "chi_entropy" and not row["chi_defined"]:
                row["mean_delta_" + metric] = math.nan
            else:
                row["mean_delta_" + metric] = float(np.mean(values))
            for _, source in group.iterrows():
                row[str(source["pdb_id"]) + "_delta_" + metric] = float(
                    source["delta_after_minus_before_" + metric]
                )
        paired_rows.append(row)
    paired = pd.DataFrame(paired_rows)
    paired.to_csv(output_dir / "7L1D_7RRG_position_replication.csv", index=False)

    correlation: Dict[str, Any] = {}
    pivot = {
        pdb_id: residues[residues["pdb_id"] == pdb_id]
        .sort_values("peptide_position")
        .reset_index(drop=True)
        for pdb_id in PDB_IDS
    }
    for metric in METRICS:
        first = pivot[PDB_IDS[0]]
        second = pivot[PDB_IDS[1]]
        mask = np.ones(len(first), dtype=bool)
        if metric == "chi_entropy":
            mask = first["chi_defined"].to_numpy(dtype=bool)
        result = spearmanr(
            first.loc[mask, "delta_after_minus_before_" + metric],
            second.loc[mask, "delta_after_minus_before_" + metric],
        )
        correlation[metric] = {
            "position_count": int(mask.sum()),
            "spearman_rho": float(result.statistic),
            "pvalue": float(result.pvalue),
        }

    central = residues[residues["peptide_position"].between(3, 8)]
    terminals = residues[residues["peptide_position"].isin([1, 2, 9])]
    def group_mean(group: pd.DataFrame, metric: str) -> float:
        selected = group
        if metric == "chi_entropy":
            selected = selected[selected["chi_defined"]]
        return float(selected["delta_after_minus_before_" + metric].mean())

    summary = {
        "comparison": (
            "complete TCR-pMHC minus pMHC subset at identical retained coordinates"
        ),
        "peptide_sequence": "ALHGGWTTK",
        "complexes": list(PDB_IDS),
        "chi_entropy_rule": (
            "Ala P1 and Gly P4/P5 are excluded from H_chi interpretation because "
            "their training labels have no defined chi torsion."
        ),
        "position_pattern_replication": correlation,
        "central_P3_P8_mean_delta": {
            metric: group_mean(central, metric)
            for metric in METRICS
        },
        "terminal_P1_P2_P9_mean_delta": {
            metric: group_mean(terminals, metric)
            for metric in METRICS
        },
        "direct_contact_counts": {
            pdb_id: int(
                residues.loc[residues["pdb_id"] == pdb_id, "direct_tcr_contact_5A"].sum()
            )
            for pdb_id in PDB_IDS
        },
    }
    (output_dir / "7L1D_7RRG_residue_biology_summary.json").write_text(
        json.dumps(json_safe(summary), indent=2, ensure_ascii=False),
        encoding="utf-8",
    )


if __name__ == "__main__":
    main()
