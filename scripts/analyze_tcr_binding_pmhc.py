"""Compare pMHC before TCR binding with the complete TCR-pMHC assembly.

The comparison keeps the experimental coordinates fixed:

* before binding: chains A+B+C (HLA-A*03:01, beta-2-microglobulin, peptide)
* after binding: chains A+B+C+D+E (the complete Biological Assembly 1)

FlexPPI predicts residue-level RMSF, side-chain entropy and SASA.  It does not
generate a new set of coordinates.  The output PDB files store predicted RMSF
or after-minus-before changes in the B-factor column for PyMOL colouring.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Dict, List, Mapping, Sequence, Tuple

import numpy as np
import torch

from analyze_tcr_pmhc_context import (
    METRICS,
    _metric_maps,
    _prediction_map,
    _range,
    _residue_key,
    _sha256,
    _write_bfactor_pdb,
    _write_csv,
    _write_json,
    extract_chains,
    predict_state,
)
from flexppi.inference import load_model


ResidueKey = Tuple[str, str]
PMHC_CHAINS = ("A", "B", "C")


def build_binding_rows(
    pdb_id: str,
    pmhc: Mapping[str, Any],
    whole: Mapping[str, Any],
) -> List[Dict[str, Any]]:
    """Create matched A+B+C rows for after-minus-before comparisons."""

    before_map = _prediction_map(pmhc)
    after_map = _prediction_map(whole)
    keys = sorted(
        set(before_map) & set(after_map),
        key=lambda item: (item[0], int("".join(ch for ch in item[1] if ch.isdigit()) or 0), item[1]),
    )
    rows: List[Dict[str, Any]] = []
    for key in keys:
        if key[0] not in PMHC_CHAINS:
            continue
        before = before_map[key]
        after = after_map[key]
        row: Dict[str, Any] = {
            "pdb_id": pdb_id,
            "chain": key[0],
            "chain_role": before["chain_role"],
            "author_residue_id": key[1],
            "amino_acid": before["amino_acid"],
        }
        for metric in METRICS:
            before_value = float(before[metric])
            after_value = float(after[metric])
            before_std = float(before[metric + "_std"])
            after_std = float(after[metric + "_std"])
            combined_std = math.sqrt(before_std**2 + after_std**2)
            delta = after_value - before_value
            row["pmhc_before_" + metric] = before_value
            row["pmhc_before_" + metric + "_std"] = before_std
            row["tcr_pmhc_after_" + metric] = after_value
            row["tcr_pmhc_after_" + metric + "_std"] = after_std
            row["delta_after_minus_before_" + metric] = delta
            row["delta_over_combined_std_" + metric] = (
                delta / combined_std if combined_std else math.nan
            )
        rows.append(row)
    return rows


def numeric_summary(values: Sequence[float]) -> Dict[str, Any]:
    array = np.asarray(values, dtype=np.float64)
    array = array[np.isfinite(array)]
    if not len(array):
        return {"n": 0}
    return {
        "n": int(len(array)),
        "mean": float(array.mean()),
        "median": float(np.median(array)),
        "minimum": float(array.min()),
        "maximum": float(array.max()),
    }


def summarize_binding(rows: Sequence[Mapping[str, Any]]) -> Dict[str, Any]:
    groups = {
        "pmhc": list(rows),
        "hla_heavy_chain": [row for row in rows if row["chain"] == "A"],
        "beta2m": [row for row in rows if row["chain"] == "B"],
        "peptide": [row for row in rows if row["chain"] == "C"],
    }
    result: Dict[str, Any] = {"delta_definition": "TCR-pMHC after - pMHC before", "groups": {}}
    for group_name, selected in groups.items():
        group: Dict[str, Any] = {"residue_count": len(selected)}
        for metric in METRICS:
            before = [float(row["pmhc_before_" + metric]) for row in selected]
            after = [float(row["tcr_pmhc_after_" + metric]) for row in selected]
            delta = [float(row["delta_after_minus_before_" + metric]) for row in selected]
            group[metric] = {
                "pmhc_before": numeric_summary(before),
                "tcr_pmhc_after": numeric_summary(after),
                "delta_after_minus_before": numeric_summary(delta),
            }
        result["groups"][group_name] = group
    return result


def cross_peptide_rows(
    rows_by_id: Mapping[str, Sequence[Mapping[str, Any]]]
) -> List[Dict[str, Any]]:
    identifiers = list(rows_by_id)
    maps = {
        pdb_id: {
            (str(row["chain"]), str(row["author_residue_id"])): row
            for row in rows
            if row["chain"] == "C"
        }
        for pdb_id, rows in rows_by_id.items()
    }
    keys = sorted(
        set(maps[identifiers[0]]) & set(maps[identifiers[1]]),
        key=lambda item: int("".join(ch for ch in item[1] if ch.isdigit()) or 0),
    )
    result: List[Dict[str, Any]] = []
    for key in keys:
        row: Dict[str, Any] = {
            "chain": key[0],
            "author_residue_id": key[1],
            "amino_acid": maps[identifiers[0]][key]["amino_acid"],
        }
        for pdb_id in identifiers:
            source = maps[pdb_id][key]
            for metric in METRICS:
                for prefix in (
                    "pmhc_before_",
                    "tcr_pmhc_after_",
                    "delta_after_minus_before_",
                    "delta_over_combined_std_",
                ):
                    field = prefix + metric
                    row[pdb_id + "_" + field] = float(source[field])
        result.append(row)
    return result


def write_binding_pymol_script(
    pdb_id: str,
    before_pdb: Path,
    after_pdb: Path,
    delta_pdb: Path,
    output_path: Path,
    session_path: Path,
    before_range: Tuple[float, float],
    after_range: Tuple[float, float],
    delta_limit: float,
) -> None:
    shared_min = min(before_range[0], after_range[0])
    shared_max = max(before_range[1], after_range[1])
    prefix = "PDB" + pdb_id
    script = f'''# {pdb_id}: pMHC before TCR binding versus complete TCR-pMHC after binding
# Coordinates are experimental and fixed; B-factor stores FlexPPI RMSF in Angstrom.
reinitialize
bg_color white
set ray_opaque_background, off
set antialias, 2
set cartoon_fancy_helices, 1
set cartoon_smooth_loops, 1
load "{before_pdb.resolve().as_posix()}", {prefix}_pMHC_before
load "{after_pdb.resolve().as_posix()}", {prefix}_TCR_pMHC_after
load "{delta_pdb.resolve().as_posix()}", {prefix}_binding_delta
hide everything, all
show cartoon, {prefix}_pMHC_before or {prefix}_TCR_pMHC_after or {prefix}_binding_delta
show sticks, ({prefix}_pMHC_before or {prefix}_TCR_pMHC_after or {prefix}_binding_delta) and chain C
show sticks, byres ({prefix}_TCR_pMHC_after and chain D+E within 5 of ({prefix}_TCR_pMHC_after and chain C))
spectrum b, blue_white_red, {prefix}_pMHC_before, minimum={shared_min:.4f}, maximum={shared_max:.4f}
spectrum b, blue_white_red, {prefix}_TCR_pMHC_after, minimum={shared_min:.4f}, maximum={shared_max:.4f}
spectrum b, blue_white_red, {prefix}_binding_delta, minimum={-delta_limit:.4f}, maximum={delta_limit:.4f}
select {prefix}_peptide_before, {prefix}_pMHC_before and chain C
select {prefix}_peptide_after, {prefix}_TCR_pMHC_after and chain C
select {prefix}_pTrp6_before, {prefix}_pMHC_before and chain C and resi 6
select {prefix}_pTrp6_after, {prefix}_TCR_pMHC_after and chain C and resi 6
set grid_mode, 1
set grid_slot, 1, {prefix}_pMHC_before
set grid_slot, 2, {prefix}_TCR_pMHC_after
set grid_slot, 3, {prefix}_binding_delta
orient {prefix}_pMHC_before and chain A+C
zoom all, 4
set seq_view, 1
scene {prefix}_TCR_binding_comparison, store
save "{session_path.resolve().as_posix()}"
'''
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(script, encoding="utf-8")


def write_report(
    path: Path,
    checkpoint: Path,
    summaries: Mapping[str, Any],
    peptide_rows: Sequence[Mapping[str, Any]],
) -> None:
    lines = [
        "# Corrected FlexPPI comparison: pMHC before versus TCR-pMHC after binding",
        "",
        "Before = A+B+C pMHC. After = A+B+C+D+E TCR-pMHC. Coordinates shared by "
        "the two states are unchanged. Delta = after - before; negative RMSF delta means "
        "the model predicts that TCR binding rigidifies the residue.",
        "",
        f"Checkpoint: `{checkpoint.resolve()}`",
        f"Checkpoint SHA256: `{_sha256(checkpoint)}`",
        "",
        "| PDB | peptide RMSF before (A) | peptide RMSF after (A) | mean delta (A) |",
        "|---|---:|---:|---:|",
    ]
    for pdb_id, summary in summaries.items():
        peptide = summary["groups"]["peptide"]["rmsf"]
        lines.append(
            f"| {pdb_id} | {peptide['pmhc_before']['mean']:.3f} | "
            f"{peptide['tcr_pmhc_after']['mean']:.3f} | "
            f"{peptide['delta_after_minus_before']['mean']:.3f} |"
        )
    lines.extend([
        "",
        "## Per-peptide-residue RMSF",
        "",
        "| pos | aa | 7L1D before | 7L1D after | 7L1D delta | 7RRG before | 7RRG after | 7RRG delta |",
        "|---:|:---:|---:|---:|---:|---:|---:|---:|",
    ])
    for row in peptide_rows:
        lines.append(
            f"| {row['author_residue_id']} | {row['amino_acid']} | "
            f"{row['7L1D_pmhc_before_rmsf']:.3f} | "
            f"{row['7L1D_tcr_pmhc_after_rmsf']:.3f} | "
            f"{row['7L1D_delta_after_minus_before_rmsf']:.3f} | "
            f"{row['7RRG_pmhc_before_rmsf']:.3f} | "
            f"{row['7RRG_tcr_pmhc_after_rmsf']:.3f} | "
            f"{row['7RRG_delta_after_minus_before_rmsf']:.3f} |"
        )
    lines.extend([
        "",
        "This is the binding comparison. The previous C+D+E versus A-E analysis is a separate "
        "MHC-removal context ablation and must not be labelled as before-versus-after TCR binding.",
        "",
    ])
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output-dir",
        default="artifacts/evaluation/7L1D_7RRG_tcr_pmhc_flexppi",
    )
    parser.add_argument("--checkpoint", default="artifacts/pretrain/best.pt")
    parser.add_argument("--device", default="auto")
    parser.add_argument("--temperature", type=float, default=300.0)
    parser.add_argument("--window-ns", type=float, default=100.0)
    args = parser.parse_args()

    output_dir = Path(args.output_dir)
    checkpoint = Path(args.checkpoint)
    prepared_dir = output_dir / "inputs" / "prepared"
    prediction_dir = output_dir / "predictions"
    analysis_dir = output_dir / "analysis"
    pymol_dir = output_dir / "pymol"

    whole_paths = {
        pdb_id: prepared_dir / f"{pdb_id}_biological_assembly_merged.pdb"
        for pdb_id in ("7L1D", "7RRG")
    }
    pmhc_paths: Dict[str, Path] = {}
    preparation: Dict[str, Any] = {}
    for pdb_id, whole_path in whole_paths.items():
        if not whole_path.exists():
            raise FileNotFoundError(
                f"Prepared complete assembly is missing: {whole_path}. "
                "Run analyze_tcr_pmhc_context.py first."
            )
        pmhc_path = prepared_dir / f"{pdb_id}_pMHC_ABC_before_TCR_binding.pdb"
        preparation[pdb_id] = extract_chains(whole_path, pmhc_path, PMHC_CHAINS)
        pmhc_paths[pdb_id] = pmhc_path

    model, model_config, device = load_model(str(checkpoint), args.device)
    predictions: Dict[str, Dict[str, Any]] = {}
    rows_by_id: Dict[str, List[Dict[str, Any]]] = {}
    summaries: Dict[str, Any] = {}
    pymol_artifacts: Dict[str, Any] = {}

    for pdb_id in whole_paths:
        before = predict_state(
            model, model_config, device, pmhc_paths[pdb_id], args.temperature, args.window_ns
        )
        after = predict_state(
            model, model_config, device, whole_paths[pdb_id], args.temperature, args.window_ns
        )
        predictions[pdb_id] = {"pmhc_before": before, "tcr_pmhc_after": after}
        _write_json(prediction_dir / pdb_id / "pmhc_before" / "prediction.json", before)
        _write_csv(
            prediction_dir / pdb_id / "pmhc_before" / "residue_predictions.csv",
            before["rows"],
        )
        # Preserve a self-contained copy of the after-binding output for this comparison.
        _write_json(prediction_dir / pdb_id / "tcr_pmhc_after" / "prediction.json", after)
        _write_csv(
            prediction_dir / pdb_id / "tcr_pmhc_after" / "residue_predictions.csv",
            after["rows"],
        )

        rows = build_binding_rows(pdb_id, before, after)
        rows_by_id[pdb_id] = rows
        summaries[pdb_id] = summarize_binding(rows)
        _write_csv(analysis_dir / f"{pdb_id}_pMHC_before_vs_TCR_pMHC_after_residues.csv", rows)
        _write_json(
            analysis_dir / f"{pdb_id}_pMHC_before_vs_TCR_pMHC_after_summary.json",
            summaries[pdb_id],
        )

        before_rmsf = _metric_maps(before, "rmsf")
        after_rmsf = _metric_maps(after, "rmsf")
        delta_rmsf: Dict[ResidueKey, float] = {
            _residue_key(row): float(row["delta_after_minus_before_rmsf"])
            for row in rows
        }
        visual_dir = prediction_dir / pdb_id / "pymol_values"
        before_visual = visual_dir / f"{pdb_id}_pMHC_before_FlexPPI_RMSF_A.pdb"
        after_visual = visual_dir / f"{pdb_id}_TCR_pMHC_after_FlexPPI_RMSF_A.pdb"
        delta_visual = visual_dir / f"{pdb_id}_TCR_binding_after_minus_before_RMSF_A.pdb"
        _write_bfactor_pdb(
            pmhc_paths[pdb_id], before_visual, before_rmsf,
            ("B-factor = FlexPPI RMSF; pMHC before TCR binding",),
        )
        _write_bfactor_pdb(
            whole_paths[pdb_id], after_visual, after_rmsf,
            ("B-factor = FlexPPI RMSF; complete TCR-pMHC after binding",),
        )
        _write_bfactor_pdb(
            pmhc_paths[pdb_id], delta_visual, delta_rmsf,
            ("B-factor = TCR-pMHC after minus pMHC before RMSF",),
        )
        delta_limit = max(abs(value) for value in delta_rmsf.values()) or 1.0
        pml = pymol_dir / f"compare_{pdb_id}_pMHC_before_vs_TCR_pMHC_after.pml"
        pse = pymol_dir / f"compare_{pdb_id}_pMHC_before_vs_TCR_pMHC_after.pse"
        write_binding_pymol_script(
            pdb_id,
            before_visual,
            after_visual,
            delta_visual,
            pml,
            pse,
            _range(before_rmsf.values()),
            _range(after_rmsf.values()),
            delta_limit,
        )
        pymol_artifacts[pdb_id] = {
            "comparison_pml": str(pml.resolve()),
            "session_after_running": str(pse.resolve()),
            "pmhc_before_rmsf_pdb": str(before_visual.resolve()),
            "tcr_pmhc_after_rmsf_pdb": str(after_visual.resolve()),
            "binding_delta_pdb": str(delta_visual.resolve()),
        }

    peptide_rows = cross_peptide_rows(rows_by_id)
    peptide_csv = analysis_dir / "7L1D_7RRG_pMHC_before_vs_TCR_pMHC_after_peptide_RMSF.csv"
    _write_csv(peptide_csv, peptide_rows)
    report_path = output_dir / "REPORT_TCR_BINDING_CORRECTED.md"
    write_report(report_path, checkpoint, summaries, peptide_rows)

    metadata_path = output_dir / "tcr_binding_metadata.json"
    metadata = {
        "analysis": "correct binding comparison: pMHC before versus TCR-pMHC after",
        "state_definition": {
            "before": "A+B+C pMHC",
            "after": "A+B+C+D+E complete TCR-pMHC Biological Assembly 1",
            "delta": "after minus before",
        },
        "checkpoint": str(checkpoint.resolve()),
        "checkpoint_sha256": _sha256(checkpoint),
        "device": str(device),
        "torch_version": torch.__version__,
        "python_version": sys.version,
        "temperature_kelvin": args.temperature,
        "window_ns_condition": args.window_ns,
        "preparation": preparation,
        "summaries": summaries,
        "peptide_rows": peptide_rows,
        "pymol": pymol_artifacts,
        "warnings": [
            "FlexPPI predicts residue properties on fixed coordinates; it does not generate structures.",
            "The before state is extracted from the bound assembly coordinates and is therefore an in-silico TCR deletion, not an experimentally solved apo pMHC conformer.",
            "Both complete assemblies exceed the pretraining max_residues=768 and are length extrapolations.",
        ],
    }
    _write_json(metadata_path, metadata)
    print(json.dumps({
        "report": str(report_path.resolve()),
        "peptide_csv": str(peptide_csv.resolve()),
        "metadata": str(metadata_path.resolve()),
        "pymol": pymol_artifacts,
    }, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()

