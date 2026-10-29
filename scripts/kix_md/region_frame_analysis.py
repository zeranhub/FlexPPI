"""Region-resolved, frame-resolved analysis of KIX MD-frame predictions.

Regions combine the pre-registered L12 loop (614-621), the whole KIX domain, and
data-derived interfaces (KIX residues in >=50 percent of sampled frames within
5 A heavy-atom contact of the cMyb or MLL peptide).  For every state, every
frame and every region the script reports absolute agreement and the paired
difference versus apo, so that a few good frames cannot be hidden by averaging.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from flexppi.evaluation import _correlation, _spearman  # noqa: E402

LOCAL_ROOT = ROOT / "artifacts/external/kix_md_20260914"
KIX_ROOT = ROOT / "artifacts/kix_md_20260914"
OUTPUT_DIRECTORY = LOCAL_ROOT / "frame_evaluation"
STATES = ("apo", "myb", "mll", "both")
CONTRASTS = (("myb", "apo"), ("mll", "apo"), ("both", "apo"))
PARTNERS = {"myb": ("A",), "mll": ("C",), "both": ("A", "C")}
TARGET_CHAIN = "B"
CONTACT_CUTOFF = 5.0
INTERFACE_FREQUENCY = 0.5
SAMPLING_STRIDE = 100
L12 = [str(number) for number in range(614, 622)]
DEAD_ZONE = {"rmsf": 0.1, "hchi": 0.1}


def compute_interfaces() -> pd.DataFrame:
    import MDAnalysis as mda

    cache = OUTPUT_DIRECTORY / "contact_frequency.csv"
    if cache.exists():
        return pd.read_csv(cache, dtype={"author_residue_id": str})
    rows = []
    for state in ("myb", "mll", "both"):
        topology = LOCAL_ROOT / state / "protein.pdb"
        segments = sorted((LOCAL_ROOT / state / "production").glob("segment_*.dcd"))
        universe = mda.Universe(str(topology), [str(segment) for segment in segments])
        kix = universe.select_atoms("segid {} and not element H".format(TARGET_CHAIN))
        kix_residues = np.asarray([atom.resid for atom in kix])
        counts = {partner: {} for partner in PARTNERS[state]}
        for partner in PARTNERS[state]:
            counts[partner] = {residue: 0 for residue in np.unique(kix_residues)}
        n_frames = 0
        for timestep in universe.trajectory[::SAMPLING_STRIDE]:
            n_frames += 1
            for partner in PARTNERS[state]:
                partner_atoms = universe.select_atoms(
                    "segid {} and not element H".format(partner))
                tree = cKDTree(partner_atoms.positions)
                close = tree.query_ball_point(kix.positions, CONTACT_CUTOFF)
                touched = set()
                for atom_index, neighbours in enumerate(close):
                    if neighbours:
                        touched.add(int(kix_residues[atom_index]))
                for residue in touched:
                    counts[partner][residue] += 1
        for partner in PARTNERS[state]:
            for residue, count in counts[partner].items():
                rows.append({"state": state, "partner": partner,
                             "author_residue_id": str(residue),
                             "frequency": count / max(n_frames, 1)})
        print("contacts", state, "frames", n_frames, flush=True)
    table = pd.DataFrame(rows)
    table.to_csv(cache, index=False)
    return table


def build_regions(contacts: pd.DataFrame) -> dict:
    def interface(state: str) -> set:
        subset = contacts[(contacts.state == state) & (contacts.frequency >= INTERFACE_FREQUENCY)]
        return set(subset.author_residue_id.astype(str))

    myb_interface = interface("myb")
    mll_interface = interface("mll")
    both_interface = interface("both")
    all_residues = sorted(set(contacts.author_residue_id.astype(str)), key=int)
    regions = {
        "whole_KIX": all_residues,
        "L12_614_621": L12,
        "myb_interface": sorted(myb_interface, key=int),
        "mll_interface": sorted(mll_interface, key=int),
        "shared_interface": sorted(myb_interface & mll_interface, key=int),
        "non_interface": sorted(set(all_residues) - (myb_interface | mll_interface), key=int),
        "both_interface": sorted(both_interface, key=int),
    }
    return regions


def per_frame_metric(prediction: pd.Series, target: pd.Series, dead_zone: float) -> dict:
    joined = pd.DataFrame({"pred": prediction, "md": target}).dropna()
    if len(joined) < 5:
        return {"n": len(joined), "spearman": None, "pearson": None, "mae": None,
                "sign_accuracy": None, "n_active": 0}
    prediction_values = joined["pred"].to_numpy(float)
    target_values = joined["md"].to_numpy(float)
    active = np.abs(target_values) > dead_zone
    return {"n": int(len(joined)), "spearman": _spearman(prediction_values, target_values),
            "pearson": _correlation(prediction_values, target_values),
            "mae": float(np.abs(prediction_values - target_values).mean()),
            "sign_accuracy": float((np.sign(prediction_values[active])
                                    == np.sign(target_values[active])).mean())
            if active.any() else None,
            "n_active": int(active.sum())}


def main() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    predictions = pd.read_csv(OUTPUT_DIRECTORY / "predictions_long.csv",
                              dtype={"author_residue_id": str})
    reference = pd.read_csv(
        KIX_ROOT / "pilot_results_20260916/comparison/matched_residues.csv",
        dtype={"author_residue_number": str})
    md = {}
    for (state, metric, residue), group in reference.groupby(
            ["state", "metric", "author_residue_number"], sort=False):
        md[(state, metric, residue)] = float(group["md_value"].iloc[0])

    contacts = compute_interfaces()
    regions = build_regions(contacts)
    (OUTPUT_DIRECTORY / "regions.json").write_text(json.dumps(regions, indent=2))

    fold_mean = (predictions.groupby(["state", "time_ns", "author_residue_id"])[
        ["pred_rmsf", "pred_hchi"]].mean().reset_index())
    absolute_rows = []
    delta_rows = []
    for region, residues in regions.items():
        residues = set(residues)
        for state in STATES:
            for metric, column in (("rmsf", "pred_rmsf"), ("hchi", "pred_hchi")):
                for time_ns, frame in fold_mean[fold_mean.state == state].groupby("time_ns"):
                    frame = frame.set_index("author_residue_id")
                    selected = frame[frame.index.isin(residues)]
                    target = pd.Series({residue: md.get((state, metric, residue), np.nan)
                                        for residue in selected.index})
                    values = per_frame_metric(selected[column], target, DEAD_ZONE[metric])
                    absolute_rows.append({"region": region, "state": state, "metric": metric,
                                          "time_ns": int(time_ns), **values})
        for left, right in CONTRASTS:
            for metric, column in (("rmsf", "pred_rmsf"), ("hchi", "pred_hchi")):
                left_frames = fold_mean[fold_mean.state == left]
                right_frames = fold_mean[fold_mean.state == right]
                for time_ns in sorted(left_frames.time_ns.unique()):
                    left_frame = left_frames[left_frames.time_ns == time_ns].set_index(
                        "author_residue_id")[column]
                    right_frame = right_frames[right_frames.time_ns == time_ns].set_index(
                        "author_residue_id")[column]
                    pred_delta = left_frame - right_frame
                    pred_delta = pred_delta[pred_delta.index.isin(residues)]
                    target = pd.Series({
                        residue: md.get((left, metric, residue), np.nan)
                        - md.get((right, metric, residue), np.nan)
                        for residue in pred_delta.index})
                    values = per_frame_metric(pred_delta, target, DEAD_ZONE[metric])
                    delta_rows.append({"region": region, "contrast": "{}-{}".format(left, right),
                                       "metric": metric, "time_ns": int(time_ns), **values})
    absolute = pd.DataFrame(absolute_rows)
    deltas = pd.DataFrame(delta_rows)
    absolute.to_csv(OUTPUT_DIRECTORY / "region_frame_absolute.csv", index=False)
    deltas.to_csv(OUTPUT_DIRECTORY / "region_frame_delta.csv", index=False)

    summary_rows = []
    for (region, contrast, metric), group in deltas.groupby(["region", "contrast", "metric"]):
        valid = group.dropna(subset=["spearman"])
        if valid.empty:
            summary_rows.append({
                "region": region, "contrast": contrast, "metric": metric, "n_frames": 0,
                "spearman_mean": np.nan, "spearman_sd": np.nan, "spearman_min": np.nan,
                "spearman_max": np.nan, "frac_positive": np.nan, "best_frame_ns": None,
                "sign_accuracy_mean": np.nan, "mae_mean": np.nan})
            continue
        summary_rows.append({
            "region": region, "contrast": contrast, "metric": metric, "n_frames": len(valid),
            "spearman_mean": valid["spearman"].mean(), "spearman_sd": valid["spearman"].std(ddof=1),
            "spearman_min": valid["spearman"].min(), "spearman_max": valid["spearman"].max(),
            "frac_positive": float((valid["spearman"] > 0).mean()),
            "best_frame_ns": int(valid.loc[valid["spearman"].idxmax(), "time_ns"]),
            "sign_accuracy_mean": valid["sign_accuracy"].mean(),
            "mae_mean": valid["mae"].mean()})
    summary = pd.DataFrame(summary_rows)
    summary.to_csv(OUTPUT_DIRECTORY / "region_frame_delta_summary.csv", index=False)

    best_frame_consistency = []
    for metric in ("rmsf", "hchi"):
        for left, right in CONTRASTS:
            contrast = "{}-{}".format(left, right)
            subset = deltas[(deltas.contrast == contrast) & (deltas.metric == metric)
                            & deltas.spearman.notna()]
            pivot = subset.pivot_table(index="region", columns="time_ns", values="spearman")
            if len(pivot) >= 2:
                best_frame_consistency.append({
                    "contrast": contrast, "metric": metric,
                    "region_pairwise_frame_correlation_mean": float(
                        pivot.T.corr().to_numpy()[np.triu_indices(len(pivot), 1)].mean()),
                    "n_regions": len(pivot)})
    consistency = pd.DataFrame(best_frame_consistency)
    consistency.to_csv(OUTPUT_DIRECTORY / "region_frame_consistency.csv", index=False)

    lines = ["# Region-resolved and frame-resolved KIX delta analysis", ""]
    lines.append("Regions: " + "; ".join(
        "{} (n={})".format(name, len(values)) for name, values in regions.items()))
    lines.append("")
    for metric in ("rmsf", "hchi"):
        lines.append("## {} delta versus apo".format(metric))
        lines.append("")
        lines.append("| contrast | region | n_frames | mean rho | SD | min | max | frac rho>0 | best frame (ns) |")
        lines.append("|---|---|---|---|---|---|---|---|---|")
        subset = summary[(summary.metric == metric)]
        for _, row in subset.iterrows():
            lines.append("| {} | {} | {} | {:+.3f} | {:.3f} | {:+.3f} | {:+.3f} | {:.2f} | {} |".format(
                row["contrast"], row["region"], row["n_frames"], row["spearman_mean"],
                row["spearman_sd"], row["spearman_min"], row["spearman_max"],
                row["frac_positive"], row["best_frame_ns"]))
        lines.append("")
    lines.append("Frame-profile correlation across regions (1.0 means a frame that is good for "
                 "one region is good for all):")
    for _, row in consistency.iterrows():
        lines.append("- {} {}: {:.3f} across {} regions".format(
            row["contrast"], row["metric"], row["region_pairwise_frame_correlation_mean"],
            row["n_regions"]))
    lines.append("")
    lines.append("Caveats: single MD repeat; input frames and reference labels come from the same "
                 "trajectories; best-frame numbers are descriptive and must not be selected as a "
                 "model claim; the frame SD together with the split-half noise file bounds what "
                 "frame choice could explain.")
    (OUTPUT_DIRECTORY / "SUMMARY_REGIONS.md").write_text("\n".join(lines) + "\n")
    print("wrote", OUTPUT_DIRECTORY / "SUMMARY_REGIONS.md", flush=True)


if __name__ == "__main__":
    main()
