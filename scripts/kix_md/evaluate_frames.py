"""Predict KIX MD frames with the frozen five-fold model and score against MD.

Inputs
------
frames/frames_manifest.json                     10 production frames per state
pilot_results_20260916/comparison/matched_residues.csv   MD per-residue RMSF/Hchi
frozen_predictions/residue_predictions.csv       the old single-input predictions
artifacts/performance_first_v2_atlas_hchi_v1/fold_*/epoch_040.pt   5 checkpoints

Outputs (frame_evaluation/)
---------------------------
predictions_long.csv      state x frame x fold x residue predictions
absolute_per_frame.csv    per-frame agreement with the MD reference
delta_metrics.csv         myb/mll/both minus apo, model versus MD
summary.json              frame-mean deltas, variance decomposition, noise floor
SUMMARY.md                compact report
"""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from flexppi.config import ModelConfig  # noqa: E402
from flexppi.data.dataset import collate_states  # noqa: E402
from flexppi.data.structures import parse_pdb  # noqa: E402
from flexppi.evaluation import _correlation, _spearman  # noqa: E402
from flexppi.inference import _structure_example  # noqa: E402
from flexppi.model.network import FlexPPI  # noqa: E402
from flexppi.trainer import move_to_device  # noqa: E402

LOCAL_ROOT = ROOT / "artifacts/external/kix_md_20260914"
KIX_ROOT = ROOT / "artifacts/kix_md_20260914"
CHECKPOINT_DIRECTORY = ROOT / "artifacts/performance_first_v2_atlas_hchi_v1"
OUTPUT_DIRECTORY = LOCAL_ROOT / "frame_evaluation"
STATES = ("apo", "myb", "mll", "both")
CONTRASTS = (("myb", "apo"), ("mll", "apo"), ("both", "apo"))
TARGET_CHAIN = "B"
L12_RESIDUES = tuple(str(number) for number in range(614, 622))
DEAD_ZONE = {"rmsf": 0.1, "hchi": 0.1}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_reference() -> pd.DataFrame:
    table = pd.read_csv(KIX_ROOT / "pilot_results_20260916/comparison/matched_residues.csv",
                        dtype={"author_residue_number": str})
    rows = []
    for (state, metric, residue), group in table.groupby(
            ["state", "metric", "author_residue_number"], sort=False):
        md_values = group["md_value"].to_numpy(float)
        if len(np.unique(md_values)) != 1:
            raise ValueError("MD value differs across folds for {} {} {}".format(
                state, metric, residue))
        rows.append({"state": state, "metric": metric, "author_residue_id": residue,
                     "amino_acid": group["amino_acid"].iloc[0], "md": md_values[0]})
    return pd.DataFrame(rows)


def load_frozen() -> pd.DataFrame:
    table = pd.read_csv(KIX_ROOT / "frozen_predictions/residue_predictions.csv",
                        dtype={"author_residue_id": str})
    table = table[table.chain == TARGET_CHAIN]
    rows = []
    for (state, residue), group in table.groupby(["state", "author_residue_id"], sort=False):
        rows.append({"state": state, "author_residue_id": residue,
                     "amino_acid": group["amino_acid"].iloc[0],
                     "frozen_rmsf": group["predicted_rmsf_angstrom"].mean(),
                     "frozen_hchi": group["predicted_hchi_nats"].mean()})
    return pd.DataFrame(rows)


def predict_frames(manifest: dict) -> pd.DataFrame:
    import torch

    torch.set_num_threads(6)
    device = torch.device("cpu")
    rows = []
    for state in STATES:
        frames = manifest["states"][state]["frames"]
        for frame in frames:
            pdb_path = ROOT / frame["pdb"]
            structure = parse_pdb(str(pdb_path))
            example = _structure_example(str(pdb_path), structure.unique_chains(), 300.0, 100.0, 0.0)
            batch = move_to_device(collate_states([example]), device)
            for fold in range(5):
                checkpoint_path = CHECKPOINT_DIRECTORY / "fold_{}".format(fold) / "epoch_040.pt"
                checkpoint = torch.load(str(checkpoint_path), map_location="cpu", weights_only=False)
                if checkpoint.get("epoch") != 40:
                    raise ValueError("{} is not epoch 40".format(checkpoint_path))
                config = ModelConfig(**checkpoint["config"]["model"])
                model = FlexPPI(config)
                model.load_state_dict(checkpoint["model"])
                model.to(device).eval()
                with torch.inference_mode():
                    output = model.forward_state(batch)
                mean = output["mean"][0].cpu().numpy()
                valid = (batch["target_mask"][0].numpy().astype(bool)
                         & batch["geometry_mask"][0].numpy().astype(bool)
                         & output["mask"][0].cpu().numpy().astype(bool))
                metric_index = {name: index for index, name in enumerate(config.metrics)}
                for index, chain in enumerate(structure.chain_ids):
                    if chain != TARGET_CHAIN:
                        continue
                    rows.append({
                        "state": state, "time_ns": frame["time_ns"], "fold": fold,
                        "chain": chain, "author_residue_id": structure.author_residue_ids[index],
                        "amino_acid": structure.sequence[index], "valid": bool(valid[index]),
                        "pred_rmsf": float(mean[index, metric_index["rmsf"]]),
                        "pred_hchi": float(mean[index, metric_index["chi_entropy"]]),
                    })
                del model
        print("predicted", state, len(frames), "frames x 5 folds", flush=True)
    return pd.DataFrame(rows)


def regression(prediction: np.ndarray, target: np.ndarray) -> dict:
    error = prediction - target
    return {"n": int(len(target)), "mae": float(np.abs(error).mean()),
            "bias": float(error.mean()), "pearson": _correlation(prediction, target),
            "spearman": _spearman(prediction, target)}


def main() -> None:
    OUTPUT_DIRECTORY.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((LOCAL_ROOT / "frames/frames_manifest.json").read_text())
    reference = load_reference()
    frozen = load_frozen()
    long_path = OUTPUT_DIRECTORY / "predictions_long.csv"
    if long_path.exists():
        predictions = pd.read_csv(long_path, dtype={"author_residue_id": str})
        print("reused", long_path, flush=True)
    else:
        predictions = predict_frames(manifest)
        predictions.to_csv(long_path, index=False)

    metric_column = {"rmsf": "pred_rmsf", "hchi": "pred_hchi"}
    frame_rows = []
    frame_mean_frames = []
    for state in STATES:
        state_reference = reference[reference.state == state]
        for metric in ("rmsf", "hchi"):
            target = state_reference[state_reference.metric == metric].set_index(
                "author_residue_id")["md"]
            per_frame = []
            for time_ns, group in predictions[predictions.state == state].groupby("time_ns"):
                fold_mean = group.groupby("author_residue_id")[metric_column[metric]].mean()
                joined = pd.DataFrame({"pred": fold_mean}).join(target, how="inner").dropna()
                values = regression(joined["pred"].to_numpy(float), joined["md"].to_numpy(float))
                frame_rows.append({"state": state, "metric": metric, "time_ns": time_ns, **values})
                per_frame.append(joined["pred"].rename(time_ns))
            frame_mean = pd.concat(per_frame, axis=1).mean(axis=1)
            frame_sd = pd.concat(per_frame, axis=1).std(axis=1, ddof=1)
            frame_mean_frames.append(pd.DataFrame({"state": state, "metric": metric,
                                                   "author_residue_id": frame_mean.index,
                                                   "pred_frame_mean": frame_mean.values,
                                                   "pred_frame_sd": frame_sd.values}))
    per_frame = pd.DataFrame(frame_rows)
    per_frame.to_csv(OUTPUT_DIRECTORY / "absolute_per_frame.csv", index=False)
    frame_means = pd.concat(frame_mean_frames, ignore_index=True)

    delta_rows = []
    for left, right in CONTRASTS:
        for metric in ("rmsf", "hchi"):
            md_left = reference[(reference.state == left) & (reference.metric == metric)].set_index(
                "author_residue_id")["md"]
            md_right = reference[(reference.state == right) & (reference.metric == metric)].set_index(
                "author_residue_id")["md"]
            pd_left = frame_means[(frame_means.state == left) & (frame_means.metric == metric)].set_index(
                "author_residue_id")["pred_frame_mean"]
            pd_right = frame_means[(frame_means.state == right) & (frame_means.metric == metric)].set_index(
                "author_residue_id")["pred_frame_mean"]
            joined = pd.DataFrame({"md_left": md_left, "md_right": md_right,
                                   "pred_left": pd_left, "pred_right": pd_right}).dropna()
            md_delta = joined["md_left"].to_numpy(float) - joined["md_right"].to_numpy(float)
            pred_delta = joined["pred_left"].to_numpy(float) - joined["pred_right"].to_numpy(float)
            active = np.abs(md_delta) > DEAD_ZONE[metric]
            delta_rows.append({"contrast": "{}-{}".format(left, right), "metric": metric,
                               **regression(pred_delta, md_delta),
                               "md_delta_abs_mean": float(np.abs(md_delta).mean()),
                               "pred_delta_abs_mean": float(np.abs(pred_delta).mean()),
                               "active_n": int(active.sum()),
                               "sign_accuracy": float((np.sign(pred_delta[active])
                                                       == np.sign(md_delta[active])).mean())
                               if active.any() else None})
    deltas = pd.DataFrame(delta_rows)
    deltas.to_csv(OUTPUT_DIRECTORY / "delta_metrics.csv", index=False)

    frozen_rows = []
    for state in STATES:
        for metric, column in (("rmsf", "frozen_rmsf"), ("hchi", "frozen_hchi")):
            target = reference[(reference.state == state) & (reference.metric == metric)].set_index(
                "author_residue_id")["md"]
            single = frozen[frozen.state == state].set_index("author_residue_id")[column]
            joined = pd.DataFrame({"pred": single}).join(target, how="inner").dropna()
            frozen_rows.append({"state": state, "metric": metric,
                                **regression(joined["pred"].to_numpy(float),
                                             joined["md"].to_numpy(float))})
    frozen_metrics = pd.DataFrame(frozen_rows)
    frozen_metrics.to_csv(OUTPUT_DIRECTORY / "frozen_single_input_metrics.csv", index=False)

    l12_rows = []
    l12_means = {}
    for state in STATES:
        for metric in ("rmsf", "hchi"):
            target = reference[(reference.state == state) & (reference.metric == metric)].set_index(
                "author_residue_id")["md"]
            model = frame_means[(frame_means.state == state) & (
                frame_means.metric == metric)].set_index("author_residue_id")["pred_frame_mean"]
            joined = pd.DataFrame({"md": target, "pred": model})
            joined = joined[joined.index.isin(L12_RESIDUES)]
            l12_means[(state, metric)] = (float(joined["md"].mean()), float(joined["pred"].mean()))
    for state in STATES:
        for metric in ("rmsf", "hchi"):
            md_mean, pred_mean = l12_means[(state, metric)]
            apo_md, apo_pred = l12_means[("apo", metric)]
            l12_rows.append({"state": state, "metric": metric,
                             "n": int(len(L12_RESIDUES)),
                             "md_mean": md_mean, "pred_mean": pred_mean,
                             "md_delta_vs_apo": None if state == "apo" else md_mean - apo_md,
                             "pred_delta_vs_apo": None if state == "apo" else pred_mean - apo_pred})
    l12 = pd.DataFrame(l12_rows)

    within_state_sd = frame_means.groupby(["state", "metric"])["pred_frame_sd"].mean()
    between_state_sd = (frame_means.pivot_table(index=["metric", "author_residue_id"],
                                                columns="state", values="pred_frame_mean")
                        .std(axis=1, ddof=1).groupby(level=0).mean())
    decomposition = {
        "mean_within_state_frame_sd": {metric: float(np.mean(
            [within_state_sd.loc[(state, metric)] for state in STATES]))
            for metric in ("rmsf", "hchi")},
        "mean_between_state_sd_of_frame_means": {metric: float(between_state_sd[metric])
                                                 for metric in ("rmsf", "hchi")}}

    split_rows = []
    for state in STATES:
        state_predictions = predictions[predictions.state == state]
        for metric, column in (("rmsf", "pred_rmsf"), ("hchi", "pred_hchi")):
            halves = {}
            for parity, name in ((0, "early"), (1, "late")):
                selected = state_predictions[(state_predictions.time_ns // 10) % 2 == parity]
                halves[name] = selected.groupby("author_residue_id")[column].mean()
            fake = (halves["early"] - halves["late"]).abs().mean()
            split_rows.append({"state": state, "metric": metric, "fake_delta_abs_mean": float(fake)})
    split = pd.DataFrame(split_rows)
    split.to_csv(OUTPUT_DIRECTORY / "split_half_model_noise.csv", index=False)

    summary = {
        "absolute_per_frame": {"{}|{}".format(*key): value for key, value in
                               per_frame.groupby(["state", "metric"])[
                                   ["spearman", "pearson", "mae"]].mean().round(4)
                               .to_dict(orient="index").items()},
        "absolute_per_frame_sd": {"{}|{}".format(*key): value for key, value in
                                  per_frame.groupby(["state", "metric"])[
                                      ["spearman", "mae"]].std(ddof=1).round(4)
                                  .to_dict(orient="index").items()},
        "frozen_single_input": frozen_metrics.round(4).to_dict(orient="records"),
        "deltas": deltas.round(4).to_dict(orient="records"),
        "variance_decomposition": decomposition,
        "split_half_model_noise": split.round(4).to_dict(orient="records"),
        "l12": l12.round(4).to_dict(orient="records"),
        "checkpoints": {str(fold): sha256(
            CHECKPOINT_DIRECTORY / "fold_{}".format(fold) / "epoch_040.pt") for fold in range(5)},
    }
    (OUTPUT_DIRECTORY / "summary.json").write_text(json.dumps(summary, indent=2, default=float))

    lines = ["# KIX MD-frame evaluation (10 frames per state, 5 folds)", ""]
    lines.append("## Absolute agreement with the MD reference (fold-mean per frame)")
    lines.append("")
    lines.append("| state | metric | frame-mean Spearman | frame SD | frame-mean MAE |")
    lines.append("|---|---|---|---|---|")
    for (state, metric), group in per_frame.groupby(["state", "metric"]):
        lines.append("| {} | {} | {:.3f} | {:.3f} | {:.3f} |".format(
            state, metric, group["spearman"].mean(), group["spearman"].std(ddof=1),
            group["mae"].mean()))
    lines.append("")
    lines.append("Old single-input (frozen) absolute Spearman: " + "; ".join(
        "{} {} {:.3f}".format(row["state"], row["metric"], row["spearman"])
        for row in summary["frozen_single_input"]))
    lines.append("")
    lines.append("## State deltas versus apo (frame-mean predictions)")
    lines.append("")
    lines.append("| contrast | metric | n | Spearman | Pearson | MAE | |MD delta| | |model delta| | sign acc |")
    lines.append("|---|---|---|---|---|---|---|---|---|")
    for row in summary["deltas"]:
        lines.append("| {} | {} | {} | {:.3f} | {:.3f} | {:.3f} | {:.3f} | {:.3f} | {} |".format(
            row["contrast"], row["metric"], row["n"], row["spearman"], row["pearson"], row["mae"],
            row["md_delta_abs_mean"], row["pred_delta_abs_mean"],
            "n/a" if row["sign_accuracy"] is None else "{:.3f}".format(row["sign_accuracy"])))
    lines.append("")
    lines.append("## Model frame-to-frame noise versus state differences")
    lines.append("")
    lines.append("Mean within-state frame SD (rmsf {:.3f} A, hchi {:.3f} nats); mean between-state SD "
                 "of frame means (rmsf {:.3f} A, hchi {:.3f} nats).".format(
                     decomposition["mean_within_state_frame_sd"]["rmsf"],
                     decomposition["mean_within_state_frame_sd"]["hchi"],
                     decomposition["mean_between_state_sd_of_frame_means"]["rmsf"],
                     decomposition["mean_between_state_sd_of_frame_means"]["hchi"]))
    lines.append("Split-half fake-delta magnitudes: " + "; ".join(
        "{} {} {:.3f}".format(row["state"], row["metric"], row["fake_delta_abs_mean"])
        for row in summary["split_half_model_noise"]))
    lines.append("")
    lines.append("Reference MD labels are single-repeat production RMSF/Hchi; model frame spread is "
                 "not MD uncertainty, and the input frames come from the same trajectories that "
                 "define the reference.")
    (OUTPUT_DIRECTORY / "SUMMARY.md").write_text("\n".join(lines) + "\n")
    print("wrote", OUTPUT_DIRECTORY, flush=True)


if __name__ == "__main__":
    main()
