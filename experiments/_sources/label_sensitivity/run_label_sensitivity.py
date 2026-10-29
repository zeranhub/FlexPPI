"""Frozen source-stratified temporal/alignment label sensitivity; CPU only.

Raw sources are read only. XTC headers can be scanned without writing source
offset files. Full trajectory reads occur only from a verified local staging
copy. Main analysis uses replica 1; three-replica means are supplementary.
"""
from __future__ import annotations

import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import gzip
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import random
import shutil
import struct
import sys
import tempfile
import time
import traceback
import types

HERE = Path(__file__).resolve().parent
CORRECTED = Path(r"E:\Project\FlexPPI\experiments\core_corrected_pipeline")
sys.path.insert(0, str(HERE / "runtime"))
sys.path.insert(0, str(HERE.parent / "single_target"))
sys.path.insert(0, str(CORRECTED / "training_payload" / "src"))
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")

# The public __init__ imports the neural model/Torch eagerly. This CPU-only
# analysis loads the original, unmodified data helper files through namespaces,
# so it does not need to install or instantiate the neural network.
for namespace_name, namespace_path in [("flexppi", CORRECTED / "training_payload" / "src" / "flexppi"),
                                       ("flexppi.data", CORRECTED / "training_payload" / "src" / "flexppi" / "data")]:
    namespace = types.ModuleType(namespace_name)
    namespace.__path__ = [str(namespace_path)]
    sys.modules[namespace_name] = namespace

import numpy as np
import pandas as pd
from completion_common import correlation, jsonl, sha256, statistics, write_json
from flexppi.constants import CHI_ATOMS, METRIC_NAMES, THREE_TO_ONE
from flexppi.data.labels import _dihedral_degrees, _entropy, _summarise_samples, save_label_file
from flexppi.data.structures import parse_pdb, target_mask

legacy_path = CORRECTED / "prepare_corrected_labels.py"
legacy_spec = importlib.util.spec_from_file_location("audited_corrected_labels", legacy_path)
legacy = importlib.util.module_from_spec(legacy_spec)
legacy_spec.loader.exec_module(legacy)


def scan_xtc_headers(path):
    """Read complete XTC frame metadata while skipping compressed coordinates."""
    times, offsets, steps = [], [], []
    byte_count = Path(path).stat().st_size
    atom_count = None
    with Path(path).open("rb") as handle:
        while handle.tell() < byte_count:
            offset = handle.tell()
            raw = handle.read(52)
            if len(raw) != 52:
                raise ValueError("Truncated XTC header")
            header = struct.unpack(">iiif9f", raw)
            magic, atoms, step, timestamp = header[:4]
            if magic != 1995 or atoms <= 0 or (atom_count is not None and atoms != atom_count):
                raise ValueError("Invalid/changing XTC atom count")
            atom_count = atoms
            coordinate_atoms = struct.unpack(">i", handle.read(4))[0]
            if coordinate_atoms != atoms:
                raise ValueError("XTC coordinate/header atom mismatch")
            if atoms <= 9:
                handle.seek(12 * atoms, 1)
            else:
                compressed_header = handle.read(36)
                if len(compressed_header) != 36:
                    raise ValueError("Truncated XTC compression header")
                payload = struct.unpack(">i", compressed_header[-4:])[0]
                if payload <= 0:
                    raise ValueError("Invalid XTC compressed byte count")
                handle.seek(4 * ((payload + 3) // 4), 1)
            if handle.tell() > byte_count:
                raise ValueError("Truncated XTC frame")
            offsets.append(offset)
            times.append(float(timestamp))
            steps.append(step)
    values = np.asarray(times)
    if len(values) < 2 or not np.isfinite(values).all():
        raise ValueError("Non-finite or insufficient trajectory timestamps")
    differences = np.diff(values)
    return dict(times_ps=times, frame_offsets=offsets, atom_count=atom_count, steps=steps,
                frames=len(times), first_ps=float(values[0]), last_ps=float(values[-1]),
                median_step_ps=float(np.median(differences[differences > 0])), min_step_ps=float(differences.min()),
                max_step_ps=float(differences.max()), duplicate_time_indices=np.flatnonzero(differences == 0).tolist(),
                backwards_time_indices=np.flatnonzero(differences < 0).tolist())


def freeze(args):
    output = args.output
    output.mkdir(parents=True, exist_ok=True)
    plan_path = output / "sensitivity_plan.json"
    if plan_path.exists():
        plan = json.loads(plan_path.read_text(encoding="utf-8"))
        assert plan["seed"] == args.seed and plan["per_source"] == args.per_source
        assert plan["analysis_code_sha256"] == sha256(__file__)
        print("Using existing frozen plan", flush=True)
        return plan
    all_records = jsonl(CORRECTED / "training_payload" / "splits" / "all_ready.jsonl")
    rng = random.Random(args.seed)
    chosen = []
    for source in ["atlas", "dynarepo", "itrace"]:
        candidates = sorted([row for row in all_records if row["dataset"] == source], key=lambda row: row["sample_id"])
        chosen += sorted(rng.sample(candidates, args.per_source), key=lambda row: row["sample_id"])
    # Persist identities before raw inspection, without looking at labels/predictions.
    write_json(output / "selected_records_before_raw_inspection.json", dict(seed=args.seed, per_source=args.per_source,
        records=[dict(sample_id=row["sample_id"], dataset=row["dataset"], fold=row["metadata"]["fold"]) for row in chosen]))
    old_manifest = {row["sample_id"]: row for row in jsonl(CORRECTED / "neural_exports" / "evaluation" / "prediction_manifest.jsonl") if row["variant"] == "full"}
    rows, stamp_rows = [], []
    for index, record in enumerate(chosen):
        sample = record["sample_id"]
        cached_path = output / "input_preflight" / (sample + ".json")
        if cached_path.exists():
            row = json.loads(cached_path.read_text(encoding="utf-8"))
            assert row["sample_id"] == sample and row["dataset"] == record["dataset"]
            assert all(Path(item["trajectory"]).stat().st_size == item["trajectory_bytes"] for item in row["trajectories"])
            for item in row["trajectories"]:
                for frame, timestamp in enumerate(item["times_ps"]):
                    stamp_rows.append(dict(sample_id=sample, dataset=row["dataset"], replica_index=item["replica_index"], frame=frame, time_ps=timestamp))
            rows.append(row)
            print("cached preflight", index + 1, len(chosen), sample, flush=True)
            continue
        audit_path = CORRECTED / "dataset" / "audit" / (sample + ".json")
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        sources = audit["trajectory_sources"]
        if record["dataset"] == "atlas":
            folder = Path("E:/data/protein_md_datasets/ATLAS_extracted") / sample
            sources = [dict(topology=str(folder / (sample + ".pdb")), trajectory=str(folder / (sample + "_R" + str(rep) + ".xtc"))) for rep in [1, 2, 3]]
        structure = CORRECTED / "dataset" / "structures" / (sample + ".pdb")
        label = CORRECTED / "dataset" / "labels" / (sample + ".npz")
        prediction = CORRECTED / "neural_exports" / "evaluation" / "full" / (sample + ".npz")
        row = dict(sample_id=sample, dataset=record["dataset"], fold=int(record["metadata"]["fold"]),
                   target_chains=record["target_chains"], structure=str(structure), old_label=str(label), prediction=str(prediction),
                   structure_sha256=sha256(structure), old_label_sha256=sha256(label), prediction_sha256=sha256(prediction),
                   old_prediction_manifest=old_manifest[sample], trajectories=[])
        assert row["structure_sha256"] == old_manifest[sample]["input_sha256"]
        assert row["old_label_sha256"] == old_manifest[sample]["label_sha256"]
        for number, source in enumerate(sources):
            top, trajectory = Path(source["topology"]), Path(source["trajectory"])
            metadata = scan_xtc_headers(trajectory)
            item = dict(replica_index=number + 1, topology=str(top), trajectory=str(trajectory),
                        topology_sha256=sha256(top), trajectory_sha256=sha256(trajectory), trajectory_bytes=trajectory.stat().st_size,
                        **metadata)
            if "trajectory_sha256" in source:
                assert item["trajectory_sha256"] == source["trajectory_sha256"]
                assert item["topology_sha256"] == source["topology_sha256"]
            for frame, timestamp in enumerate(metadata["times_ps"]):
                stamp_rows.append(dict(sample_id=sample, dataset=record["dataset"], replica_index=number + 1, frame=frame, time_ps=timestamp))
            row["trajectories"].append(item)
        rows.append(row)
        write_json(output / "input_preflight" / (sample + ".json"), row)
        print("preflight", index + 1, len(chosen), sample, flush=True)
    # Define eligibility using timestamps only, before computing alternative labels.
    selected_count=len(rows)
    excluded=[]
    eligible=[]
    for row in rows:
        reasons=[]
        if any(item['last_ps']<10000.0 for item in row['trajectories']):reasons.append('less than 10 ns available; no useful long-window comparison')
        if any(item['backwards_time_indices'] for item in row['trajectories']):reasons.append('backwards trajectory timestamp; cannot assign a unique shared time axis')
        if reasons:excluded.append(dict(sample_id=row['sample_id'],dataset=row['dataset'],reasons=reasons))
        else:eligible.append(row)
    write_json(output/'timestamp_exclusions.json',dict(selected=selected_count,eligible=len(eligible),excluded=excluded,
        rule='At least 10 ns in each saved trajectory and monotone timestamps; set before alternative label calculation; no replacement.'))
    original_rows=rows
    rows=eligible
    all_trajectories = [item for row in rows for item in row["trajectories"]]
    write_json(output / "all_selected_inputs_before_window_definition.json", dict(seed=args.seed, records=original_rows))
    assert all(item["first_ps"] == 0 for item in all_trajectories)
    # A short input changes the window for everybody; it is never replaced based on score.
    common_end = min(100000.0, min(item["last_ps"] for item in all_trajectories))
    # ITrace stores 200 ps frames, so 100 ps cannot be a shared observed grid.
    step_ps=200.0
    end_ps = np.floor(common_end / step_ps) * step_ps
    grid = np.arange(0.0, end_ps + 0.1, step_ps)
    if len(grid) < 2:
        raise ValueError("No useful common window")
    for item in all_trajectories:
        times = np.asarray(item["times_ps"])
        nearest = np.searchsorted(times, grid)
        nearest = np.clip(nearest, 0, len(times) - 1)
        previous = np.maximum(nearest - 1, 0)
        choose_previous = abs(times[previous] - grid) < abs(times[nearest] - grid)
        nearest[choose_previous] = previous[choose_previous]
        difference = abs(times[nearest] - grid)
        if difference.max() > 0.05 or len(np.unique(nearest)) != len(grid):
            raise ValueError("An eligible input cannot support the shared 200 ps grid; no interpolation allowed")
        item["selected_indices"] = nearest.tolist()
        item["max_grid_time_error_ps"] = float(difference.max())
    plan = dict(protocol="common_window_chain_Kabsch_replica1_v1", seed=args.seed, per_source=args.per_source,
                selected_records=selected_count,eligible_records=len(rows),timestamp_exclusions=excluded,
                sampling="sorted source records; Python random.Random seed; uniform without replacement; no score/label based selection",
                common_start_ps=0.0, common_end_ps=float(end_ps), common_step_ps=step_ps, common_frames=len(grid),
                main="replica1 only for every source", supplementary="available3 mean for ATLAS/ITrace; one trajectory for DynaRepo",
                source_hashes=dict(legacy_label_code=sha256(legacy_path), freeze=sha256(CORRECTED / "training_payload" / "training_freeze.json")),
                analysis_code_sha256=sha256(__file__), records=rows, output=str(output.resolve()),
                rmsf="all CA atoms per chain; each frame Kabsch to first selected frame; fluctuation around aligned mean; Å",
                hchi="all defined chi atoms required; fixed [-180,-60),[-60,60),[60,180] bins; joint 3**q codes; natural-log entropy; nats",
                exclusion="No replacement; timestamp-only >=10 ns eligibility before labels; mapping failures reported; identical paired valid intersection; no outlier trimming",
                duplicate_time_rule="searchsorted left chooses first frame at duplicate time; no interpolation, backwards timestamps excluded",
                boundary="Equal estimator/time/sampling/primary replica count, not equal MD forcefields, temperatures or independent families")
    pd.DataFrame(stamp_rows).to_csv(output / "all_selected_trajectory_timestamps.csv.gz", index=False, compression="gzip")
    write_json(plan_path, plan)
    time_summary = []
    for row in rows:
        for item in row["trajectories"]:
            time_summary.append(dict(sample_id=row["sample_id"], dataset=row["dataset"], replica_index=item["replica_index"],
                                     frames=item["frames"], first_ps=item["first_ps"], last_ps=item["last_ps"], step_ps=item["median_step_ps"],
                                     selected_frames=len(item["selected_indices"]), max_grid_time_error_ps=item["max_grid_time_error_ps"]))
    pd.DataFrame(time_summary).to_csv(output / "trajectory_timestamp_audit.csv", index=False)
    print(json.dumps(dict(records=len(rows), trajectories=len(all_trajectories), common_end_ns=end_ps / 1000, frames=len(grid)), indent=2), flush=True)
    return plan


def extract_uniform(reference, top, trajectory, indices, times):
    import MDAnalysis as mda
    replica = parse_pdb(str(top))
    universe = mda.Universe(str(top), str(trajectory))
    try:
        assert len(universe.trajectory) == len(times)
        residues = [residue for residue in universe.residues if residue.resname in THREE_TO_ONE and "CA" in residue.atoms.names]
        keys = [(str(residue.atoms.chainIDs[0]).strip() or "_", str(residue.resid) + str(residue.icode).strip(), str(residue.resname)) for residue in residues]
        assert keys == list(zip(replica.chain_ids, replica.author_residue_ids, replica.residue_names))
        ca, quads, chis = [], [], {}
        for residue_index, residue in enumerate(residues):
            lookup = {atom.name: atom.index for atom in residue.atoms}
            ca.append(lookup["CA"])
            definitions = CHI_ATOMS.get(residue.resname, [])
            if definitions and all(all(atom in lookup for atom in quad) for quad in definitions):
                chis[residue_index] = list(range(len(quads), len(quads) + len(definitions)))
                quads += [[lookup[atom] for atom in quad] for quad in definitions]
        quad_array = np.asarray(quads, dtype=int)
        coordinates, states, bad_bonds = [], [], 0
        for frame_index in indices:
            frame = universe.trajectory[frame_index]
            assert abs(frame.time - times[frame_index]) <= 0.05
            coordinates.append(frame.positions[ca].copy())
            if len(quads):
                points = frame.positions[quad_array]
                # Diagnose raw geometry without silently filtering frames or entropy labels.
                bad_bonds += int((np.linalg.norm(np.diff(points, axis=1), axis=-1) > 5.0).sum())
                angles = _dihedral_degrees(points)
                states.append(np.clip(np.floor((angles + 180.0) / 120.0).astype(np.int8), 0, 2))
        sample = np.full((replica.num_residues, len(METRIC_NAMES)), np.nan)
        sample[:, 0] = legacy.chain_fit(np.asarray(coordinates), replica.chain_ids)
        states = np.asarray(states)
        for residue_index, columns in chis.items():
            codes = sum(states[:, column].astype(np.int32) * 3 ** exponent for exponent, column in enumerate(columns))
            sample[residue_index, 1] = _entropy(codes)
        mapping, mapping_stats = legacy.replica_map(reference, replica)
        mapped = np.full((reference.num_residues, len(METRIC_NAMES)), np.nan)
        for target_index, source_index in mapping:
            mapped[target_index] = sample[source_index]
        return mapped, dict(mapping=mapping_stats, frames=len(indices), original_frames=len(times),
                            chi_residues=len(chis), chi_bond_instances_gt5A=bad_bonds,
                            residue_coverage=len(mapping) / reference.num_residues)
    finally:
        universe.trajectory.close()


def process_record(arguments):
    row, output_str, plan_sha = arguments
    output = Path(output_str)
    audit_path = output / "record_audit" / (row["sample_id"] + ".json")
    if audit_path.exists():
        existing = json.loads(audit_path.read_text(encoding="utf-8"))
        if existing.get("status") == "complete" and existing.get("plan_sha256") == plan_sha:
            return existing
    started = time.monotonic()
    report = dict(sample_id=row["sample_id"], dataset=row["dataset"], fold=row["fold"], plan_sha256=plan_sha)
    try:
        assert sha256(row["structure"]) == row["structure_sha256"]
        assert sha256(row["old_label"]) == row["old_label_sha256"]
        assert sha256(row["prediction"]) == row["prediction_sha256"]
        reference = parse_pdb(row["structure"])
        staging = output / "_staging"
        staging.mkdir(exist_ok=True)
        samples, audits = [], []
        for item in row["trajectories"]:
            with tempfile.TemporaryDirectory(prefix="selected_trajectory_", dir=staging) as folder:
                folder = Path(folder)
                assert folder.resolve().is_relative_to(staging.resolve())
                staged = folder / "trajectory.xtc"
                # Copy/hash once; MDAnalysis may write offset caches only to this temporary copy.
                digest = hashlib.sha256()
                with Path(item["trajectory"]).open("rb") as source, staged.open("wb") as destination:
                    for block in iter(lambda: source.read(4 * 1024 * 1024), b""):
                        destination.write(block)
                        digest.update(block)
                assert digest.hexdigest() == item["trajectory_sha256"]
                assert sha256(item["topology"]) == item["topology_sha256"]
                top = folder / "canonical.pdb"
                legacy.canonical_pdb(item["topology"], top)
                sample, audit = extract_uniform(reference, top, staged, item["selected_indices"], item["times_ps"])
                samples.append(sample)
                audits.append(audit | dict(replica_index=item["replica_index"]))
        label_paths = {}
        for variant, selected in [("uniform_rep1", samples[:1]), ("uniform_available_replicas", samples)]:
            values, variances, mask = _summarise_samples(selected, reference.num_residues)
            mask &= reference.frame_mask[:, None] & target_mask(reference, row["target_chains"])[:, None]
            path = output / "labels" / variant / (row["sample_id"] + ".npz")
            save_label_file(str(path), values, variances, mask)
            label_paths[variant] = dict(path=str(path), sha256=sha256(path), valid_counts=mask[:, :2].sum(axis=0).tolist())
        report.update(status="complete", replica_audits=audits, labels=label_paths, replicas=len(samples), elapsed_seconds=time.monotonic() - started)
    except Exception as error:
        report.update(status="failed", error=repr(error), traceback=traceback.format_exc(), elapsed_seconds=time.monotonic() - started)
    write_json(audit_path, report)
    return report


def compute(args, plan):
    assert sha256(__file__) == plan["analysis_code_sha256"]
    plan_sha = sha256(args.output / "sensitivity_plan.json")
    work = [(row, str(args.output), plan_sha) for row in plan["records"]]
    results = []
    with ProcessPoolExecutor(max_workers=args.workers) as pool:
        for index, result in enumerate(pool.map(process_record, work)):
            results.append(result)
            print("compute", index + 1, len(work), result["sample_id"], result["status"], flush=True)
    write_json(args.output / "compute_audit.json", dict(status="complete" if all(row["status"] == "complete" for row in results) else "partial",
               records=len(results), completed=sum(row["status"] == "complete" for row in results),
               failures=[row for row in results if row["status"] != "complete"], plan_sha256=plan_sha))


def analyze(args, plan):
    metrics, label_changes, masks = [], [], []
    for row in plan["records"]:
        report = json.loads((args.output / "record_audit" / (row["sample_id"] + ".json")).read_text(encoding="utf-8"))
        if report["status"] != "complete":
            continue
        with np.load(row["old_label"], allow_pickle=False) as old:
            old_values, old_mask = old["values"].astype(float), old["mask"].astype(bool)
        with np.load(row["prediction"], allow_pickle=False) as cache:
            prediction, original_valid = cache["mean"].astype(float), cache["valid"].astype(bool)
        identity = dict(sample_id=row["sample_id"], dataset=row["dataset"], fold=row["fold"])
        for variant, source in report["labels"].items():
            with np.load(source["path"], allow_pickle=False) as new:
                new_values, new_mask = new["values"].astype(float), new["mask"].astype(bool)
            for metric_index, metric in enumerate(METRIC_NAMES[:2]):
                common = original_valid[:, metric_index] & old_mask[:, metric_index] & new_mask[:, metric_index]
                common &= np.isfinite(old_values[:, metric_index]) & np.isfinite(new_values[:, metric_index]) & np.isfinite(prediction[:, metric_index])
                masks.append(identity | dict(variant=variant, metric=metric, old_valid=int(original_valid[:, metric_index].sum()),
                          new_valid=int(new_mask[:, metric_index].sum()), common=int(common.sum()),
                          old_only=int((original_valid[:, metric_index] & ~new_mask[:, metric_index]).sum()),
                          new_only=int((new_mask[:, metric_index] & ~original_valid[:, metric_index]).sum())))
                if not common.any():
                    continue
                y_old, y_new, p = old_values[common, metric_index], new_values[common, metric_index], prediction[common, metric_index]
                delta = y_new - y_old
                label_changes.append(identity | dict(variant=variant, metric=metric, count=len(delta),
                   label_mae_shift=float(abs(delta).mean()), label_bias_shift=float(delta.mean()),
                   label_rmse_shift=float(np.sqrt(np.mean(delta ** 2))), label_old_new_spearman=correlation(y_old, y_new, True)))
                before, after = statistics(y_old, p), statistics(y_new, p)
                metric_row = identity | dict(variant=variant, metric=metric, count=len(delta))
                for name in ["mae", "rmse", "spearman", "pearson", "bias"]:
                    metric_row["old_" + name] = before[name]
                    metric_row["new_" + name] = after[name]
                    metric_row["new_minus_old_" + name] = after[name] - before[name]
                metrics.append(metric_row)
    for name, values in [("label_changes_per_record", label_changes), ("prediction_effect_per_record", metrics), ("valid_mask_audit", masks)]:
        pd.DataFrame(values).to_csv(args.output / (name + ".csv"), index=False)
    frame = pd.DataFrame(metrics)
    numbers = [column for column in frame if column.startswith(("old_", "new_"))]
    summary = frame.groupby(["variant", "dataset", "metric"])[numbers].agg(["mean", "std"])
    summary.columns = [name + "_" + scope for name, scope in summary.columns]
    summary.reset_index().to_csv(args.output / "prediction_effect_record_mean_sample_sd.csv", index=False)
    pd.DataFrame(label_changes).groupby(["variant", "dataset", "metric"]).mean(numeric_only=True).to_csv(args.output / "label_change_macro.csv")
    write_json(args.output / "analysis_audit.json", dict(status="complete", selected=len(plan["records"]),
        completed_samples=int(frame.sample_id.nunique()), effective_samples_by_variant_source_metric=frame.groupby(["variant", "dataset", "metric"]).size().reset_index(name="n").to_dict(orient="records"),
        plan_sha256=sha256(args.output / "sensitivity_plan.json"), code_sha256=sha256(__file__),
        interpretation="Frozen OOF means re-scored under alternate labels on identical paired valid intersections; no retraining and no claim of uniform simulation conditions",
        variability="Sample SD across selected records within each source; not five-fold SD or SE"))
    print(summary.to_string(), flush=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--phase", choices=["freeze", "compute", "analyze", "all"], default="all")
    parser.add_argument("--output", type=Path, default=HERE)
    parser.add_argument("--per-source", type=int, default=30)
    parser.add_argument("--seed", type=int, default=20261010)
    parser.add_argument("--workers", type=int, default=4)
    args = parser.parse_args()
    args.output = args.output.resolve()
    plan = freeze(args)
    if args.phase in ["compute", "all"]:
        compute(args, plan)
    if args.phase in ["analyze", "all"]:
        analyze(args, plan)


if __name__ == "__main__":
    main()
