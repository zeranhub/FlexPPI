"""Run the unchanged FastProtFlex model on frozen complete structures.

Smoke runs exactly one median-size representative from each of the three sources.
Full runs require a separate explicit invocation. No fitting or label calibration.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import time
from datetime import datetime, timezone

import numpy as np
import pandas as pd
import psutil

HERE = Path(__file__).resolve().parent
COHORT = HERE.parent / "cohort"
UPSTREAM = Path(r"E:\Project\protein_flexibility_models\FastProtFlex")
RSCRIPT = Path(r"C:\Program Files\R\R-4.6.0\bin\Rscript.exe")
R_LIBRARIES = [r"E:/Project/protein_flexibility_models/R_library", r"E:/Project/protein_flexibility_models/r_libs"]
AA = dict(zip(
    "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
    "ARNDCQEGHILKMFPSTWYV"))


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2), encoding="utf-8")


def author_id(number, insertion):
    # Numeric PDB residue number plus original insertion code, if present.
    number = str(number).strip()
    try:
        number = str(int(float(number)))
    except ValueError:
        pass
    return number + str(insertion).strip()


def read_atoms(path):
    atoms = []
    for line in Path(path).read_text(encoding="ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = line[17:20].strip()
        element = line[76:78].strip() if len(line) >= 78 else ""
        if residue not in AA or element.upper() in {"H", "D"}:
            continue
        atoms.append({
            "atom": line[12:16].strip(), "residue": residue,
            "chain": line[21:22].strip() or "_",
            "author_residue_id": author_id(line[22:26], line[26:27]),
            "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54]),
        })
    return pd.DataFrame(atoms)


def identity_tuples(frame):
    return list(zip(frame.chain.astype(str), frame.author_residue_id.astype(str), frame.amino_acid.astype(str)))


def wait_for_memory(n_atoms, sample_id):
    # Measured single-process-tree peak: 2,944.76 MiB at 6,620 protein atoms.
    # The unchanged official all-at-once routine allocates dense distance/adjacency
    # arrays, so its peak grows approximately with the square of atom count.
    # This is a scheduling rule only: input, branch threshold, and model math stay exact.
    mib = 1024**2
    estimated_rss = (128.0 + (2944.76171875 - 128.0) * (n_atoms / 6620.0)**2) * mib
    required_available = int(1.15 * estimated_rss + 512 * mib)
    waiting_started = time.perf_counter()
    last_report = -30.0
    while psutil.virtual_memory().available < required_available:
        waited = time.perf_counter() - waiting_started
        if waited - last_report >= 30.0:
            status = {
                "status": "waiting_for_available_RAM", "sample_id": sample_id,
                "available_mib": psutil.virtual_memory().available / mib,
                "required_available_mib": required_available / mib,
                "estimated_R_process_tree_rss_mib": estimated_rss / mib,
                "wait_seconds": waited,
            }
            write_json(HERE / "waiting_for_memory.json", status)
            print(json.dumps(status), flush=True)
            last_report = waited
        time.sleep(5.0)
    return {
        "memory_wait_seconds": time.perf_counter() - waiting_started,
        "estimated_R_process_tree_rss_mib": estimated_rss / mib,
        "required_available_RAM_mib": required_available / mib,
        "available_RAM_before_launch_mib": psutil.virtual_memory().available / mib,
    }


def validate_input(record, truth):
    assert sha256(record.structure_path) == record.structure_sha256, record.sample_id
    atoms = read_atoms(record.structure_path)
    ca = atoms[atoms.atom == "CA"].copy()
    ca["amino_acid"] = ca.residue.map(AA)
    assert len(ca) == record.n_residues == len(truth), (record.sample_id, len(ca), len(truth))
    assert identity_tuples(ca) == identity_tuples(truth), (record.sample_id, "input CA identities differ")
    assert not ca.duplicated(["chain", "author_residue_id"]).any(), record.sample_id
    return atoms, ca


def launch_model(record, truth, run_mode):
    output_dir = HERE / "records" / record.sample_id
    output_dir.mkdir(parents=True, exist_ok=True)
    success_file = output_dir / "validation.json"
    if success_file.exists():
        old = json.loads(success_file.read_text(encoding="utf-8"))
        if old.get("status") == "complete" and old.get("input_sha256") == record.structure_sha256:
            assert sha256(output_dir / "prediction.npz") == old["prediction_sha256"]
            if (output_dir / "exception.json").exists():
                (output_dir / "exception.json").replace(output_dir / "exception_before_success.json")
            return old
    atoms, ca = validate_input(record, truth)
    for name in ("predict.r", "FUNCTION_GDV.r"):
        shutil.copy2(HERE / "upstream" / name, output_dir / name)
        assert sha256(output_dir / name) == sha256(UPSTREAM / name)
    shutil.copy2(record.structure_path, output_dir / "protein.pdb")
    shutil.copy2(HERE / "capture_full_precision.R", output_dir / "capture_full_precision.R")
    env = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                 "RCPP_PARALLEL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        env[name] = "1"
    env["R_LIBS_USER"] = ";".join(R_LIBRARIES)
    command = [str(RSCRIPT), "--vanilla", "capture_full_precision.R"]
    memory_schedule = wait_for_memory(len(atoms), record.sample_id)
    started = datetime.now(timezone.utc).isoformat()
    start = time.perf_counter()
    peak_rss = 0
    peak_threads = 0
    cpu_seconds = 0.0
    observed_pids = set()
    with open(output_dir / "stdout.log", "w", encoding="utf-8") as stdout, open(output_dir / "stderr.log", "w", encoding="utf-8") as stderr:
        process = subprocess.Popen(command, cwd=output_dir, env=env, stdout=stdout, stderr=stderr)
        monitor = psutil.Process(process.pid)
        affinity_cpu = psutil.cpu_count() - 1
        while process.poll() is None:
            try:
                processes = [monitor] + monitor.children(recursive=True)
                rss = 0
                thread_count = 0
                for item in processes:
                    try:
                        # Enforce one logical CPU for this subprocess tree as well as library thread caps.
                        item.cpu_affinity([affinity_cpu])
                        observed_pids.add(item.pid)
                        rss += item.memory_info().rss
                        thread_count += item.num_threads()
                        times = item.cpu_times()
                        cpu_seconds = max(cpu_seconds, times.user + times.system)
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                peak_rss = max(peak_rss, rss)
                peak_threads = max(peak_threads, thread_count)
            except psutil.NoSuchProcess:
                pass
            time.sleep(0.1)
        return_code = process.wait()
    elapsed = time.perf_counter() - start
    performance = {
        "sample_id": record.sample_id, "dataset": record.dataset, "mode": run_mode,
        "input_sha256": record.structure_sha256, "started_utc": started,
        "wall_seconds": elapsed, "observed_max_process_cpu_seconds": cpu_seconds,
        "peak_process_tree_rss_bytes": peak_rss, "peak_process_tree_rss_mib": peak_rss / 1024**2,
        "observed_process_ids": sorted(observed_pids), "cpu_affinity_logical_cpu": affinity_cpu,
        "thread_library_caps": 1, "peak_os_threads_including_runtime_threads": peak_threads,
        "return_code": return_code, "input_protein_heavy_atoms": len(atoms),
        "n_residues": int(record.n_residues), "n_chains": int(record.n_chains),
        "official_branch": "GDVallatonce" if len(atoms) <= 10000 else "GDVperpartes_new2",
        "official_atom_threshold": 10000,
        **memory_schedule,
    }
    write_json(output_dir / "execution_performance.json", performance)
    if return_code:
        performance["status"] = "failed_model_execution"
        write_json(success_file, performance)
        return performance
    raw = pd.read_csv(output_dir / "atom_predictions_full_precision.csv", keep_default_na=False)
    raw["chain"] = raw.chain.astype(str).replace("", "_")
    raw["author_residue_id"] = [author_id(a, b) for a, b in zip(raw["resno"], raw["insert"])]
    assert len(raw) == len(atoms), (record.sample_id, "all-atom row count mismatch")
    for key, upstream_key in (("atom", "elety"), ("residue", "resid"), ("chain", "chain"), ("author_residue_id", "author_residue_id")):
        assert atoms[key].astype(str).tolist() == raw[upstream_key].astype(str).tolist(), (record.sample_id, key)
    assert np.max(np.abs(atoms[["x", "y", "z"]].to_numpy() - raw[["x", "y", "z"]].to_numpy())) < 1e-9
    assert np.isfinite(raw.prediction.to_numpy(dtype=float)).all(), (record.sample_id, "nonfinite raw atoms")
    raw_ca = raw[raw.elety == "CA"].copy()
    raw_ca["amino_acid"] = raw_ca.resid.map(AA)
    assert identity_tuples(raw_ca) == identity_tuples(truth), (record.sample_id, "predicted CA identities differ")
    predictions = raw_ca.prediction.to_numpy(dtype=np.float64)
    np.savez_compressed(
        output_dir / "prediction.npz", prediction=predictions,
        target=truth.target_rmsf.to_numpy(dtype=np.float64),
        valid=truth.valid_rmsf.to_numpy(dtype=bool),
        chain=truth.chain.astype(str).to_numpy(dtype=str),
        author_residue_id=truth.author_residue_id.astype(str).to_numpy(dtype=str),
        amino_acid=truth.amino_acid.astype(str).to_numpy(dtype=str),
        residue_index=truth.residue_index.to_numpy(dtype=np.int64),
    )
    performance.update({
        "status": "complete", "ca_positions_verified": len(predictions),
        "all_atom_identities_and_coordinates_verified": True,
        "ca_identities_equal_frozen_truth_in_original_order": True,
        "all_ca_predictions_finite": True,
        "output_quantity": "normalized_relative_flexibility_CA_atom_score",
        "absolute_RMSF_Angstrom": False, "Hchi_supported": False,
        "reference_used_only_for_identity_and_evaluation": True,
        "prediction_sha256": sha256(output_dir / "prediction.npz"),
        "full_precision_atoms_sha256": sha256(output_dir / "atom_predictions_full_precision.csv"),
        "rds_exact_output_sha256": sha256(output_dir / "atom_predictions_exact.rds"),
        "official_pdb_output_sha256": sha256(output_dir / "proteinRMSF.pdb"),
    })
    write_json(success_file, performance)
    if (output_dir / "exception.json").exists():
        (output_dir / "exception.json").replace(output_dir / "exception_before_success.json")
    return performance


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--mode", choices=["smoke", "full"], required=True)
    args = parser.parse_args()
    assert RSCRIPT.is_file()
    (HERE / "upstream").mkdir(exist_ok=True)
    upstream_hashes = {}
    for name in ("predict.r", "FUNCTION_GDV.r", "LICENSE", "README.md"):
        destination = HERE / "upstream" / name
        if not destination.exists():
            shutil.copy2(UPSTREAM / name, destination)
        assert sha256(destination) == sha256(UPSTREAM / name)
        upstream_hashes[name] = sha256(destination)
    records = pd.read_csv(COHORT / "records.csv", keep_default_na=False)
    if args.mode == "full":
        sizes = pd.read_csv(HERE / "preflight_sizes.csv", keep_default_na=False)
        assert set(sizes.sample_id) == set(records.sample_id)
        records = records.merge(sizes[["sample_id", "protein_heavy_atoms"]], on="sample_id", validate="one_to_one")
        # Only scheduling changes: every frozen complete structure remains included.
        records = records.sort_values(["protein_heavy_atoms", "sample_id"], kind="stable").reset_index(drop=True)
        schedule = records[["sample_id", "dataset", "main_group_candidate", "n_chains", "n_residues", "protein_heavy_atoms"]].copy()
        schedule.insert(0, "schedule_order", range(1, len(schedule) + 1))
        schedule.to_csv(HERE / "full_schedule.csv", index=False)
    if args.mode == "smoke":
        selected = []
        for dataset in ("atlas", "dynarepo", "itrace"):
            part = records[records.dataset == dataset]
            if dataset == "dynarepo":
                part = part[part.n_chains > 1]
            selected.append(part.loc[(part.n_residues - part.n_residues.median()).abs().idxmin()])
        records = pd.DataFrame(selected)
    records.to_csv(HERE / (args.mode + "_inputs.csv"), index=False)
    truth_path = records.iloc[0].reference_table_path
    truth = pd.read_csv(truth_path, keep_default_na=False, low_memory=False,
                        dtype={"chain": str, "author_residue_id": str, "amino_acid": str})
    truth_groups = dict(tuple(truth.groupby("sample_id", sort=False)))
    write_json(HERE / "run_protocol.json", {
        "model": "FastProtFlex_official_GDV_linear_coefficients",
        "Rscript": str(RSCRIPT), "R_libraries": R_LIBRARIES,
        "upstream_directory": str(UPSTREAM), "upstream_sha256": upstream_hashes,
        "wrapper_sha256": sha256(__file__),
        "capture_wrapper_sha256": sha256(HERE / "capture_full_precision.R"),
        "complete_structure_input": True, "independent_chain_inference": False,
        "hydrogen_nonprotein_filter": "unchanged_official_bio3d_filters",
        "branch_threshold_atoms": 10000, "branch_rule": "official <=10000 all-at-once; >10000 official GDVperpartes_new2",
        "regression_coefficients_unchanged": True, "no_model_fitting": True,
        "normalization": "unchanged official atom-level scale(log(GDV+1)); score scale over all complete-structure atoms",
        "residue_output": "select CA atom prediction, no residue averaging or additional normalization",
        "raw_exact_output": "RDS preserves exact numeric arrays; CSV retained before two-decimal PDB B-factor rounding",
        "output_units": "normalized relative score, not Angstrom RMSF; no Hchi",
        "parallel_processes": 1, "logical_cpus_per_process_tree": 1, "library_threads": 1,
        "RAM_scheduling_only": "Before each structure wait for 1.15 * smoke-measured quadratic R peak plus 512 MiB reserve. No change to model branch rules.",
        "smoke_selection": "source-median residue count: ATLAS single-chain, DynaRepo multi-chain, ITrace TCR-pMHC",
        "input_cohort_manifest_sha256": sha256(COHORT / "records.csv"),
        "schedule_rule": "full: complete structures ordered by protein-heavy-atom count ascending, then sample_id; cached completed records are verified and reused",
        "full_schedule_sha256": sha256(HERE / "full_schedule.csv") if args.mode == "full" else None,
        "truth_sha256": sha256(truth_path),
    })
    results = []
    for position, record in enumerate(records.itertuples(index=False), 1):
        try:
            result = launch_model(record, truth_groups[record.sample_id], args.mode)
        except Exception as error:
            result = {"sample_id": record.sample_id, "dataset": record.dataset, "status": "failed_validation", "error": repr(error)}
            target = HERE / "records" / record.sample_id
            target.mkdir(parents=True, exist_ok=True)
            write_json(target / "exception.json", result)
        results.append(result)
        pd.DataFrame(results).to_csv(HERE / (args.mode + "_results.csv"), index=False)
        print(json.dumps({"progress": f"{position}/{len(records)}", **result}), flush=True)
    write_json(HERE / (args.mode + "_audit.json"), {
        "mode": args.mode, "records_requested": len(records),
        "records_complete": sum(item["status"] == "complete" for item in results),
        "status": "complete" if all(item["status"] == "complete" for item in results) else "some_failures",
        "results": results,
    })


if __name__ == "__main__":
    main()
