"""Prepared8058-atom resource probe. Default inspect only; never claims equivalence.

Actual --run requires the separate real three-source exact-output audit. Its
launch gate extrapolates measured resources, retaining15 percent plus256MiB,
and additionally requires system commit above predicted peak plus1024MiB,
and refuses launch unless current available RAM and competing-process checks
pass. All outputs are isolated; the original records/protocol are untouched.
"""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import run_local_lifetime_v2 as runner

SAMPLE = "7zh9_A"
N_ATOMS = 8058


def inspect_extreme():
    provenance, schedule, records = runner.inspect_sources()
    record = records[SAMPLE]
    assert runner.sha256(record["structure_path"]) == record["structure_sha256"]
    atoms = [line for line in Path(record["structure_path"]).read_text(encoding="ascii").splitlines()
             if line.startswith(("ATOM  ", "HETATM")) and line[17:20].strip() in runner.AA
             and line[76:78].strip().upper() not in {"H", "D"}]
    assert len(atoms) == N_ATOMS
    assert sum(line[12:16].strip() == "CA" for line in atoms) == int(record["n_residues"]) == 1024
    assert next(row for row in schedule if row["sample_id"] == SAMPLE)["protein_heavy_atoms"] == str(N_ATOMS)
    return provenance, record


def wait_extrapolated_resource_gate(calibration, maximum_seconds, psutil, output):
    """Extrapolate resource demand only, then require actual current RAM."""
    requirements = runner.resource_requirements(calibration,N_ATOMS)
    started = time.perf_counter()
    while True:
        available = psutil.virtual_memory().available / 1024**2
        commit = runner.windows_commit_status()
        jobs = runner.competing_jobs(psutil)
        elapsed = time.perf_counter()-started
        ready = available >= requirements["required_available_RAM_MiB"] and commit["available_system_commit_MiB"] is not None and commit["available_system_commit_MiB"] >= requirements["required_available_commit_MiB"] and not jobs
        gate = {"status":"ready" if ready else "waiting_for_RAM_commit_or_existing_R", "n_atoms":N_ATOMS,
                **requirements, **commit,
                "available_RAM_before_launch_MiB":available, "elapsed_wait_seconds":elapsed,
                "source_calibration_sha256":calibration["calibration_sha256"],
                "launch_gate_kind":"resource_N_squared_extrapolation_from_three_actual_peaks_only",
                "blocking_jobs":jobs}
        if ready:
            return gate
        if elapsed >= maximum_seconds:
            gate["status"] = "blocked_finite_memory_wait_timeout"
            return gate
        runner.write_json(output / "waiting_for_memory.json", gate)
        print(json.dumps(gate), flush=True)
        time.sleep(min(30, maximum_seconds-elapsed))


def load_one_truth(record, pd):
    found = []
    for chunk in pd.read_csv(record["reference_table_path"], keep_default_na=False, low_memory=False, chunksize=50000,
                             usecols=runner.TRUTH_COLUMNS, dtype={"chain":str, "author_residue_id":str, "amino_acid":str}):
        selected = chunk.loc[chunk.sample_id.astype(str) == SAMPLE]
        if len(selected):
            found.append(selected.copy())
    assert found
    truth = pd.concat(found, ignore_index=True)
    assert len(truth) == int(record["n_residues"]) == 1024
    return truth


def execute_R(directory, timeout_seconds, psutil):
    environment = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                 "RCPP_PARALLEL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        environment[name] = "1"
    environment["R_LIBS_USER"] = ";".join(runner.R_LIBRARIES)
    command = [str(runner.RSCRIPT), "--vanilla", "capture_full_precision.R"]
    peak_rss = peak_threads = 0
    peak_private_bytes = minimum_commit_MiB = minimum_physical_MiB = None
    known_private_samples = unknown_private_samples = known_commit_samples = unknown_commit_samples = 0
    cpu_times, owned = {}, {}
    affinity = psutil.cpu_count()-1
    timed_out = False
    started = time.perf_counter()
    with (directory / "stdout.log").open("w", encoding="utf-8") as stdout, (directory / "stderr.log").open("w", encoding="utf-8") as stderr:
        process = subprocess.Popen(command, cwd=directory, env=environment, stdout=stdout, stderr=stderr)
        monitor = psutil.Process(process.pid)
        while process.poll() is None:
            try:
                rss = threads = private_bytes = 0
                private_known = True
                for child in [monitor]+monitor.children(recursive=True):
                    try:
                        owned[child.pid] = child
                        child.cpu_affinity([affinity])
                        memory = child.memory_info()
                        rss += memory.rss
                        private = getattr(memory,"private",None)
                        if private is None:
                            private_known = False
                        else:
                            private_bytes += private
                        threads += child.num_threads()
                        usage = child.cpu_times()
                        cpu_times[child.pid] = max(cpu_times.get(child.pid, 0.0), usage.user+usage.system)
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        private_known = False
                peak_rss, peak_threads = max(peak_rss, rss), max(peak_threads, threads)
                if private_known:
                    known_private_samples += 1
                    peak_private_bytes = max(peak_private_bytes or 0,private_bytes)
                else:
                    unknown_private_samples += 1
                commit = runner.windows_commit_status()["available_system_commit_MiB"]
                if commit is None:
                    unknown_commit_samples += 1
                else:
                    known_commit_samples += 1
                    minimum_commit_MiB = commit if minimum_commit_MiB is None else min(minimum_commit_MiB,commit)
                available = psutil.virtual_memory().available/1024**2
                minimum_physical_MiB = available if minimum_physical_MiB is None else min(minimum_physical_MiB,available)
                if time.perf_counter()-started > timeout_seconds:
                    timed_out = True
                    for child in reversed(list(owned.values())):
                        try:
                            child.terminate()
                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            pass
                    process.terminate()
                    break
            except psutil.NoSuchProcess:
                pass
            time.sleep(.1)
        try:
            code = process.wait(timeout=10)
        except subprocess.TimeoutExpired:
            for child in owned.values():
                try:
                    child.kill()
                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    pass
            process.kill()
            code = process.wait(timeout=10)
    return {"command":command, "return_code":code, "timed_out":timed_out,
            "wall_seconds_R_capture_only":time.perf_counter()-started,
            "peak_R_process_tree_RSS_bytes":peak_rss, "peak_R_process_tree_RSS_MiB":peak_rss/1024**2,
            "peak_R_process_tree_private_bytes":peak_private_bytes,
            "private_bytes_counter_status":"unknown" if not known_private_samples else "available" if not unknown_private_samples else "partially_available",
            "known_private_samples":known_private_samples, "unknown_private_samples":unknown_private_samples,
            "minimum_system_available_commit_MiB":minimum_commit_MiB,
            "system_commit_sampling_status":"unknown" if not known_commit_samples else "available" if not unknown_commit_samples else "partially_available",
            "known_commit_samples":known_commit_samples, "unknown_commit_samples":unknown_commit_samples,
            "minimum_system_available_physical_MiB":minimum_physical_MiB,
            "observed_process_tree_CPU_seconds":sum(cpu_times.values()), "observed_process_ids":sorted(cpu_times),
            "cpu_affinity_logical_cpu":affinity, "library_thread_caps":1, "sampling_interval_seconds":.1,
            "peak_OS_threads_including_runtime_helpers":peak_threads}


def verify_outputs(directory, atoms, truth, np, pd):
    raw = pd.read_csv(directory / "atom_predictions_full_precision.csv", keep_default_na=False, float_precision="round_trip")
    raw["chain"] = raw.chain.astype(str).replace("", "_")
    raw["author_residue_id"] = [runner.author_id(a,b) for a,b in zip(raw.resno,raw["insert"])]
    assert len(raw) == len(atoms) == N_ATOMS
    for input_key, raw_key in (("atom","elety"),("residue","resid"),("chain","chain"),("author_residue_id","author_residue_id")):
        assert atoms[input_key].astype(str).tolist() == raw[raw_key].astype(str).tolist()
    difference = float(np.max(np.abs(atoms[["x","y","z"]].to_numpy()-raw[["x","y","z"]].to_numpy())))
    assert difference < 1e-9
    assert np.isfinite(raw.prediction.to_numpy(dtype=np.float64)).all()
    ca = raw.loc[raw.elety == "CA"].copy()
    ca["amino_acid"] = ca.resid.map(runner.AA)
    assert runner.identities(ca) == runner.identities(truth)
    assert not ca.duplicated(["chain", "author_residue_id"]).any()
    prediction = ca.prediction.to_numpy(dtype=np.float64)
    assert len(prediction) == 1024 and np.isfinite(prediction).all()
    np.savez_compressed(directory / "prediction.npz", prediction=prediction,
                        target=truth.target_rmsf.to_numpy(dtype=np.float64), valid=truth.valid_rmsf.to_numpy(dtype=bool),
                        chain=truth.chain.astype(str).to_numpy(dtype=str),
                        author_residue_id=truth.author_residue_id.astype(str).to_numpy(dtype=str),
                        amino_acid=truth.amino_acid.astype(str).to_numpy(dtype=str),
                        residue_index=truth.residue_index.to_numpy(dtype=np.int64))
    integrity = {"n_atoms":N_ATOMS, "n_CA":len(prediction), "all_atom_identities_equal_input_in_order":True,
                 "all_atom_coordinates_verified":True, "coordinate_maximum_absolute_difference":difference,
                 "coordinate_verification_tolerance":1e-9, "all_atom_predictions_finite":True,
                 "ca_identities_equal_frozen_truth_in_order":True, "all_ca_predictions_finite":True,
                 "no_duplicate_CA_identity":True, "output_quantity":"normalized_relative_flexibility_CA_atom_score",
                 "absolute_RMSF_Angstrom":False, "Hchi_supported":False,
                 "scope":"Resource and frozen-input/output-integrity verification only; no original8058-atom precision output exists, so no equivalence comparison was made."}
    runner.write_json(directory / "integrity.json", integrity)
    return integrity


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--validation-audit", type=Path)
    parser.add_argument("--run-tag", default="extreme_01")
    parser.add_argument("--memory-wait-seconds", type=float, default=600)
    parser.add_argument("--case-timeout-seconds", type=float, default=180)
    args = parser.parse_args()
    assert re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,63}", args.run_tag)
    assert args.memory_wait_seconds >= 0 and args.case_timeout_seconds > 0
    provenance, record = inspect_extreme()
    reference_available = (runner.HERE / "records" / SAMPLE / "atom_predictions_exact.rds").exists()
    calibration = runner.load_actual_calibration(args.validation_audit, provenance) if args.validation_audit else None
    if not args.run:
        print(json.dumps({"status":"static_inspection_only_no_R_no_original_changes", "sample_id":SAMPLE, "n_atoms":N_ATOMS,
                          "input_sha256":record["structure_sha256"], "actual_three_source_audit_supplied":calibration is not None,
                          "this_invocation_no_extreme_measurement":True, "original_full_precision_reference_available":reference_available,
                          "resource_policy_version":runner.RESOURCE_POLICY_VERSION,
                          "launch_requirement":"Actual3-source exact audit; no active Fast/R; current physical>=1.15peak+256MiB and system commit>=peak+1024MiB; unknown commit blocks launch."}, indent=2))
        return 0
    assert calibration is not None, "--run requires --validation-audit from an actual successful three-source exact-output validation"
    assert not reference_available, "An official8058-atom precision output now exists; stop and ask root to review the resource-only no-reference protocol"
    import numpy as np
    import pandas as pd
    import psutil
    from validate_candidate import session_versions
    assert runner.RSCRIPT.is_file()
    assert not runner.competing_jobs(psutil), "No active original Fast runner, lifetime-v2 runner, validation/probe, or R process is allowed"
    output = HERE / "resource_probes" / args.run_tag
    assert not output.exists(), "Choose a new run-tag; never overwrite a resource probe"
    output.mkdir(parents=True)
    directory = output / SAMPLE
    audit = {"status":"prepared_not_launched", "created_utc":datetime.now(timezone.utc).isoformat(), "sample_id":SAMPLE, "n_atoms":N_ATOMS,
             "input_sha256":record["structure_sha256"], "cohort_manifest_sha256":runner.sha256(runner.COHORT / "records.csv"),
             "schedule_sha256":runner.sha256(runner.HERE / "full_schedule.csv"),
             "candidate_function_sha256":provenance["candidate_function_sha256"], "official_function_sha256":provenance["official_function_sha256"],
             "official_predict_sha256":provenance["official_predict_sha256"], "capture_wrapper_sha256":provenance["capture_wrapper_sha256"],
             "numerical_equivalence_audit_path":calibration["validation_audit_path"], "numerical_equivalence_audit_sha256":calibration["validation_audit_sha256"],
             "base_three_source_resource_calibration":calibration,
             "probe_script_path":str(Path(__file__).resolve()), "probe_script_sha256":runner.sha256(__file__),
             "resume_entrypoint_path":str(Path(runner.__file__).resolve()), "resume_entrypoint_sha256":runner.sha256(runner.__file__),
             "resource_policy_version":runner.RESOURCE_POLICY_VERSION, "policy_change_reason":runner.POLICY_CHANGE_REASON,
             "original_full_precision_reference_available":False, "original_equivalence_compared":False,
             "claim_limit":"No original precision-equivalence claim for7zh9_A; three-source equivalence evidence is separate and mandatory."}
    exit_code = 2
    try:
        truth = load_one_truth(record, pd)
        atoms = runner.validate_input(record, truth, pd)
        assert len(atoms) == N_ATOMS
        gate = wait_extrapolated_resource_gate(calibration, args.memory_wait_seconds, psutil, output)
        audit["gate"] = gate
        if gate["status"] != "ready":
            audit["status"] = gate["status"]
            exit_code = 3
        else:
            directory.mkdir()
            for source, filename, expected_sha in ((runner.UPSTREAM / "predict.r", "predict.r", provenance["official_predict_sha256"]),
                                     (HERE / "candidate/FUNCTION_GDV.r", "FUNCTION_GDV.r", provenance["candidate_function_sha256"]),
                                     (runner.HERE / "capture_full_precision.R", "capture_full_precision.R", provenance["capture_wrapper_sha256"]),
                                     (Path(record["structure_path"]), "protein.pdb", record["structure_sha256"])):
                assert runner.sha256(source) == expected_sha
                shutil.copy2(source, directory / filename)
                assert runner.sha256(directory / filename) == expected_sha
            audit.update(execute_R(directory, args.case_timeout_seconds, psutil))
            assert audit["return_code"] == 0 and audit["timed_out"] is False
            assert audit["peak_R_process_tree_RSS_bytes"] > 0 and audit["observed_process_tree_CPU_seconds"] > 0
            actual_versions = session_versions(directory / "R_session_info.txt")
            reference_versions = session_versions(runner.HERE / "records/itrace_1qrn/R_session_info.txt")
            audit["R_and_package_versions_match_baseline"] = actual_versions == reference_versions
            audit["R_and_package_versions"] = actual_versions
            assert audit["R_and_package_versions_match_baseline"]
            audit["integrity"] = verify_outputs(directory, atoms, truth, np, pd)
            audit["combined_resource_envelope"] = {"version":runner.EXTENDED_GATE_VERSION,
                **runner.fit_resource_envelope(calibration["measured_points"] + [{"sample_id":SAMPLE, "n_atoms":N_ATOMS,
                                                                                "measured_peak_MiB":audit["peak_R_process_tree_RSS_MiB"]}])}
            audit["status"] = "extreme_resource_and_integrity_check_complete_no_reference_equivalence"
            exit_code = 0
    except Exception as error:
        audit["status"] = "failed_resource_or_integrity_do_not_adopt"
        audit["error"] = repr(error)
    audit["output_files"] = [{"relative_path":path.relative_to(output).as_posix(), "size_bytes":path.stat().st_size,
                              "sha256":runner.sha256(path)} for path in sorted(output.rglob("*")) if path.is_file()]
    audit["completed_utc"] = datetime.now(timezone.utc).isoformat()
    runner.write_json(output / "resource_extension_audit.json", audit)
    print(json.dumps({"status":audit["status"], "audit_path":str(output / "resource_extension_audit.json"),
                      "actual_peak_MiB":audit.get("peak_R_process_tree_RSS_MiB"), "original_equivalence_compared":False}), flush=True)
    return exit_code


if __name__ == "__main__":
    raise SystemExit(main())
