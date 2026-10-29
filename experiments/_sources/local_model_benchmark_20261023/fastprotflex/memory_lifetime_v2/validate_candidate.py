"""Explicit, isolated, serial three-source validation. Default is inspect only.

No R process is launched without --run. An active original FastProtFlex wrapper
or any existing R process blocks launch; the original RAM gate is retained.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time

HERE = Path(__file__).resolve().parent
BASELINE = HERE.parent
RSCRIPT = Path("C:/Program Files/R/R-4.6.0/bin/Rscript.exe")
R_LIBRARIES = ["E:/Project/protein_flexibility_models/R_library", "E:/Project/protein_flexibility_models/r_libs"]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for data in iter(lambda: handle.read(1024**2), b""):
            digest.update(data)
    return digest.hexdigest()


def write_json(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2), encoding="utf-8")
    temporary.replace(path)


def static_checks():
    from prepare_static import numeric_program_lines
    provenance = json.loads((HERE / "lifecycle_provenance.json").read_text(encoding="utf-8"))
    original = HERE / "original/FUNCTION_GDV.r"
    candidate = HERE / "candidate/FUNCTION_GDV.r"
    assert sha256(original) == provenance["official_function_sha256"]
    assert sha256(Path(provenance["official_function_path"])) == provenance["official_function_sha256"]
    assert sha256(candidate) == provenance["candidate_function_sha256"]
    # read_bytes preserves the original CRLF; the comparison only removes rm/gc.
    assert numeric_program_lines(original.read_bytes().decode("utf-8")) == numeric_program_lines(candidate.read_bytes().decode("utf-8"))
    assert sha256(HERE / "wrapper/predict.r") == provenance["official_predict_sha256"]
    assert sha256(HERE / "wrapper/capture_full_precision.R") == provenance["capture_wrapper_sha256"]
    assert sha256(BASELINE / "capture_full_precision.R") == provenance["capture_wrapper_sha256"]
    for case in provenance["cases"]:
        assert sha256(Path(case["baseline_RDS"])) == case["baseline_RDS_sha256"]
        assert sha256(Path(case["baseline_validation"])) == case["baseline_validation_sha256"]
        pdb = BASELINE / "records" / case["sample_id"] / "protein.pdb"
        assert sha256(pdb) == case["input_sha256"]
    return provenance


def competing_jobs(psutil):
    blocked = []
    for process in psutil.process_iter(["pid", "name", "cmdline"]):
        try:
            name = (process.info["name"] or "").lower()
            command = " ".join(process.info["cmdline"] or []).replace("\\", "/").lower()
            if name in {"r.exe", "rscript.exe", "rterm.exe", "rgui.exe"}:
                blocked.append({"pid": process.pid, "name": name, "reason": "existing_R_process"})
            elif "python" in name and "/fastprotflex/run_local.py" in command:
                blocked.append({"pid": process.pid, "name": name, "reason": "original_FastProtFlex_wrapper_still_active"})
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return blocked


def wait_gate(psutil, n_atoms, maximum_seconds):
    mib = 1024**2
    estimated_original_peak = 128.0 + (2944.76171875 - 128.0) * (n_atoms / 6620.0)**2
    required = int((1.15 * estimated_original_peak + 512.0) * mib)
    started = time.perf_counter()
    while True:
        jobs = competing_jobs(psutil)
        available = psutil.virtual_memory().available
        ready = not jobs and available >= required
        status = {"ready": ready, "available_MiB": available / mib, "required_available_MiB": required / mib,
                  "estimated_original_R_peak_MiB": estimated_original_peak, "blocking_jobs": jobs,
                  "elapsed_wait_seconds": time.perf_counter() - started, "maximum_wait_seconds": maximum_seconds}
        if ready or status["elapsed_wait_seconds"] >= maximum_seconds:
            return status
        print(json.dumps({"status": "waiting_for_original_memory_gate_and_no_competing_R", **status}), flush=True)
        time.sleep(min(30, maximum_seconds - status["elapsed_wait_seconds"]))


def session_versions(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    versions = sorted(set(re.findall(r"\b[A-Za-z][A-Za-z0-9.]*_[0-9][^\s]*", text)))
    first_line = text.splitlines()[0]
    return {"R_version_line": first_line, "package_version_tokens": versions}


def run_case(case, directory, args, psutil):
    gate = wait_gate(psutil, case["n_atoms"], args.memory_wait_seconds)
    if not gate["ready"]:
        result = {"status": "blocked_before_R_launch", "sample_id": case["sample_id"], "gate": gate}
        write_json(directory / "validation.json", result)
        return result
    directory.mkdir(parents=True, exist_ok=True)
    for source, target in ((HERE / "candidate/FUNCTION_GDV.r", "FUNCTION_GDV.r"),
                           (HERE / "wrapper/predict.r", "predict.r"),
                           (HERE / "wrapper/capture_full_precision.R", "capture_full_precision.R"),
                           (HERE / "capture_and_validate.R", "capture_and_validate.R"),
                           (BASELINE / "records" / case["sample_id"] / "protein.pdb", "protein.pdb")):
        shutil.copy2(source, directory / target)
        assert sha256(source) == sha256(directory / target)
    environment = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
                 "RCPP_PARALLEL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        environment[name] = "1"
    environment["R_LIBS_USER"] = ";".join(R_LIBRARIES)
    command = [str(RSCRIPT), "--vanilla", "capture_and_validate.R", case["baseline_RDS"]]
    started_utc = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter()
    peak_rss = 0
    peak_threads = 0
    cpu_times = {}
    recorded_children = {}
    timeout = False
    affinity_cpu = psutil.cpu_count() - 1
    with (directory / "stdout.log").open("w", encoding="utf-8") as stdout, (directory / "stderr.log").open("w", encoding="utf-8") as stderr:
        process = subprocess.Popen(command, cwd=directory, env=environment, stdout=stdout, stderr=stderr)
        monitor = psutil.Process(process.pid)
        while process.poll() is None:
            try:
                tree = [monitor] + monitor.children(recursive=True)
                rss = 0
                threads = 0
                for child in tree:
                    try:
                        child.cpu_affinity([affinity_cpu])
                        recorded_children[child.pid] = child
                        rss += child.memory_info().rss
                        threads += child.num_threads()
                        usage = child.cpu_times()
                        cpu_times[child.pid] = max(cpu_times.get(child.pid, 0.0), usage.user + usage.system)
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                peak_rss = max(peak_rss, rss)
                peak_threads = max(peak_threads, threads)
                if time.perf_counter() - started > args.case_timeout_seconds:
                    timeout = True
                    for child in reversed(list(recorded_children.values())):
                        try:
                            child.terminate()
                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            pass
                    process.terminate()
                    break
            except psutil.NoSuchProcess:
                pass
            time.sleep(0.1)
        return_code = process.wait(timeout=10)
    wall_seconds = time.perf_counter() - started
    comparison_path = directory / "comparison.json"
    comparison = json.loads(comparison_path.read_text(encoding="utf-8")) if comparison_path.exists() else {}
    actual_versions = session_versions(directory / "R_session_info.txt") if (directory / "R_session_info.txt").exists() else None
    baseline_versions = session_versions(Path(case["baseline_RDS"]).parent / "R_session_info.txt")
    versions_match = actual_versions == baseline_versions
    passed = return_code == 0 and not timeout and comparison.get("passed") is True and versions_match
    result = {"status": "complete_exact_all_atom_match" if passed else "failed_do_not_adopt",
              "sample_id": case["sample_id"], "started_utc": started_utc, "gate": gate,
              "command": command, "return_code": return_code, "timed_out": timeout,
              "wall_seconds_including_capture_and_RDS_comparison": wall_seconds,
              "observed_process_tree_CPU_seconds": sum(cpu_times.values()),
              "observed_max_single_process_CPU_seconds": max(cpu_times.values(), default=0.0),
              "peak_R_process_tree_RSS_bytes": peak_rss, "peak_R_process_tree_RSS_MiB": peak_rss / 1024**2,
              "reference_baseline_peak_MiB": case["baseline_R_process_peak_MiB"],
              "measured_peak_difference_MiB": peak_rss / 1024**2 - case["baseline_R_process_peak_MiB"],
              "cpu_affinity_logical_cpu": affinity_cpu, "library_thread_caps": 1,
              "sampling_interval_seconds": 0.1, "peak_OS_threads_including_runtime_helpers": peak_threads,
              "all_atom_RDS_comparison": comparison, "R_and_package_versions_match_baseline": versions_match,
              "actual_R_and_package_versions": actual_versions, "reference_R_and_package_versions": baseline_versions,
              "candidate_function_sha256": sha256(directory / "FUNCTION_GDV.r"),
              "input_sha256": sha256(directory / "protein.pdb"), "reference_RDS_sha256": sha256(Path(case["baseline_RDS"])),
              "candidate_RDS_sha256": sha256(directory / "atom_predictions_exact.rds") if (directory / "atom_predictions_exact.rds").exists() else None,
              "limitation": "Historical baseline peak was measured in a different run; final memory-gate changes require deliberate review of verified new measurements, never automatic bypass."}
    write_json(directory / "validation.json", result)
    return result


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true", help="Explicit future authorization to execute the three serial R validations; default never runs R")
    parser.add_argument("--run-tag", default="smoke_01")
    parser.add_argument("--memory-wait-seconds", type=float, default=600.0)
    parser.add_argument("--case-timeout-seconds", type=float, default=180.0)
    args = parser.parse_args()
    assert re.fullmatch(r"[a-z0-9][a-z0-9_-]{0,63}", args.run_tag)
    assert args.memory_wait_seconds >= 0 and args.case_timeout_seconds > 0
    provenance = static_checks()
    if not args.run:
        print(json.dumps({"status": "static_inspection_passed_R_not_executed", "candidate_sha256": provenance["candidate_function_sha256"],
                          "cases": [case["sample_id"] for case in provenance["cases"]],
                          "execution_requires": "--run after explicit root instruction; active original FastProtFlex wrapper/any R process and original RAM gate block launch"}, indent=2))
        return 0
    import psutil
    assert RSCRIPT.is_file()
    output = HERE / "runs" / args.run_tag
    assert not output.exists(), "Choose a new run-tag; never overwrite an existing validation measurement"
    output.mkdir(parents=True)
    results = []
    exit_code = 0
    for case in provenance["cases"]:
        result = run_case(case, output / case["sample_id"], args, psutil)
        results.append(result)
        print(json.dumps({"sample_id": case["sample_id"], "status": result["status"], "peak_MiB": result.get("peak_R_process_tree_RSS_MiB")}), flush=True)
        if result["status"] != "complete_exact_all_atom_match":
            exit_code = 3 if result["status"] == "blocked_before_R_launch" else 2
            break
    completed = len(results) == 3 and all(row["status"] == "complete_exact_all_atom_match" for row in results)
    write_json(output / "validation_audit.json", {"status": "three_sources_exact_all_atom_match" if completed else "blocked_or_failed_do_not_adopt",
               "created_utc": datetime.now(timezone.utc).isoformat(), "provenance": provenance,
               "validation_script_sha256": sha256(__file__), "R_validation_wrapper_sha256": sha256(HERE / "capture_and_validate.R"),
               "cases_requested": 3, "cases_attempted": len(results), "results": results,
               "candidate_adoption": "Not automatic. Even exact output equality and lower measured peak require root review before creating a distinct inference version or RAM schedule."})
    return 0 if completed else exit_code or 2


if __name__ == "__main__":
    sys.exit(main())
