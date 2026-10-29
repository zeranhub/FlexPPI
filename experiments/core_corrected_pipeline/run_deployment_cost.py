"""Serial, append-only three-round benchmark after both frozen queues finish.

No training, score selection, resume-cache timing, old-weight substitution,
or GPU operations occur. Failed rounds/logs are retained as incomplete runs.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import sys
import time
import traceback
import uuid

HERE = Path(__file__).resolve().parent


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, ensure_ascii=True), encoding="utf-8")


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def percentile(values, fraction):
    ordered = sorted(values)
    position = (len(ordered) - 1) * fraction
    lo, hi = math.floor(position), math.ceil(position)
    return ordered[lo] + (position - lo) * (ordered[hi] - ordered[lo])


def readiness(plan):
    if plan["status"] != "ready_for_benchmark" or plan["missing_final_weights"]:
        raise ValueError("Final frozen weights pending; deployment benchmark deferred")
    for path, expected in plan["files_sha256"].items():
        if sha(path) != expected:
            raise ValueError("Prepared input/model/code changed: " + path)
    for path in plan["completion_gates"]:
        if not Path(path).is_file() or load(path).get("status") != "complete":
            raise ValueError("Both frozen training queues and external exports must complete first: " + path)
    if platform.system() != "Linux" or not Path("/proc/cpuinfo").is_file():
        raise ValueError("The final common-host benchmark must run on hit-gpu, not this host")
    if "Gold 6342" not in Path("/proc/cpuinfo").read_text():
        raise ValueError("The predeclared Xeon Gold 6342 host is required")
    for executable in plan["python_executables"].values():
        if not Path(executable).is_file():
            raise ValueError("Declared method environment absent: " + executable)


def write_csv(path, rows):
    with Path(path).open("x", encoding="utf-8", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def summarize(run, plan, worker_reports):
    raw, summaries, medians = [], [], []
    inputs = {r["sample_id"]: r for r in plan["records"]}
    hosts = {(w["host_before"]["hostname"], tuple(w["host_before"].get("cpu_model_names", []))) for w in worker_reports}
    if len(hosts) != 1:
        raise ValueError("All timing workers must use the same actual host/CPU")
    for protocol, names in [("primary", plan["primary_methods"]), ("strict_single", plan["strict_methods"])]:
        for method in names:
            calls = []
            round_totals = []
            initializations = []
            for repetition in range(3):
                directory = run / protocol / f"round_{repetition}" / method
                state = load(directory / "worker_state.json")
                if state["status"] != "complete_round" or state["measured_records"] != 32:
                    raise ValueError("Partial or smoke outputs cannot enter final costs")
                rows = [json.loads(x) for x in (directory / "calls.jsonl").read_text().splitlines()]
                if [r["sample_id"] for r in rows] != list(inputs) or len(rows) != 32:
                    raise ValueError("Timing round does not cover the complete fixed input order")
                for r in rows:
                    if r["cache_hit"] or r["phase"] != "measurement" or r["method"] != method or r["round"] != repetition:
                        raise ValueError("Cached, misidentified or partial timing row")
                    if not math.isfinite(r["seconds"]) or r["seconds"] <= 0:
                        raise ValueError("Invalid recorded runtime")
                    if sha(r["output_path"]) != r["output_sha256"]:
                        raise ValueError("Timed output changed after collection")
                    input_record = inputs[r["sample_id"]]
                    if r["residues"] != input_record["metadata"]["residue_count"] or r["chains"] != input_record["chain_count"]:
                        raise ValueError("Timing input size changed")
                    calls.append(r)
                    raw.append(dict(protocol=protocol, method=method, round=repetition,
                        sample_id=r["sample_id"], residues=r["residues"], chains=r["chains"],
                        seconds=r["seconds"], cpu_threads=r["cpu_threads"], device=r["device"],
                        cache_hit=False, output_bytes=r["output_bytes"], output_sha256=r["output_sha256"],
                        output_metrics=";".join(r["output_metrics"]), output_path=r["output_path"]))
                round_totals.append(sum(r["seconds"] for r in rows))
                initializations.append({"round": repetition, "dependency_import_seconds": state["dependency_import_seconds"],
                    "model_load_seconds": state["model_load_seconds"], "warmup_seconds": state["warmup"]["seconds"],
                    "worker_state": str(directory / "worker_state.json")})
            record_medians = []
            for sample_id, r in inputs.items():
                times = [x["seconds"] for x in calls if x["sample_id"] == sample_id]
                if len(times) != 3:
                    raise ValueError("Exactly three complete measurements required per record")
                median = statistics.median(times)
                record_medians.append(median)
                medians.append(dict(protocol=protocol, method=method, sample_id=sample_id,
                    residues=r["metadata"]["residue_count"], chains=r["chain_count"],
                    round_0_seconds=times[0], round_1_seconds=times[1], round_2_seconds=times[2],
                    median_seconds=median, maximum_seconds=max(times)))
            weight_paths = plan["methods"][method].get("checkpoint_paths", [])
            metadata = load(run / protocol / "round_0" / method / "worker_state.json")["adapter_metadata"]
            if method == "backflip2":
                weight_paths = [plan["methods"][method]["checkpoint_path"]]
            elif plan["methods"][method]["kind"] == "static":
                weight_paths = [p for member in plan["methods"][method]["static_paths"] for p in member.values()]
            summaries.append(dict(protocol=protocol, method=method, records=32, repetitions=3,
                output_metrics=plan["methods"][method]["output_metrics"],
                saved_checkpoint_files_bytes=sum(Path(p).stat().st_size for p in weight_paths),
                checkpoint_bytes_scope="Serialized training checkpoint(s), including optimizer/config when present; not pure parameters or inference RAM",
                pure_parameter_bytes=(metadata["parameter_bytes"] if method == "backflip2" else
                    sum(metadata["parameter_bytes_per_member"]) if plan["methods"][method]["kind"] == "neural" else None),
                record_median_seconds=statistics.median(record_medians),
                record_p95_seconds=percentile(record_medians, .95),
                maximum_record_median_seconds=max(record_medians),
                longest_actual_call_seconds=max(r["seconds"] for r in calls),
                all_96_calls_p95_seconds=percentile([r["seconds"] for r in calls], .95),
                sum_of_record_medians_seconds=sum(record_medians),
                round_total_seconds=round_totals, median_round_total_seconds=statistics.median(round_totals),
                initializations=initializations, independent_family_ci=False, p_value=None))
    write_csv(run / "calls_all_rounds.csv", raw)
    write_csv(run / "record_median_cost.csv", medians)
    write_json(run / "cost_summary.json", summaries)
    return summaries


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, default=HERE / "deployment_cost/runs")
    args = parser.parse_args()
    plan = load(args.plan)
    # Readiness is read-only, before claiming/starting a benchmark.
    readiness(plan)
    root = args.output_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    lock_path = root / "benchmark.lock"
    run = root / (datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8])
    report = dict(status="running", started_utc=datetime.now(timezone.utc).isoformat(),
        plan_sha256=sha(args.plan), runner_sha256=sha(__file__), new_fits=0,
        accuracy_evaluation=False, final_cost_results=False, completed_workers=[])
    wall_start = time.perf_counter()
    lock_fd = os.open(lock_path, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    try:
        with os.fdopen(lock_fd, "w", encoding="utf-8") as handle:
            json.dump(dict(pid=os.getpid(), host=platform.node(), started_utc=datetime.now(timezone.utc).isoformat()), handle)
        run.mkdir(exist_ok=False)
        write_json(run / "plan.json", plan)
        write_json(run / "run_state.json", report)
        workers = []
        for protocol, names in [("primary", plan["primary_methods"]), ("strict_single", plan["strict_methods"])]:
            for repetition in range(3):
                order = names[repetition:] + names[:repetition]
                for method in order:
                    directory = run / protocol / f"round_{repetition}" / method
                    directory.parent.mkdir(parents=True, exist_ok=True)
                    log = directory.parent / (method + ".log")
                    environment = plan["methods"][method]["software_environment"]
                    executable = plan["python_executables"][environment]
                    command = [executable, "-u", str(HERE / "deployment_cost_worker.py"),
                        "--plan", str(run / "plan.json"), "--method", method,
                        "--round", str(repetition), "--out", str(directory)]
                    print(protocol, repetition, method, "starting", flush=True)
                    with log.open("x", encoding="utf-8") as output:
                        result = subprocess.run(command, stdout=output, stderr=subprocess.STDOUT, check=False)
                    if result.returncode:
                        raise RuntimeError("Timing worker failed; retained its log and partial round: " + str(log))
                    state = load(directory / "worker_state.json")
                    workers.append(state)
                    report["completed_workers"].append(dict(protocol=protocol, round=repetition, method=method,
                        worker_state=str(directory / "worker_state.json"), exit_code=result.returncode))
                    write_json(run / "run_state.json", report)
        summaries = summarize(run, plan, workers)
        report.update(status="complete", final_cost_results=True,
            completed_utc=datetime.now(timezone.utc).isoformat(), run_wall_seconds=time.perf_counter() - wall_start,
            measured_methods=len(summaries), actual_host=platform.node(),
            scope=plan["scope"], training_accuracy_complete=False,
            statistics="All three rounds retained; per-record median plus full round totals. No CI/p value or fastest-call selection.")
        write_json(run / "run_state.json", report)
        hashes = {str(p.relative_to(run)): sha(p) for p in run.rglob("*") if p.is_file()}
        write_json(run / "output_manifest.json", {"files_sha256": hashes})
        print("Complete genuine deployment-cost run:", run, flush=True)
    except (Exception, KeyboardInterrupt) as error:
        report.update(status="failed_incomplete_cost_run", error=repr(error), traceback=traceback.format_exc(),
            final_cost_results=False, completed_utc=datetime.now(timezone.utc).isoformat())
        if run.exists():
            write_json(run / "run_state.json", report)
        else:
            write_json(root / (run.name + "_initialization_failure.json"), report)
        raise
    finally:
        # This is the lock created by this process, never another user's file.
        lock_path.unlink()


if __name__ == "__main__":
    main()
