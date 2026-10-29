"""One CPU-only method/round in a fresh process; never train or score labels."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback

# Enforce, rather than inherit, the common deployment thread/device settings.
os.environ["CUDA_VISIBLE_DEVICES"] = ""
os.environ["WANDB_MODE"] = "disabled"
for key in ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS"]:
    os.environ[key] = "2"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, data):
    Path(path).write_text(json.dumps(data, indent=2, ensure_ascii=True), encoding="utf-8")


def host_snapshot():
    data = dict(utc=datetime.now(timezone.utc).isoformat(), platform=platform.platform(),
        machine=platform.machine(), hostname=platform.node(), processor=platform.processor(),
        logical_cpus=os.cpu_count(), python=sys.version, executable=sys.executable,
        thread_environment={k: os.environ.get(k) for k in
            ["OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS", "NUMEXPR_NUM_THREADS", "CUDA_VISIBLE_DEVICES"]})
    if hasattr(os, "getloadavg"):
        data["host_load_1_5_15_min"] = list(os.getloadavg())
    if Path("/proc/cpuinfo").exists():
        lines = Path("/proc/cpuinfo").read_text().splitlines()
        data["cpu_model_names"] = sorted({x.split(":", 1)[1].strip() for x in lines if x.startswith("model name")})
    if Path("/proc/meminfo").exists():
        data["memory_kib"] = {x.split(":", 1)[0]: x.split(":", 1)[1].strip()
            for x in Path("/proc/meminfo").read_text().splitlines() if x.startswith(("MemTotal:", "MemAvailable:"))}
    try:
        result = subprocess.run(["nvidia-smi", "--query-gpu=name,utilization.gpu,memory.used,memory.total",
            "--format=csv,noheader"], capture_output=True, text=True, timeout=10, check=False)
        data["background_gpu_readonly_snapshot"] = dict(exit_code=result.returncode, output=result.stdout.strip())
    except (OSError, subprocess.TimeoutExpired) as error:
        data["background_gpu_readonly_snapshot"] = {"unavailable": str(error)}
    return data


def preflight(plan, method):
    if plan["cpu_threads"] != 2 or plan["repetitions"] != 3 or plan["batch_size"] != 1:
        raise ValueError("Deployment cost plan contract changed")
    if len(plan["records"]) != 32 or plan["family_components"] != 1:
        raise ValueError("External deployment cohort changed")
    if plan["records"] != sorted(plan["records"], key=lambda r: r["sample_id"]):
        raise ValueError("Input order must be fixed before timing")
    for r in plan["records"]:
        if set(r) != {"sample_id", "structure_path", "target_chains", "metadata", "input_sha256", "chain_count"}:
            raise ValueError("Deployment input must expose no label/trajectory/features")
        if set(r["metadata"]) != {"residue_count"}:
            raise ValueError("Unexpected deployment metadata")
    for path, expected in plan["files_sha256"].items():
        if sha(path) != expected:
            raise ValueError("Deployment preflight hash changed: " + path)
    spec = plan["methods"][method]
    for path in spec.get("checkpoint_paths", []):
        if path not in plan["files_sha256"]:
            raise ValueError("Final checkpoint absent from hashed plan: " + path)
    return spec


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--method", required=True)
    parser.add_argument("--round", type=int, default=0)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--smoke-only", action="store_true")
    parser.add_argument("--audit-only", action="store_true")
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=False)
    started_utc = datetime.now(timezone.utc).isoformat()
    report = dict(status="running", method=args.method, round=args.round,
        started_utc=started_utc, worker_sha256=sha(__file__),
        final_cost_results=False, smoke_only=args.smoke_only, audit_only=args.audit_only,
        new_fits=0, accuracy_evaluation=False, ground_truth_used=False, cache_hit=False,
        timer="time.perf_counter", device="cpu", host_before=host_snapshot())
    write_json(out / "worker_state.json", report)
    try:
        plan = json.loads(args.plan.read_text(encoding="utf-8"))
        report.update(plan_sha256=sha(args.plan), scope=plan["scope"],
            serialization_scope=plan["serialization_scope"],
            backflip2_identity_contract=plan["backflip2_identity_contract"])
        check_start = time.perf_counter()
        spec = preflight(plan, args.method)
        report["file_preflight_seconds"] = time.perf_counter() - check_start
        os.environ["FLEXPPI_DEPLOYMENT_PAYLOAD_ROOT"] = plan["primary_root"]
        if spec["software_environment"] == "neural":
            for root in reversed(plan.get("dependency_roots", [])):
                sys.path.insert(0, root)
        import_start = time.perf_counter()
        import numpy as np
        import torch
        from threadpoolctl import threadpool_limits, threadpool_info
        torch.set_num_threads(2)
        torch.set_num_interop_threads(1)
        if spec["kind"] == "backflip2":
            module = importlib.import_module("deployment_cost_backflip_adapter")
            dependencies = module.import_dependencies(spec["source_root"])
        else:
            module = importlib.import_module("deployment_cost_neural_static_adapter")
        report["dependency_import_seconds"] = time.perf_counter() - import_start
        context_start = time.perf_counter()
        if spec["kind"] == "backflip2":
            context = {**spec, "dependencies": dependencies, "records": plan["identities"]}
        else:
            context = module.prepare_context(args.method,
                checkpoint_paths=spec.get("checkpoint_paths"), config_paths=spec.get("config_paths"),
                static_paths=spec.get("static_paths"), schema_path=spec.get("schema_path"))
        report["configuration_preflight_seconds"] = time.perf_counter() - context_start
        versions = {name: importlib.metadata.version(name) for name in
            ["numpy", "torch", "scipy", "threadpoolctl"]}
        for name in ["scikit-learn", "joblib", "omegaconf", "biotite"]:
            try:
                versions[name] = importlib.metadata.version(name)
            except importlib.metadata.PackageNotFoundError:
                versions[name] = None
        report["package_versions"] = versions
        report["torch_threads"] = torch.get_num_threads()
        report["torch_interop_threads"] = torch.get_num_interop_threads()
        if args.audit_only:
            report.update(status="preflight_complete_no_model_execution", host_after=host_snapshot())
            write_json(out / "worker_state.json", report)
            print(report["status"], args.method, flush=True)
            return
        with threadpool_limits(limits=2):
            load_start = time.perf_counter()
            adapter = (module.BackFlipAdapter(context) if spec["kind"] == "backflip2"
                       else module.NeuralStaticAdapter(args.method, context))
            report["model_load_seconds"] = time.perf_counter() - load_start
            report["adapter_metadata"] = adapter.metadata() if callable(adapter.metadata) else adapter.metadata
            report["actual_threadpools"] = threadpool_info()
            if torch.get_num_threads() != 2 or torch.get_num_interop_threads() != 1:
                raise ValueError("Runtime thread limit was not applied")
            pred_dir = out / "predictions"
            pred_dir.mkdir()
            chain_log = out / "backflip2_chain_diagnostics.jsonl"
            def call(record, phase):
                path = pred_dir / (("warmup_" if phase == "warmup" else "") + record["sample_id"] + ".npz")
                if path.exists():
                    raise ValueError("Cost outputs cannot be overwritten or used as cache hits")
                begun = time.perf_counter()
                values = adapter.predict_record(record)
                identity = plan["identities"][record["sample_id"]]
                n = record["metadata"]["residue_count"]
                for key in ["amino_acid", "chain", "author_residue_id", "geometry_mask", "target_mask"]:
                    if not np.array_equal(np.asarray(values[key]), np.asarray(identity[key])):
                        raise ValueError("Timed prediction changed residue identity/mask: " + key)
                mean = np.asarray(values["mean"])
                if mean.shape != (n, len(spec["output_metrics"])) or not np.isfinite(mean).all():
                    raise ValueError("Incomplete or nonfinite deployment output")
                mask = np.asarray(values["output_mask"], dtype=bool)
                if mask.shape != (n,):
                    raise ValueError("Output geometry mask has wrong length")
                if "model_variance" in values:
                    var = np.asarray(values["model_variance"])
                    if var.shape != mean.shape or not np.isfinite(var).all() or (var < 0).any():
                        raise ValueError("Invalid uncalibrated predictive variance")
                arrays = {k: np.asarray(values[k]) for k in
                    ["amino_acid", "chain", "author_residue_id", "geometry_mask", "target_mask", "output_mask"]}
                arrays["rmsf_angstrom"] = mean[:, 0]
                if len(spec["output_metrics"]) == 2:
                    arrays["chi_entropy_nats"] = mean[:, 1]
                # The common timed export saves predicted metrics and identities;
                # uncertainty is computed/converted but not serialized here.
                np.savez_compressed(path, **arrays)
                elapsed = time.perf_counter() - begun
                row = dict(method=args.method, round=args.round, phase=phase,
                    sample_id=record["sample_id"], residues=n, chains=record["chain_count"],
                    seconds=elapsed, cpu_threads=2, device="cpu", cache_hit=False,
                    output_metrics=spec["output_metrics"], output_path=str(path),
                    output_sha256=sha(path), output_bytes=path.stat().st_size,
                    serialization_keys=list(arrays), output_dtypes={k: str(v.dtype) for k, v in arrays.items()})
                if spec["kind"] == "backflip2":
                    with chain_log.open("a", encoding="utf-8") as handle:
                        handle.write(json.dumps({"phase": phase, "sample_id": record["sample_id"],
                            "chains": adapter.last_chain_timings}) + "\n")
                return row
            warmup = call(plan["records"][0], "warmup")
            report["warmup"] = warmup
            rows = []
            records = plan["records"][:1] if args.smoke_only else plan["records"]
            with (out / "calls.jsonl").open("x", encoding="utf-8") as handle:
                for record in records:
                    row = call(record, "smoke" if args.smoke_only else "measurement")
                    rows.append(row)
                    handle.write(json.dumps(row) + "\n")
                    handle.flush()
                    print(args.method, args.round, record["sample_id"], "complete", flush=True)
            report.update(status="smoke_complete_not_benchmark" if args.smoke_only else "complete_round",
                measured_records=len(rows), measured_seconds_sum=sum(r["seconds"] for r in rows),
                final_cost_results=False, host_after=host_snapshot(),
                completed_utc=datetime.now(timezone.utc).isoformat())
            write_json(out / "worker_state.json", report)
    except (Exception, KeyboardInterrupt) as error:
        report.update(status="failed", error=repr(error), traceback=traceback.format_exc(),
            completed_utc=datetime.now(timezone.utc).isoformat(), host_after=host_snapshot())
        write_json(out / "worker_state.json", report)
        raise


if __name__ == "__main__":
    main()
