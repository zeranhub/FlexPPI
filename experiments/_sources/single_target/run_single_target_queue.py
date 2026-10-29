"""Serial, resumable single-GPU queue for the ten prespecified final-epoch fits."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from completion_common import ARM_METRICS, sha256, write_json


def valid_final(path, job):
    import torch
    import yaml
    checkpoint = torch.load(path, map_location="cpu", weights_only=False)
    config = checkpoint["config"]
    expected = yaml.safe_load(Path(job["config"]).read_text(encoding="utf-8"))
    for key in ["stage", "seed", "device", "output_dir", "init_checkpoint"]:
        if config[key] != expected[key]:
            raise ValueError("Existing checkpoint differs in " + key)
    for section in ["data", "model", "train"]:
        for key, value in expected[section].items():
            if config[section][key] != value:
                raise ValueError("Existing checkpoint differs in " + section + "." + key)
    if checkpoint["epoch"] != 40 or config["seed"] != 7:
        raise ValueError("Existing final checkpoint has the wrong epoch/seed")
    if config["init_checkpoint"] is not None or config["model"]["metrics"] != ARM_METRICS[job["arm"]]:
        raise ValueError("Existing checkpoint is not the requested from-scratch target")
    if config["model"]["information_ablation"] != "full" or config["train"]["sampling_strategy"] != "source_system_balanced":
        raise ValueError("Existing checkpoint differs from matched full protocol")
    if Path(config["output_dir"]).resolve() != Path(job["output_dir"]).resolve():
        raise ValueError("Existing checkpoint output directory mismatch")
    return sha256(path)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--job-limit", type=int, default=0, help="0 runs every outstanding job")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    plan = json.loads(args.plan.read_text(encoding="utf-8"))
    training, output = Path(plan["training_root"]), Path(plan["output_root"])
    for relative, expected in {**plan["input_file_hashes"], **plan["source_hashes"]}.items():
        if sha256(training / relative) != expected:
            raise ValueError("Frozen training artifact differs: " + relative)
    env = dict(os.environ)
    env["PYTHONPATH"] = str(training / "src")
    env["OMP_NUM_THREADS"] = "2"
    env["MKL_NUM_THREADS"] = "2"
    state_path = output / "single_target_queue_state.json"
    state = dict(protocol=plan["protocol"], started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), jobs=[])
    completed_now = 0
    for job in plan["jobs"]:
        if sha256(job["config"]) != job["config_sha256"]:
            raise ValueError("Completion configuration changed: " + job["config"])
        final = Path(job["final_checkpoint"])
        row = dict(job)
        command = [sys.executable, "-u", "-m", "flexppi", "train", "--config", job["config"]]
        print(" ".join(command), flush=True)
        if args.dry_run:
            continue
        if final.exists():
            row.update(status="complete_existing", checkpoint_sha256=valid_final(final, job))
        else:
            if args.job_limit and completed_now >= args.job_limit:
                row["status"] = "pending"
                state["jobs"].append(row)
                continue
            log = output / "logs" / ("{}_fold{}.log".format(job["arm"], job["fold"]))
            log.parent.mkdir(parents=True, exist_ok=True)
            row.update(status="running", started_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
            state["jobs"].append(row)
            write_json(state_path, state)
            with log.open("ab") as handle:
                process = subprocess.run(command, cwd=training, env=env, stdout=handle, stderr=subprocess.STDOUT)
            row.update(exit_code=process.returncode, finished_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()))
            if process.returncode or not final.exists():
                row["status"] = "failed"
                write_json(state_path, state)
                raise RuntimeError("Fit failed; inspect " + str(log))
            row.update(status="complete", checkpoint_sha256=valid_final(final, job))
            completed_now += 1
            write_json(state_path, state)
            continue
        state["jobs"].append(row)
        write_json(state_path, state)
    if not args.dry_run:
        state["status"] = "complete" if all(row["status"].startswith("complete") for row in state["jobs"]) else "partial"
        write_json(state_path, state)


if __name__ == "__main__":
    main()
