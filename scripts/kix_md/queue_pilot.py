"""Four-state, replicate-1-only queue. Never schedules replicate 2 or above.

Uses the original MD runner/physical parameters and original checkpoint files.
The old queue must release the same queue.lock before this entry can proceed.
"""
import argparse
import datetime as dt
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys

STATES = ("apo", "myb", "mll", "both")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--consume-transition-stop", action="store_true",
                        help="Archive only the pilot transition stop, while holding queue lock")
    args = parser.parse_args()
    root = args.root.resolve(strict=True)
    import fcntl
    from run_md import run_job, write_json, now, sha, PARAMS
    lock = (root / "queue.lock").open("a+")
    fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
    policy_path = root / "PILOT_POLICY.json"
    policy = json.loads(policy_path.read_text())
    assert policy["allowed_repeats"] == [1] and policy["states"] == list(STATES)
    assert policy["additional_repeats_require_user_approval"] is True
    assert json.loads((root / "protocol.json").read_text()) == PARAMS
    assert (root / "smoke" / "PASSED.json").is_file()
    stamp = dt.datetime.now(dt.timezone.utc).strftime("%Y%m%dT%H%M%S%fZ")
    history = root / "queue_history" / stamp
    history.mkdir(parents=True)
    for name in ("queue.json", "progress.json"):
        if (root / name).exists():
            shutil.copy2(root / name, history / name)
    for state in STATES:
        previous_job = root / "runs" / state / "rep01" / "job.json"
        if previous_job.exists():
            shutil.copy2(previous_job, history / f"{state}_rep01_job.json")
    stop = root / "STOP_AFTER_SEGMENT"
    if stop.exists():
        if not args.consume_transition_stop:
            raise RuntimeError("Stop flag exists; explicit transition required, no MD launched")
        if sha(stop) != policy["transition_stop_sha256"]:
            raise RuntimeError("Stop flag is not the audited pilot transition; refusing to consume")
        # Recoverable move, no deletion; done ONLY after exclusive lock acquired.
        stop.rename(history / "STOP_AFTER_SEGMENT.consumed")
    elif args.consume_transition_stop:
        raise RuntimeError("Expected transition stop absent; inspect instead of guessing")
    queue = dict(pid=os.getpid(), started_utc=now(), stage="pilot_rep01_only",
                 cuda_visible_devices=os.environ.get("CUDA_VISIBLE_DEVICES"),
                 order=[f"{state}/rep01" for state in STATES],
                 allowed_repeats=[1], deferred_repeats=[2, 3, 4, 5],
                 worker_count=1, code_sha256=sha(__file__), md_runner_sha256=sha(root / "code" / "run_md.py"),
                 policy_sha256=sha(policy_path), transition_history=str(history),
                 original_full_protocol_preserved=True,
                 active_total_production_ns=400, active_total_including_equilibration_ns=440)
    write_json(root / "queue.json", queue)
    print(json.dumps(queue), flush=True)
    for state in STATES:
        if stop.exists():
            raise SystemExit(75)
        write_json(root / "progress.json", dict(status="starting_or_resuming", stage="pilot_rep01_only",
                   state=state, repeat=1, updated_utc=now()))
        run_job(root, state, 1)
        folder = root / "runs" / state / "rep01"
        done = json.loads((folder / "DONE.json").read_text())
        assert done["repeat"] == 1 and done["production_ns"] == 100 and done["final_step"] == 55000000
        summary = folder / "analysis" / "summary.json"
        if not summary.exists():
            write_json(root / "progress.json", dict(status="analyzing_single_repeat", state=state,
                       repeat=1, stage="pilot_rep01_only", updated_utc=now()))
            subprocess.run([sys.executable, str(root / "code" / "analyze_job.py"),
                            "--protein-pdb", str(root / "prepared" / state / "protein.pdb"),
                            "--job-dir", str(folder)], check=True)
        info = json.loads(summary.read_text())
        assert info["n_frames"] == 100000 and info["state"] == state and info["replicate"] == "rep01"
        assert info["output_csv_sha256"] == sha(folder / "analysis" / "residue_metrics.csv")
        # Update descriptive model comparisons as each state becomes available.
        subprocess.run([sys.executable, str(root / "code" / "compare_pilot.py"), "--root", str(root)], check=True)
    write_json(root / "progress.json", dict(status="pilot_complete_awaiting_review", stage="pilot_rep01_only",
               finished_utc=now(), jobs=4, production_ns=400, total_including_equilibration_ns=440,
               deferred_repeats=[2, 3, 4, 5], additional_repeats_require_user_approval=True))
    print("PILOT COMPLETE. All four first repeats analyzed. No additional MD is scheduled.", flush=True)


if __name__ == "__main__":
    main()
