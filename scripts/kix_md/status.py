"""Read-only status snapshot; never launches, terminates, or alters a simulation."""
import argparse
import csv
import datetime as dt
import json
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()
    policy_file = root / "PILOT_POLICY.json"
    policy = json.loads(policy_file.read_text()) if policy_file.exists() else None
    allowed_repeats = policy["allowed_repeats"] if policy else list(range(1, 6))
    rows = []
    for repeat in range(1, 6):
        for state in ["apo", "myb", "mll", "both"]:
            folder = root / "runs" / state / f"rep{repeat:02d}"
            row = dict(state=state, repeat=repeat, status="pending" if repeat in allowed_repeats else "deferred_until_review",
                       equilibration_ns_committed=0, production_ns_committed=0)
            if folder.exists():
                row["status"] = "started"
                for phase in ["equilibration", "production"]:
                    row[f"{phase}_ns_committed"] = len(list((folder / phase).glob("segment_*.json")))
                if (folder / "DONE.json").exists():
                    row["status"] = "md_completed"
                row["analysis_completed"] = (folder / "analysis" / "summary.json").exists()
                logs = sorted(folder.glob("*/segment_*.csv"), key=lambda p: p.stat().st_mtime_ns)
                if logs:
                    with logs[-1].open() as stream:
                        entries = list(csv.DictReader(stream))
                    if entries:
                        row["latest_energy_log"] = dict(file=str(logs[-1].relative_to(root)), values=entries[-1])
            rows.append(row)
    queue_path = root / "queue.json"
    queue = json.loads(queue_path.read_text()) if queue_path.exists() else None
    running = False
    if queue:
        process_command = Path(f"/proc/{queue['pid']}/cmdline")
        if process_command.exists():
            command = process_command.read_bytes().replace(b"\x00", b" ").decode()
            running = str(root) in command and ("run_md.py queue" in command or "queue_pilot.py" in command)
    result = dict(snapshot_utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                  root=str(root), queue=queue, queue_process_verified_running=running,
                  active_policy=policy,
                  preparation_completed=(root / "prepared" / "READY.json").exists(),
                  smoke_tests_passed=(root / "smoke" / "PASSED.json").exists(),
                  completed_jobs=sum(r["status"] == "md_completed" for r in rows),
                  production_ns_committed=sum(r["production_ns_committed"] for r in rows),
                  jobs=rows)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
