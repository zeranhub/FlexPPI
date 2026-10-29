"""Complete all four audited absolute predictions, one process at a time."""
from pathlib import Path
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parent
for pdb in ["3hg1", "4jff", "4jfd", "4jfe"]:
    output = ROOT / "raw_outputs/proaffinity" / pdb
    output.mkdir(parents=True, exist_ok=True)
    audit = output / "audit.json"
    if audit.exists():
        record = json.loads(audit.read_text(encoding="utf-8"))
        assert record["upstream_bug_fix"] and len(record["graph_checks"]) == 3
        print("Reuse verified", pdb, record["prediction_as_reported_pKa"], flush=True)
        continue
    with (output / "run.log").open("w", encoding="utf-8") as log:
        result = subprocess.run([sys.executable, "-u", "-X", "utf8",
                                 str(ROOT / "run_proaffinity.py"), pdb],
                                stdout=log, stderr=subprocess.STDOUT, timeout=900)
    if result.returncode:
        print((output / "run.log").read_text(encoding="utf-8"), flush=True)
        raise SystemExit(result.returncode)
    record = json.loads(audit.read_text(encoding="utf-8"))
    print("Complete", pdb, record["prediction_as_reported_pKa"], flush=True)
