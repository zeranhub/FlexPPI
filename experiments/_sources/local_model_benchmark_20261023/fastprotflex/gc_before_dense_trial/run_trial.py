"""One isolated6620-atom trial: only one extra gc before first dense matrix."""
import argparse
from datetime import datetime, timezone
import difflib
import hashlib
import importlib.util
import json
from pathlib import Path
import re
import shutil
from types import SimpleNamespace

HERE = Path(__file__).resolve().parent
BASELINE = HERE.parent
V2 = BASELINE / "memory_lifetime_v2"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run",action="store_true")
    args = parser.parse_args()
    provenance = json.loads((V2 / "lifecycle_provenance.json").read_text(encoding="utf-8"))
    case = next(row for row in provenance["cases"] if row["sample_id"] == "itrace_1qrn")
    assert case["n_atoms"] == 6620 and sha(case["baseline_RDS"]) == case["baseline_RDS_sha256"]
    source = V2 / "candidate/FUNCTION_GDV.r"
    original = source.read_bytes()
    assert sha(source) == provenance["candidate_function_sha256"]
    needle = b"  adj<-as.matrix(DD)"
    newline = b"\r\n" if b"\r\n" in original else b"\n"
    trial = original.replace(needle,b"  gc()"+newline+needle,1)
    original_lines,trial_lines = original.decode("utf-8").splitlines(),trial.decode("utf-8").splitlines()
    differences = list(difflib.ndiff(original_lines,trial_lines))
    assert [line for line in differences if line.startswith("+ ")] == ["+   gc()"]
    assert not any(line.startswith("- ") for line in differences)
    strip = lambda lines:[line for line in lines if not re.fullmatch(r"\s*(?:rm\([^;]*\)|gc\(\))\s*",line)]
    assert strip(original_lines) == strip(trial_lines)
    assert not (HERE / "candidate").exists(), "One trial only; do not overwrite earlier sources/results"
    (HERE / "candidate").mkdir()
    (HERE / "control").mkdir()
    (HERE / "wrapper").mkdir()
    (HERE / "candidate/FUNCTION_GDV.r").write_bytes(trial)
    (HERE / "control/FUNCTION_GDV_v2.r").write_bytes(original)
    for filename in ("predict.r","capture_full_precision.R"):
        shutil.copy2(V2 / "wrapper" / filename,HERE / "wrapper" / filename)
    shutil.copy2(V2 / "capture_and_validate.R",HERE / "capture_and_validate.R")
    (HERE / "FUNCTION_GDV_only_pre_dense_gc.diff").write_text("".join(difflib.unified_diff(original.decode("utf-8").splitlines(True),trial.decode("utf-8").splitlines(True),fromfile="v2/FUNCTION_GDV.r",tofile="trial/FUNCTION_GDV.r")),encoding="utf-8")
    protected = [source,V2 / "lifecycle_provenance.json",BASELINE / "run_local_lifetime_v2.py",V2 / "probe_extreme_resources.py",
                 BASELINE / "run_local.py",BASELINE / "capture_full_precision.R",Path(provenance["official_function_path"]),
                 V2 / "runs/smoke_01/validation_audit.json"]
    before = {str(path):sha(path) for path in protected}
    preparation = {"status":"isolated_single_case_prepared", "created_utc":datetime.now(timezone.utc).isoformat(),
                   "sample_id":"itrace_1qrn","n_atoms":6620,"only_change":"Insert gc() immediately before first adj<-as.matrix(DD)",
                   "control_v2_source_sha256":sha(source),"trial_source_sha256":sha(HERE / "candidate/FUNCTION_GDV.r"),
                   "run_script_sha256":sha(__file__),"original_reference_RDS_sha256":case["baseline_RDS_sha256"],
                   "numeric_program_text_identical":True,"protected_files_sha256":before,
                   "resource_gate":"Unchanged original6620 smoke gate:1.15*2944.76171875+512MiB=3898.4759759902954MiB", "R_runs_requested":1}
    (HERE / "trial_preparation.json").write_text(json.dumps(preparation,indent=2),encoding="utf-8")
    if not args.run:
        print(json.dumps(preparation,indent=2))
        return 0
    # Reuse the existing exact-RDS comparison/serial monitor without editing it.
    spec = importlib.util.spec_from_file_location("unchanged_single_case_validator",V2 / "validate_candidate.py")
    validator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(validator)
    validator.HERE,validator.BASELINE = HERE,BASELINE
    import psutil
    result = validator.run_case(case,HERE / "output/itrace_1qrn",SimpleNamespace(memory_wait_seconds=180.0,case_timeout_seconds=180.0),psutil)
    assert before == {str(path):sha(path) for path in protected}
    smoke = json.loads((V2 / "runs/smoke_01/validation_audit.json").read_text(encoding="utf-8"))
    control_peak = next(row["peak_R_process_tree_RSS_MiB"] for row in smoke["results"] if row["sample_id"] == case["sample_id"])
    peak = result.get("peak_R_process_tree_RSS_MiB")
    summary = {"created_utc":datetime.now(timezone.utc).isoformat(), "preparation":preparation,"result":result,
               "control_v2_actual_peak_MiB":control_peak,"trial_peak_MiB":peak,
               "difference_from_control_v2_MiB":None if peak is None else peak-control_peak,
               "below2500MiB":None if peak is None else peak<2500,
               "all_protected_files_unchanged":True,"scope":"One isolated measured trial; no candidate adoption or new three-source equivalence claim."}
    (HERE / "trial_result.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    print(json.dumps(summary,indent=2),flush=True)
    return 0 if result["status"] == "complete_exact_all_atom_match" else 3 if result["status"] == "blocked_before_R_launch" else 2


if __name__ == "__main__":
    raise SystemExit(main())
