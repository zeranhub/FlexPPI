"""Freeze a visually checked binary benchmark and its classification evidence."""
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import shutil

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT.parent
DESTINATION = SOURCE / "snapshots/completed_14_protein_complex"


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def main():
    qa_path = ROOT / "QA_REPORT.json"
    qa = json.loads(qa_path.read_text(encoding="utf-8"))
    assert qa["status"] == "PASS_CURRENT_SNAPSHOT_WITH_DOCUMENTED_VALIDATOR_LIMITATIONS"
    assert qa["passed_custom_checks"] == qa["total_custom_checks"]
    assert all(row["status"] == "inspected_current_hash_with_notes" for row in qa["visual_reviews"])
    scoring = json.loads((ROOT / "evaluation/scoring_audit.json").read_text(encoding="utf-8"))
    assert len(scoring["complete_methods"]) == 14 and scoring["group_counts"] == {"protein": 2203, "protein_complex": 860}
    assert scoring["unresolved_records"] == 0
    for path, digest in qa["input_sha256"].items():
        assert sha(path) == digest, f"Changed after QA: {path}"
    classification = SOURCE / "grouping_protein_complex_20261024"
    membership_audit = json.loads((classification / "membership_freeze_audit.json").read_text(encoding="utf-8"))
    assert membership_audit["membership_sha256"] == scoring["membership_sha256"] == sha(classification / "membership_frozen.csv")
    for name, digest in {**membership_audit["source_hashes"], **membership_audit["code_hashes"]}.items():
        assert sha(classification / name) == digest, f"Classification evidence changed: {name}"

    files = {}
    for path in ROOT.rglob("*"):
        if path.is_file() and not any(part in ("__pycache__", "runtime", ".git") for part in path.relative_to(ROOT).parts):
            files[str(path.relative_to(ROOT))] = path
    for path in classification.rglob("*"):
        relative = path.relative_to(classification)
        if path.is_file() and not any(part in ("__pycache__", "runtime", ".git") for part in relative.parts):
            files[str(Path("classification_evidence") / relative)] = path
    supporting = [
        "evaluate_benchmark.py", "evaluate_exposure_sensitivity.py",
        "TRAINING_EXPOSURE_SENSITIVITY_zh.md", "controller_state.json", "controller_phase_results.json",
        "references/training_exposure_audit_zh.md", "references/model_input_output_capabilities_zh.md",
        "references/model_input_output_capabilities_zh.csv", "references/journal_benchmark_figure_references_zh.md",
        "fastprotflex/full_audit.json",
        "fastprotflex/full_audit.lifetime_v2.json", "fastprotflex/run_local_lifetime_v2.py",
        "fastprotflex/memory_lifetime_v2/candidate/FUNCTION_GDV.r",
        "fastprotflex/memory_lifetime_v2/runs/smoke_01/validation_audit.json",
        "fastprotflex/memory_lifetime_v2/working_set_trim_events.json",
    ]
    for name in supporting:
        path = SOURCE / name
        assert path.is_file(), name
        files[str(Path("completed_inference_evidence") / name)] = path
    # Include all persisted small resource-calibration evidence, not model arrays.
    for path in (SOURCE / "fastprotflex/memory_lifetime_v2").rglob("*.json"):
        relative = path.relative_to(SOURCE)
        files[str(Path("completed_inference_evidence") / relative)] = path
    inventory = [{"snapshot_path": name, "source_path": str(path), "sha256": sha(path), "bytes": path.stat().st_size}
        for name, path in sorted(files.items())]
    if DESTINATION.exists():
        manifest = json.loads((DESTINATION / "SNAPSHOT_MANIFEST.json").read_text(encoding="utf-8"))
        assert manifest["files"] == inventory, "Snapshot exists with different contents; do not overwrite."
    else:
        DESTINATION.mkdir(parents=True)
        for row in inventory:
            destination = DESTINATION / row["snapshot_path"]
            destination.parent.mkdir(parents=True, exist_ok=True)
            shutil.copyfile(row["source_path"], destination)
            assert sha(destination) == row["sha256"]
        manifest = {"created_utc": datetime.now(timezone.utc).isoformat(),
            "scope": "Fourteen completed method/configuration entries, frozen biological proteins/protein-complex membership, common scores, verified exports and source evidence. Heavy full-cohort models remain deferred.",
            "records": 3063, "residues": 961321, "group_counts": scoring["group_counts"],
            "membership_sha256": scoring["membership_sha256"], "qa_status": qa["status"],
            "files": inventory, "raw_model_outputs": "Original output arrays and weights remain at the paths in the frozen source audits; this snapshot copies small provenance and completed scoring tables, not all model environments/weights.",
            "historical_controller_status": "Historical controller visual_review_pending is preserved verbatim. The new two-group QA and visual hashes establish this later review; no historical status was rewritten."}
        (DESTINATION / "SNAPSHOT_MANIFEST.json").write_text(json.dumps(manifest, indent=2, ensure_ascii=False), encoding="utf-8")
    assert all(sha(row["source_path"]) == row["sha256"] == sha(DESTINATION / row["snapshot_path"]) for row in inventory)
    print(json.dumps({"status": "frozen_and_verified", "snapshot": str(DESTINATION), "files": len(inventory),
        "megabytes": sum(row["bytes"] for row in inventory) / 1048576}, indent=2), flush=True)


if __name__ == "__main__":
    main()
