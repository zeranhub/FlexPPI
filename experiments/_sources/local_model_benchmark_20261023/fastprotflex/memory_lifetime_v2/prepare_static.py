"""Create a lifecycle-only candidate; Python-only, never executes R."""
from datetime import datetime, timezone
import difflib
import hashlib
import json
from pathlib import Path
import re

HERE = Path(__file__).resolve().parent
OFFICIAL = Path("E:/Project/protein_flexibility_models/FastProtFlex")
BASELINE = HERE.parent


def sha(data):
    return hashlib.sha256(data).hexdigest()


def numeric_program_lines(text):
    # The only allowed removed/added executable lines are rm(...) and gc().
    return [line for line in text.splitlines() if not re.fullmatch(r"\s*(?:rm\([^;]*\)|gc\(\))\s*", line)]


def main():
    original = (OFFICIAL / "FUNCTION_GDV.r").read_bytes()
    original_text = original.decode("utf-8-sig")
    lines = original_text.splitlines()
    # Original one-based line numbers were audited before making this copy.
    additions = {
        14: ["  rm(DD)", "  gc()"],
        16: ["  rm(adj)", "  gc()"],
        19: ["  rm(g)", "  gc()"],
        22: ["  rm(DDsave)", "  gc()"],
        23: ["  rm(DD)", "  gc()"],
    }
    assert lines[13].strip() == "adj<-as.matrix(DD)"
    assert lines[15].strip() == 'g<-graph_from_adjacency_matrix(adj,mode="undirected",weighted=NULL)'
    assert lines[17].strip() == "rm(adj,DD)"
    assert lines[18].strip() == "go<-count_orbits_per_node(g,max_graphlet_size=4) # orbit"
    assert lines[21].lstrip().startswith("DD[DDsave<=2.1]<-1")
    assert lines[22].strip() == "adj<-as.matrix(DD)"
    candidate_lines = []
    for line_number, line in enumerate(lines, 1):
        if line_number != 18:
            candidate_lines.append(line)
        candidate_lines.extend(additions.get(line_number, []))
    original_newline = "\r\n" if "\r\n" in original_text else "\n"
    candidate = original_newline.join(candidate_lines) + (original_newline if original_text.endswith("\n") else "")
    assert numeric_program_lines(original_text) == numeric_program_lines(candidate)
    diff = list(difflib.unified_diff(original_text.splitlines(keepends=True), candidate.splitlines(keepends=True),
                                   fromfile="official/FUNCTION_GDV.r", tofile="memory_lifetime_v2/FUNCTION_GDV.r"))
    changed = [line for line in diff if line.startswith(("+", "-")) and not line.startswith(("+++", "---"))]
    assert len(changed) == 11
    assert all(re.fullmatch(r"[+-]\s*(?:rm\([^;]*\)|gc\(\))\s*", line) for line in changed)
    for directory in ("original", "candidate", "wrapper"):
        (HERE / directory).mkdir(parents=True, exist_ok=True)
    (HERE / "original/FUNCTION_GDV.r").write_bytes(original)
    (HERE / "candidate/FUNCTION_GDV.r").write_text(candidate, encoding="utf-8", newline="")
    (HERE / "FUNCTION_GDV_lifecycle_only.diff").write_text("".join(diff), encoding="utf-8", newline="")
    for filename, source in (("predict.r", OFFICIAL / "predict.r"),
                             ("capture_full_precision.R", BASELINE / "capture_full_precision.R")):
        (HERE / "wrapper" / filename).write_bytes(source.read_bytes())
    cases = []
    for sample_id, atoms in (("2o70_F", 1385), ("A007V", 3072), ("itrace_1qrn", 6620)):
        metadata_path = BASELINE / "records" / sample_id / "validation.json"
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        rds = metadata_path.parent / "atom_predictions_exact.rds"
        assert metadata["status"] == "complete" and metadata["input_protein_heavy_atoms"] == atoms
        assert sha(rds.read_bytes()) == metadata["rds_exact_output_sha256"]
        cases.append({"sample_id": sample_id, "n_atoms": atoms, "baseline_validation": str(metadata_path),
                      "baseline_validation_sha256": sha(metadata_path.read_bytes()),
                      "baseline_RDS": str(rds), "baseline_RDS_sha256": sha(rds.read_bytes()),
                      "input_sha256": metadata["input_sha256"],
                      "baseline_R_process_peak_MiB": metadata["peak_process_tree_rss_mib"],
                      "baseline_wall_seconds": metadata["wall_seconds"]})
    audit = {
        "status": "prepared_static_only_R_not_executed", "created_utc": datetime.now(timezone.utc).isoformat(),
        "official_function_path": str(OFFICIAL / "FUNCTION_GDV.r"), "official_function_sha256": sha(original),
        "candidate_function_sha256": sha(candidate.encode("utf-8")),
        "official_predict_sha256": sha((OFFICIAL / "predict.r").read_bytes()),
        "capture_wrapper_sha256": sha((BASELINE / "capture_full_precision.R").read_bytes()),
        "diff_sha256": sha((HERE / "FUNCTION_GDV_lifecycle_only.diff").read_bytes()),
        "non_cleanup_program_text_identical_in_order": True,
        "all_diff_executable_lines_are_only_rm_or_gc": True,
        "removed_cleanup_original_line": 18,
        "insert_cleanup_after_original_lines": additions,
        "numerical_changes": "none; thresholds, functions, smoothing order, normalization and all 15 coefficients unchanged",
        "wrapper_source_resolution": "Byte-identical predict.r and capture_full_precision.R are copied into isolated per-case working directories. Their unchanged relative source calls resolve only to the candidate file in that same directory.",
        "cases": cases,
        "acceptance": {
            "same_pdb_bytes": True, "same_R_libraries_and_single_CPU": True,
            "all_atom_table_identical_num_eq_FALSE": True,
            "all_atom_identity_columns_identical_num_eq_FALSE": True,
            "all_atom_xyz_identical_num_eq_FALSE": True,
            "all_atom_prediction_vector_identical_num_eq_FALSE": True,
            "maximum_absolute_prediction_difference_must_equal": 0,
            "finite_predictions_and_same_atom_count": True,
            "never_use_rounded_Bfactor_PDB_as_reference": True,
            "no_candidate_adoption_on_any_validation_failure": True,
        },
        "limits": "This preparation does not prove any measured RSS reduction. Original function, running wrapper, prediction outputs, and original memory gate are untouched. Run only after explicit root instruction and when no other R process exists.",
    }
    (HERE / "lifecycle_provenance.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(json.dumps({key: audit[key] for key in ("status", "official_function_sha256", "candidate_function_sha256", "non_cleanup_program_text_identical_in_order", "all_diff_executable_lines_are_only_rm_or_gc")}, indent=2))


if __name__ == "__main__":
    main()
