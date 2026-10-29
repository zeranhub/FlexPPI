"""Independent read-only batch coverage review plus ten preselected CPU repeats.

Reads structure identities and existing prediction files, never SPR measurements.
No existing pipeline script, original prediction, or input structure is modified.
"""
from collections import Counter
import csv
import ctypes
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import threading
import time

REPRO = Path(__file__).resolve().parent
RUNTIME = REPRO.parent
ROOT = RUNTIME.parent
sys.path.insert(0, str(RUNTIME))
from local_cpu_inference import import_runtime, load_audited_model, sha256, write_json


class MemoryCounters(ctypes.Structure):
    _fields_ = [("cb", ctypes.c_ulong), ("PageFaultCount", ctypes.c_ulong)] + [
        (name, ctypes.c_size_t) for name in ["PeakWorkingSetSize", "WorkingSetSize",
            "QuotaPeakPagedPoolUsage", "QuotaPagedPoolUsage", "QuotaPeakNonPagedPoolUsage",
            "QuotaNonPagedPoolUsage", "PagefileUsage", "PeakPagefileUsage", "PrivateUsage"]]


def memory_bytes():
    fn = ctypes.windll.psapi.GetProcessMemoryInfo
    fn.argtypes = [ctypes.c_void_p, ctypes.POINTER(MemoryCounters), ctypes.c_ulong]
    fn.restype = ctypes.c_int
    process = ctypes.windll.kernel32.GetCurrentProcess
    process.restype = ctypes.c_void_p
    data = MemoryCounters()
    data.cb = ctypes.sizeof(data)
    if not fn(process(), ctypes.byref(data), data.cb):
        raise ctypes.WinError()
    return {"working_set": int(data.WorkingSetSize), "process_peak_working_set": int(data.PeakWorkingSetSize),
        "private_usage": int(data.PrivateUsage)}


class MemorySampler:
    def __enter__(self):
        self.stop = threading.Event()
        self.samples = [memory_bytes()]
        def sample():
            while not self.stop.wait(0.01):
                self.samples.append(memory_bytes())
        self.worker = threading.Thread(target=sample, daemon=True)
        self.worker.start()
        return self

    def __exit__(self, *args):
        self.stop.set()
        self.worker.join()
        self.samples.append(memory_bytes())

    def report(self):
        return {"sampling_interval_seconds": 0.01, "sample_count": len(self.samples),
            "working_set_before_bytes": self.samples[0]["working_set"],
            "sampled_max_working_set_bytes": max(row["working_set"] for row in self.samples),
            "sampled_working_set_increase_bytes": max(row["working_set"] for row in self.samples) - self.samples[0]["working_set"],
            "process_lifetime_peak_working_set_after_bytes": self.samples[-1]["process_peak_working_set"],
            "sampled_max_private_usage_bytes": max(row["private_usage"] for row in self.samples),
            "interpretation": "Working-set samples cover this forward at 10ms; transient shorter peaks may be missed. The process-lifetime peak includes earlier imports, audit and forwards, not a reset per-forward peak."}


def main():
    started = time.perf_counter()
    runtime = import_runtime(threads=2)
    np, torch, parse_pdb = runtime["np"], runtime["torch"], runtime["parse_pdb"]
    batch_dir = ROOT / "predictions/fixed_backbone"
    batch_audit_path = batch_dir / "run_audit.json"
    manifest_path = ROOT / "structures/fixed_backbone_manifest.json"
    batch = json.loads(batch_audit_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    original_batch_hash = sha256(batch_audit_path)
    checkpoint = ROOT / "model/epoch_040.pt"
    expected_sha = "b3eed060919df15f075329bc892f959d9355232be5feb561417edd49ea82b39f"
    assert sha256(checkpoint) == expected_sha == batch["checkpoint_sha256"] == batch["expected_checkpoint_sha256"]
    assert batch["status"] == "passed" and batch["device"] == "cpu"
    assert batch["threads"] == 2 and not batch["remote_compute"] and not batch["esm_used"]
    assert not batch["prepare_only"]
    assert batch["sequence_backend"] == "learned"
    assert sha256(manifest_path) == batch["manifest_sha256"]
    assert len(batch["candidates"]) == len(manifest["candidates"]) == 76
    ids = [row["candidate_id"] for row in manifest["candidates"]]
    actual_ids = [row["candidate_id"] for row in batch["candidates"]]
    assert len(set(ids)) == 76 and set(ids) == set(actual_ids)
    assert {p.name for p in batch_dir.iterdir() if p.is_dir()} == set(ids)
    source_checks = []
    for name, expected in batch["source_sha256"].items():
        p = Path(batch["source_root"]) / name
        assert sha256(p) == expected, str(p)
        source_checks.append({"path": str(p), "sha256": expected, "passed": True})
    refs = {row["candidate_id"]: row for row in batch["candidates"]}
    template_rows = {row["template_pdb"]: row for row in manifest["templates"]}
    originals = {k: parse_pdb(row["raw_pdb"]) for k, row in template_rows.items()}
    template_counts = Counter(row["template_pdb"] for row in manifest["candidates"])
    assert template_counts == {"2BNQ": 14, "2BNR": 14, "1AO7": 16, "1QSE": 16, "1QSF": 16}
    assert sum(row["is_panel_parent"] for row in manifest["candidates"]) == 5
    assert sum(row["primary"] for row in manifest["candidates"]) == 30
    assert sum(row["sensitivity"] for row in manifest["candidates"]) == 46
    for template in template_rows:
        in_template = [row for row in manifest["candidates"] if row["template_pdb"] == template]
        assert sum(row["is_panel_parent"] for row in in_template) == 1
        panel = in_template[0]["panel_id"]
        panel_ids = {row["biological_candidate_id"] for row in manifest["candidates"] if row["panel_id"] == panel}
        assert {row["biological_candidate_id"] for row in in_template} == panel_ids
    original_output_hashes = {}
    coverage = []
    total_residues = 0
    for row in manifest["candidates"]:
        candidate_id = row["candidate_id"]
        ref = refs[candidate_id]
        source_input = Path(row["prepared_input_pdb"])
        actual_input = Path(ref["input_pdb"])
        source_sha, actual_sha = sha256(source_input), sha256(actual_input)
        assert source_sha == actual_sha == row["prepared_input_sha256"] == ref["input_sha256"]
        assert actual_sha == ref["preparation"]["prepared_pdb_sha256"]
        source = parse_pdb(str(source_input))
        actual = parse_pdb(str(actual_input))
        template = originals[row["template_pdb"]]
        for name in ["aatype", "backbone", "backbone_mask", "chain_index", "residue_index"]:
            assert np.array_equal(getattr(source, name), getattr(actual, name)), (candidate_id, name)
        for name in ["backbone", "backbone_mask", "chain_index", "residue_index"]:
            assert np.array_equal(getattr(actual, name), getattr(template, name)), (candidate_id, name)
        assert source.sequence == actual.sequence
        assert actual.chain_sequence(row["peptide_chain"]) == row["peptide_sequence"]
        assert actual.chain_ids == template.chain_ids
        assert actual.author_residue_ids == template.author_residue_ids
        for chain in template.unique_chains():
            if chain != row["peptide_chain"]:
                assert actual.chain_sequence(chain) == template.chain_sequence(chain)
        folder = batch_dir / candidate_id
        npz_path, csv_path = folder / "prediction.npz", folder / "prediction.csv"
        for path in [source_input, actual_input, npz_path, csv_path, folder / "candidate_audit.json", folder / "input_preparation_audit.json"]:
            original_output_hashes[str(path)] = sha256(path)
        with np.load(npz_path, allow_pickle=False) as data:
            n = actual.num_residues
            assert data["mean"].shape == data["log_variance"].shape == (n, 2)
            assert data["metrics"].tolist() == ["rmsf", "chi_entropy"]
            assert str(data["checkpoint_sha256"]) == expected_sha
            assert np.isfinite(data["mean"]).all() and np.isfinite(data["log_variance"]).all()
            assert (data["mean"] >= 0).all()
            assert (data["log_variance"] >= -8).all() and (data["log_variance"] <= 6).all()
            assert np.array_equal(data["chain"], np.asarray(actual.chain_ids))
            assert np.array_equal(data["author_residue_id"], np.asarray(actual.author_residue_ids))
            assert np.array_equal(data["amino_acid"], np.asarray(list(actual.sequence)))
            assert np.array_equal(data["geometry_valid"], actual.frame_mask)
            assert data["target_mask"].all(), "Expected complete-complex output"
            expected_hchi = np.asarray([bool(runtime["CHI_ATOMS"].get(name)) for name in actual.residue_names])
            assert np.array_equal(data["hchi_chemically_defined"], expected_hchi)
            assert actual.frame_mask[actual.indices_for_chain(row["peptide_chain"])].all()
            assert int(data["geometry_valid"].sum()) == ref["geometry_valid_count"]
            assert int((data["geometry_valid"] & expected_hchi).sum()) == ref["hchi_chemically_defined_count"]
            with csv_path.open(encoding="utf-8", newline="") as handle:
                table = list(csv.DictReader(handle))
            assert len(table) == n
            for i, csv_row in enumerate(table):
                assert csv_row["candidate_id"] == candidate_id
                assert int(csv_row["residue_index"]) == i
                assert csv_row["chain"] == actual.chain_ids[i]
                assert csv_row["author_residue_id"] == actual.author_residue_ids[i]
                assert csv_row["amino_acid"] == actual.sequence[i]
                for column, key, j in [("rmsf", "mean", 0), ("chi_entropy", "mean", 1),
                    ("rmsf_log_variance", "log_variance", 0), ("chi_entropy_log_variance", "log_variance", 1)]:
                    assert float(csv_row[column]) == float(data[key][i, j])
                for column, key in [("geometry_valid", "geometry_valid"), ("target", "target_mask"),
                    ("hchi_chemically_defined", "hchi_chemically_defined")]:
                    assert (csv_row[column] == "True") == bool(data[key][i])
            coverage.append({"candidate_id": candidate_id, "template_pdb": row["template_pdb"],
                "is_panel_parent": row["is_panel_parent"], "residues": n,
                "geometry_valid_count": int(data["geometry_valid"].sum()),
                "hchi_chemically_defined_count": int(expected_hchi.sum()),
                "input_sha256": actual_sha, "source_fixed_and_actual_input_identical": True,
                "AA_NCAO_and_masks_verified": True, "all_outputs_finite": True,
                "npz_sha256": sha256(npz_path), "csv_sha256": sha256(csv_path), "status": "passed"})
            total_residues += n
    print(json.dumps({"coverage": "passed", "inputs": len(coverage), "residue_rows": total_residues}), flush=True)
    selected = []
    for template in ["2BNQ", "2BNR", "1AO7", "1QSE", "1QSF"]:
        rows = [row for row in manifest["candidates"] if row["template_pdb"] == template]
        parent = next(row for row in rows if row["is_panel_parent"])
        single = next(row for row in rows if sum(a != b for a, b in zip(row["peptide_sequence"], row["panel_reference_peptide"])) == 1)
        selected.extend([{**parent, "repeat_selection": "fixed panel parent"},
                         {**single, "repeat_selection": "first single substitution in original manifest order, relative to panel parent"}])
    write_json(REPRO / "preselected_repeats.json", [{k: row[k] for k in ["candidate_id", "template_pdb", "repeat_selection"]} for row in selected])
    load_start = time.perf_counter()
    model, config = load_audited_model(checkpoint, runtime)
    model_load_seconds = time.perf_counter() - load_start
    assert all(parameter.device.type == "cpu" for parameter in model.parameters())
    repeats = []
    for row in selected:
        candidate_id = row["candidate_id"]
        source_input = Path(row["prepared_input_pdb"])
        ex = runtime["_structure_example"](str(source_input), None, 300.0, 200.0, 2.0)
        collated = runtime["collate_states"]([ex])
        assert all(value.device.type == "cpu" for value in collated.values() if torch.is_tensor(value))
        with MemorySampler() as mem:
            forward_start = time.perf_counter()
            with torch.inference_mode():
                prediction = model.forward_state(collated)
            forward_seconds = time.perf_counter() - forward_start
        mean = prediction["mean"][0].numpy()
        log_variance = prediction["log_variance"][0].numpy()
        mask = prediction["mask"][0].numpy()
        with np.load(batch_dir / candidate_id / "prediction.npz", allow_pickle=False) as ref:
            mean_err = np.abs(mean - ref["mean"])
            var_err = np.abs(log_variance - ref["log_variance"])
            assert np.array_equal(mask, ref["geometry_valid"])
            assert float(mean_err.max()) < 2e-5
            assert float(var_err.max()) < 2e-5
            mean_bit_identical = np.array_equal(mean, ref["mean"])
            logvar_bit_identical = np.array_equal(log_variance, ref["log_variance"])
        output = REPRO / (candidate_id + "__repeat.npz")
        np.savez_compressed(output, mean=mean, log_variance=log_variance, geometry_valid=mask,
            checkpoint_sha256=np.asarray(expected_sha), input_sha256=np.asarray(sha256(source_input)))
        result = {"candidate_id": candidate_id, "template_pdb": row["template_pdb"],
            "selection": row["repeat_selection"], "source_input_sha256": sha256(source_input),
            "device": "cpu", "threads": torch.get_num_threads(), "forward_seconds": forward_seconds,
            "max_absolute_mean_error": float(mean_err.max()), "mean_absolute_mean_error": float(mean_err.mean()),
            "max_absolute_log_variance_error": float(var_err.max()),
            "mean_bit_identical": mean_bit_identical, "log_variance_bit_identical": logvar_bit_identical,
            "mask_identical": True, "all_outputs_finite": bool(np.isfinite(mean).all() and np.isfinite(log_variance).all()),
            "memory": mem.report(), "repro_output": str(output), "repro_output_sha256": sha256(output), "status": "passed"}
        repeats.append(result)
        write_json(REPRO / "repeat_results.json", repeats)
        print(json.dumps({"repeat": candidate_id, "seconds": forward_seconds,
            "max_error": result["max_absolute_mean_error"]}), flush=True)
    for path, digest in original_output_hashes.items():
        assert sha256(path) == digest, "Original input/output changed during review: " + path
    assert sha256(batch_audit_path) == original_batch_hash
    assert sha256(manifest_path) == batch["manifest_sha256"]
    result = {"status": "passed", "reviewer": "/root/case_figure_design_20261021",
        "reviewed_at_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Independent inference coverage/provenance review and ten fixed CPU repeat forwards; no affinity measurements or ranking metrics read.",
        "read_SPR_measurements": False, "remote_compute": False,
        "existing_scripts_modified": False, "batch_results_or_structures_modified": False,
        "review_script": str(Path(__file__)), "review_script_sha256": sha256(Path(__file__)),
        "manifest_path": str(manifest_path), "manifest_sha256": sha256(manifest_path),
        "batch_audit_path": str(batch_audit_path), "batch_audit_sha256": original_batch_hash,
        "checkpoint": str(checkpoint), "checkpoint_sha256": expected_sha,
        "model_identity": batch["model_identity"], "model_config": config,
        "device": "cpu", "torch": torch.__version__, "threads": torch.get_num_threads(),
        "source_files_verified": source_checks,
        "coverage_summary": {"inference_candidates": 76, "unique_biological_candidates": 30,
            "per_template": dict(template_counts), "panel_parents": 5, "primary_inputs": 30,
            "sensitivity_inputs": 46, "total_residue_rows_verified": total_residues,
            "predicted_mean_values_verified": total_residues * 2,
            "log_variance_values_verified": total_residues * 2,
            "all_76_source_fixed_PDB_SHA_equal_actual_prediction_input_SHA": True,
            "all_AA_NCAO_and_masks_verified": True, "all_NPZ_CSV_values_identical": True,
            "all_outputs_finite": True, "all_full_complex_targets_retained": True},
        "coverage": coverage, "repeat_selection_rule": "Per template, one fixed panel parent then first panel-relative single substitution in original manifest order; fixed before repeat outputs.",
        "repeat_count": len(repeats), "repeats": repeats,
        "repeat_model_load_seconds": model_load_seconds,
        "repeat_forward_seconds_total": sum(row["forward_seconds"] for row in repeats),
        "repeat_max_absolute_mean_error": max(row["max_absolute_mean_error"] for row in repeats),
        "repeat_max_absolute_log_variance_error": max(row["max_absolute_log_variance_error"] for row in repeats),
        "all_repeats_bit_identical": all(row["mean_bit_identical"] and row["log_variance_bit_identical"] for row in repeats),
        "all_original_file_hashes_rechecked_unchanged_after_review": True,
        "elapsed_seconds": time.perf_counter() - started,
        "interpretation": "This verifies reproducibility and provenance of the historical strict_full 2431 CPU outputs, not their affinity-ranking accuracy or equivalence to the manuscript2430 model."}
    write_json(RUNTIME / "inference_final_review.json", result)
    print(json.dumps({"status": "passed", "repeat_n": len(repeats),
        "max_mean_error": result["repeat_max_absolute_mean_error"],
        "all_bit_identical": result["all_repeats_bit_identical"],
        "report": str(RUNTIME / "inference_final_review.json")}), flush=True)


if __name__ == "__main__":
    main()
