"""Isolated resume entry point; default inspect only, never invents calibration.

Existing official_v1 predictions are hash/identity verified and never recomputed.
Missing records may use the lifecycle-only function only after an actual exact
three-source validation. The calibrated RAM gate is never extrapolated beyond
the largest measured input; larger missing records stop for resource review.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import time

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
COHORT = ROOT / "cohort"
CANDIDATE_ROOT = HERE / "memory_lifetime_v2"
UPSTREAM = Path("E:/Project/protein_flexibility_models/FastProtFlex")
RSCRIPT = Path("C:/Program Files/R/R-4.6.0/bin/Rscript.exe")
R_LIBRARIES = ["E:/Project/protein_flexibility_models/R_library", "E:/Project/protein_flexibility_models/r_libs"]
GATE_VERSION = "lifetime_v2_three_source_rss115_reserve256_commit1024_v2"
EXTENDED_GATE_VERSION = "lifetime_v2_three_source_plus8058_rss115_reserve256_commit1024_v2"
RESOURCE_POLICY_VERSION = "lifetime_v2_rss115_reserve256_commit1024_v2"
POLICY_CHANGE_REASON = "Explicit root authorization after actual three-source bitwise equality and measured upper envelope: replace discretionary512MiB physical reserve with256MiB while retaining15percent, and add independent system-commit headroom of predicted peak plus1024MiB. No further automatic reserve reduction. Historical512MiB audits remain unchanged."
TRUTH_COLUMNS = ["sample_id", "residue_index", "amino_acid", "chain", "author_residue_id", "target_rmsf", "valid_rmsf"]
AA = dict(zip("ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(), "ARNDCQEGHILKMFPSTWYV"))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024**2), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path, value):
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False), encoding="utf-8")
    temporary.replace(path)


def json_read(path):
    return json.loads(Path(path).read_text(encoding="utf-8-sig"))


def csv_read(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def non_cleanup_lines(path):
    return [line for line in Path(path).read_bytes().decode("utf-8-sig").splitlines()
            if not re.fullmatch(r"\s*(?:rm\([^;]*\)|gc\(\))\s*", line)]


def inspect_sources():
    provenance = json_read(CANDIDATE_ROOT / "lifecycle_provenance.json")
    assert sha256(UPSTREAM / "FUNCTION_GDV.r") == provenance["official_function_sha256"]
    assert sha256(CANDIDATE_ROOT / "candidate/FUNCTION_GDV.r") == provenance["candidate_function_sha256"]
    assert non_cleanup_lines(UPSTREAM / "FUNCTION_GDV.r") == non_cleanup_lines(CANDIDATE_ROOT / "candidate/FUNCTION_GDV.r")
    assert sha256(UPSTREAM / "predict.r") == provenance["official_predict_sha256"]
    assert sha256(HERE / "capture_full_precision.R") == provenance["capture_wrapper_sha256"]
    schedule = csv_read(HERE / "full_schedule.csv")
    records = csv_read(COHORT / "records.csv")
    lookup = {record["sample_id"]: record for record in records}
    assert len(schedule) == len(records) == len(lookup) == 3063
    assert {row["sample_id"] for row in schedule} == set(lookup)
    assert [int(row["schedule_order"]) for row in schedule] == list(range(1, 3064))
    assert [int(row["protein_heavy_atoms"]) for row in schedule] == sorted(int(row["protein_heavy_atoms"]) for row in schedule)
    for row in schedule:
        for key in ("n_chains", "n_residues", "main_group_candidate"):
            assert str(row[key]) == str(lookup[row["sample_id"]][key])
    return provenance, schedule, lookup


def fit_resource_envelope(points):
    """Monotone N^2 envelope over actual measurements, before RAM reserves."""
    xs = [float(point["n_atoms"])**2 for point in points]
    ys = [float(point["measured_peak_MiB"]) for point in points]
    assert all(value > 0 and value < float("inf") for value in ys)
    mean_x, mean_y = sum(xs)/len(xs), sum(ys)/len(ys)
    slope = max(0.0, sum((x-mean_x)*(y-mean_y) for x, y in zip(xs, ys)) / sum((x-mean_x)**2 for x in xs))
    intercept = max(0.0, max(y - slope*x for x, y in zip(xs, ys)))
    assert all(intercept + slope*x + 1e-8 >= y for x, y in zip(xs, ys))
    return {"measured_points": points, "quadratic_intercept_MiB": intercept,
            "quadratic_slope_MiB_per_atom_squared": slope,
            "maximum_measured_atoms": max(point["n_atoms"] for point in points),
            "minimum_measured_atoms": min(point["n_atoms"] for point in points),
            "minimum_peak_floor_MiB": min(ys), "safety_multiplier": 1.15,
            "additional_available_RAM_reserve_MiB": 256.0,
            "additional_available_commit_reserve_MiB":1024.0,
            "prior_additional_available_RAM_reserve_MiB":512.0,
            "resource_policy_version":RESOURCE_POLICY_VERSION, "policy_change_reason":POLICY_CHANGE_REASON}


def resource_requirements(calibration, n_atoms):
    estimate = max(calibration["minimum_peak_floor_MiB"], calibration["quadratic_intercept_MiB"] + calibration["quadratic_slope_MiB_per_atom_squared"]*n_atoms**2)
    return {"estimated_R_process_tree_RSS_MiB":estimate, "required_available_RAM_MiB":1.15*estimate+256.0,
            "previous_required_available_RAM_MiB":1.15*estimate+512.0, "required_available_commit_MiB":estimate+1024.0,
            "safety_multiplier":1.15, "additional_available_RAM_reserve_MiB":256.0,
            "additional_available_commit_reserve_MiB":1024.0, "resource_policy_version":RESOURCE_POLICY_VERSION,
            "policy_change_reason":POLICY_CHANGE_REASON}


def windows_commit_status():
    """System commit counters, independent from resident physical RAM."""
    try:
        import ctypes
        class PerformanceInformation(ctypes.Structure):
            _fields_ = [("cb",ctypes.c_ulong)] + [(name,ctypes.c_size_t) for name in
                ("CommitTotal","CommitLimit","CommitPeak","PhysicalTotal","PhysicalAvailable","SystemCache","KernelTotal","KernelPaged","KernelNonpaged","PageSize")] + [(name,ctypes.c_ulong) for name in ("HandleCount","ProcessCount","ThreadCount")]
        info = PerformanceInformation()
        info.cb = ctypes.sizeof(info)
        assert ctypes.WinDLL("psapi").GetPerformanceInfo(ctypes.byref(info),ctypes.sizeof(info))
        return {"available_system_commit_MiB":(info.CommitLimit-info.CommitTotal)*info.PageSize/1024**2,
                "system_commit_limit_MiB":info.CommitLimit*info.PageSize/1024**2, "system_commit_counter":"Windows_GetPerformanceInfo", "system_commit_counter_error":None}
    except Exception as error:
        return {"available_system_commit_MiB":None, "system_commit_limit_MiB":None,
                "system_commit_counter":"unknown", "system_commit_counter_error":repr(error)}


def bind_calibration_hash(calibration):
    calibration.pop("calibration_sha256", None)
    calibration["calibration_sha256"] = hashlib.sha256(json.dumps(calibration, sort_keys=True).encode()).hexdigest()
    return calibration


def load_actual_calibration(path, provenance):
    """Require real outputs and measurements; fit a monotone upper envelope."""
    path = Path(path).resolve(strict=True)
    audit = json_read(path)
    assert audit["status"] == "three_sources_exact_all_atom_match"
    assert audit["cases_requested"] == audit["cases_attempted"] == 3
    assert audit["provenance"]["candidate_function_sha256"] == provenance["candidate_function_sha256"]
    assert audit["provenance"]["official_function_sha256"] == provenance["official_function_sha256"]
    assert audit["validation_script_sha256"] == sha256(CANDIDATE_ROOT / "validate_candidate.py")
    assert audit["R_validation_wrapper_sha256"] == sha256(CANDIDATE_ROOT / "capture_and_validate.R")
    expected = {case["sample_id"]: case for case in provenance["cases"]}
    assert {row["sample_id"] for row in audit["results"]} == set(expected)
    points = []
    required_comparisons = ("passed", "same_atom_count", "same_columns", "all_atom_identities_identical",
                           "all_atom_coordinates_identical", "all_atom_table_identical", "all_predictions_finite",
                           "all_atom_prediction_vector_identical_num_eq_FALSE")
    for result in audit["results"]:
        case = expected[result["sample_id"]]
        assert result["status"] == "complete_exact_all_atom_match"
        assert result["return_code"] == 0 and result["timed_out"] is False
        assert result["R_and_package_versions_match_baseline"] is True
        assert result["library_thread_caps"] == 1 and result["sampling_interval_seconds"] <= .1
        assert result["candidate_function_sha256"] == provenance["candidate_function_sha256"]
        assert result["input_sha256"] == case["input_sha256"]
        comparison = result["all_atom_RDS_comparison"]
        assert all(comparison[key] is True for key in required_comparisons)
        assert comparison["maximum_absolute_prediction_difference"] == 0
        assert comparison["reference_atom_count"] == comparison["candidate_atom_count"] == case["n_atoms"]
        assert result["reference_RDS_sha256"] == case["baseline_RDS_sha256"]
        assert sha256(case["baseline_RDS"]) == case["baseline_RDS_sha256"]
        actual_directory = path.parent / case["sample_id"]
        assert json_read(actual_directory / "validation.json") == result
        assert sha256(actual_directory / "atom_predictions_exact.rds") == result["candidate_RDS_sha256"]
        assert sha256(actual_directory / "FUNCTION_GDV.r") == provenance["candidate_function_sha256"]
        assert sha256(actual_directory / "predict.r") == provenance["official_predict_sha256"]
        assert sha256(actual_directory / "capture_full_precision.R") == provenance["capture_wrapper_sha256"]
        assert sha256(actual_directory / "capture_and_validate.R") == audit["R_validation_wrapper_sha256"]
        assert sha256(actual_directory / "protein.pdb") == case["input_sha256"]
        assert json_read(actual_directory / "comparison.json") == comparison
        peak_mib = float(result["peak_R_process_tree_RSS_bytes"]) / 1024**2
        assert peak_mib > 0 and abs(peak_mib - float(result["peak_R_process_tree_RSS_MiB"])) < 1e-8
        assert result["observed_process_tree_CPU_seconds"] > 0
        points.append({"sample_id": case["sample_id"], "n_atoms": case["n_atoms"], "measured_peak_MiB": peak_mib})
    calibration = {
        "version": GATE_VERSION, "validation_audit_path": str(path), "validation_audit_sha256": sha256(path),
        **fit_resource_envelope(points),
        "fitting_policy": "Nonnegative N^2 slope; intercept covers all actual peaks; physical1.15peak+256MiB and system commit peak+1024MiB are separate launch requirements.",
        "unmeasured_extremes_policy": "Never launch a missing record above maximum measured atom count. Stop for root review; this script has no override to extrapolate or lower the gate for an unvalidated8058-atom input.",
    }
    calibration["policy_change_thresholds_at8058"] = resource_requirements(calibration,8058)
    return bind_calibration_hash(calibration)


def load_resource_extension(path, calibration, provenance, records):
    """Only actual8058-resource evidence extends measured range; no numeric claim."""
    path = Path(path).resolve(strict=True)
    audit = json_read(path)
    assert audit["status"] == "extreme_resource_and_integrity_check_complete_no_reference_equivalence"
    assert audit["resource_policy_version"] == RESOURCE_POLICY_VERSION
    assert audit["sample_id"] == "7zh9_A" and audit["n_atoms"] == 8058
    assert audit["original_full_precision_reference_available"] is False
    assert audit["original_equivalence_compared"] is False
    assert Path(audit["numerical_equivalence_audit_path"]).resolve() == Path(calibration["validation_audit_path"]).resolve()
    assert sha256(audit["numerical_equivalence_audit_path"]) == audit["numerical_equivalence_audit_sha256"] == calibration["validation_audit_sha256"]
    assert audit["candidate_function_sha256"] == provenance["candidate_function_sha256"]
    assert audit["official_function_sha256"] == provenance["official_function_sha256"]
    assert audit["official_predict_sha256"] == provenance["official_predict_sha256"]
    assert audit["capture_wrapper_sha256"] == provenance["capture_wrapper_sha256"]
    probe = CANDIDATE_ROOT / "probe_extreme_resources.py"
    assert Path(audit["probe_script_path"]).resolve() == probe.resolve()
    assert sha256(probe) == audit["probe_script_sha256"]
    assert Path(audit["resume_entrypoint_path"]).resolve() == Path(__file__).resolve()
    assert sha256(__file__) == audit["resume_entrypoint_sha256"]
    record = records["7zh9_A"]
    assert audit["input_sha256"] == record["structure_sha256"] == sha256(record["structure_path"])
    assert audit["cohort_manifest_sha256"] == sha256(COHORT / "records.csv")
    assert audit["schedule_sha256"] == sha256(HERE / "full_schedule.csv")
    assert audit["return_code"] == 0 and audit["timed_out"] is False
    assert audit["R_and_package_versions_match_baseline"] is True
    assert audit["library_thread_caps"] == 1 and audit["sampling_interval_seconds"] <= .1
    assert audit["observed_process_tree_CPU_seconds"] > 0
    assert audit["gate"]["status"] == "ready"
    assert audit["gate"]["available_RAM_before_launch_MiB"] >= audit["gate"]["required_available_RAM_MiB"]
    assert audit["gate"]["available_system_commit_MiB"] >= audit["gate"]["required_available_commit_MiB"]
    assert audit["gate"]["resource_policy_version"] == RESOURCE_POLICY_VERSION
    assert audit["gate"]["source_calibration_sha256"] == calibration["calibration_sha256"]
    assert audit["base_three_source_resource_calibration"] == calibration
    for key,value in resource_requirements(calibration,8058).items():
        assert audit["gate"][key] == value
    integrity = audit["integrity"]
    for key in ("all_atom_identities_equal_input_in_order", "all_atom_coordinates_verified", "all_atom_predictions_finite",
                "ca_identities_equal_frozen_truth_in_order", "all_ca_predictions_finite", "no_duplicate_CA_identity"):
        assert integrity[key] is True
    assert integrity["n_atoms"] == 8058 and integrity["n_CA"] == int(record["n_residues"]) == 1024
    manifest = {entry["relative_path"]: entry for entry in audit["output_files"]}
    assert len(manifest) == len(audit["output_files"])
    for relative, entry in manifest.items():
        output = (path.parent / relative).resolve(strict=True)
        assert path.parent in output.parents
        assert output.stat().st_size == entry["size_bytes"] and sha256(output) == entry["sha256"]
    required = ("FUNCTION_GDV.r", "predict.r", "capture_full_precision.R", "protein.pdb", "prediction.npz",
                "atom_predictions_exact.rds", "atom_predictions_full_precision.csv", "proteinRMSF.pdb",
                "R_session_info.txt", "stdout.log", "stderr.log", "integrity.json")
    required_paths = {"7zh9_A/"+name for name in required}
    assert required_paths <= set(manifest) <= required_paths | {"waiting_for_memory.json"}
    actual_files = {item.relative_to(path.parent).as_posix() for item in path.parent.rglob("*") if item.is_file() and item != path}
    assert actual_files == set(manifest)
    for name, digest in (("FUNCTION_GDV.r", provenance["candidate_function_sha256"]), ("predict.r", provenance["official_predict_sha256"]),
                         ("capture_full_precision.R", provenance["capture_wrapper_sha256"]), ("protein.pdb", record["structure_sha256"])):
        assert manifest["7zh9_A/"+name]["sha256"] == digest
    assert json_read(path.parent / "7zh9_A/integrity.json") == integrity
    peak = audit["peak_R_process_tree_RSS_bytes"] / 1024**2
    assert peak > 0 and abs(peak-audit["peak_R_process_tree_RSS_MiB"]) < 1e-8
    envelope = fit_resource_envelope(calibration["measured_points"] + [{"sample_id":"7zh9_A", "n_atoms":8058, "measured_peak_MiB":peak}])
    assert audit["combined_resource_envelope"] == {"version":EXTENDED_GATE_VERSION, **envelope}
    extended = {**calibration, **envelope,
                "version": EXTENDED_GATE_VERSION, "resource_extension_audit_path":str(path), "resource_extension_audit_sha256":sha256(path),
                "resource_extension_probe_path":str(probe.resolve()), "resource_extension_probe_sha256":sha256(probe),
                "base_three_source_resource_calibration":dict(calibration),
                "extreme_probe_claim":"Measured resources and input/output integrity only; original full-precision equivalence was not tested for7zh9_A.",
                "unmeasured_extremes_policy":"Never launch missing input above8058 actual measured atoms. Original three-source numerical-equivalence evidence remains mandatory."}
    extended["policy_change_thresholds_at8058"] = resource_requirements(extended,8058)
    return bind_calibration_hash(extended)


def author_id(number, insertion):
    number = str(number).strip()
    try:
        number = str(int(float(number)))
    except ValueError:
        pass
    return number + str(insertion).strip()


def read_atoms(path, pd):
    atoms = []
    for line in Path(path).read_text(encoding="ascii").splitlines():
        if not line.startswith(("ATOM  ", "HETATM")):
            continue
        residue = line[17:20].strip()
        if residue not in AA or line[76:78].strip().upper() in {"H", "D"}:
            continue
        atoms.append({"atom": line[12:16].strip(), "residue": residue, "chain": line[21:22].strip() or "_",
                      "author_residue_id": author_id(line[22:26], line[26:27]),
                      "x": float(line[30:38]), "y": float(line[38:46]), "z": float(line[46:54])})
    return pd.DataFrame(atoms)


def identities(frame):
    return list(zip(frame.chain.astype(str), frame.author_residue_id.astype(str), frame.amino_acid.astype(str)))


def validate_input(record, truth, pd):
    assert sha256(record["structure_path"]) == record["structure_sha256"]
    atoms = read_atoms(record["structure_path"], pd)
    ca = atoms[atoms.atom == "CA"].copy()
    ca["amino_acid"] = ca.residue.map(AA)
    assert len(ca) == int(record["n_residues"])
    if truth is not None:
        assert len(ca) == len(truth)
        assert identities(ca) == identities(truth)
    assert not ca.duplicated(["chain", "author_residue_id"]).any()
    return atoms


def cached_prediction(record, atoms, provenance, np):
    directory = HERE / "records" / record["sample_id"]
    success = directory / "validation.json"
    if not success.exists():
        return None
    old = json_read(success)
    if old.get("status") != "complete":
        return None
    assert old["input_sha256"] == record["structure_sha256"]
    assert old["ca_identities_equal_frozen_truth_in_original_order"] is True
    assert old["all_atom_identities_and_coordinates_verified"] is True
    output = directory / "prediction.npz"
    assert sha256(output) == old["prediction_sha256"]
    version = old.get("algorithm_lifecycle_version", "official_v1")
    assert version in {"official_v1", "lifetime_v2"}
    if version == "lifetime_v2":
        assert old["candidate_function_sha256"] == provenance["candidate_function_sha256"]
        assert old["official_algorithm_source_sha256"] == provenance["official_function_sha256"]
        assert sha256(old["lifecycle_validation_audit_path"]) == old["lifecycle_validation_audit_sha256"]
        assert sha256(old["execution_entrypoint_path"]) == old["execution_entrypoint_sha256"]
        if "resource_extension_audit_path" in old:
            assert sha256(old["resource_extension_audit_path"]) == old["resource_extension_audit_sha256"]
            assert sha256(old["resource_extension_probe_path"]) == old["resource_extension_probe_sha256"]
    expected_function_sha = provenance["candidate_function_sha256"] if version == "lifetime_v2" else provenance["official_function_sha256"]
    for filename, expected_sha in (("FUNCTION_GDV.r", expected_function_sha),
                                   ("predict.r", provenance["official_predict_sha256"]),
                                   ("capture_full_precision.R", provenance["capture_wrapper_sha256"]),
                                   ("protein.pdb", record["structure_sha256"]),
                                   ("atom_predictions_full_precision.csv", old["full_precision_atoms_sha256"]),
                                   ("atom_predictions_exact.rds", old["rds_exact_output_sha256"]),
                                   ("proteinRMSF.pdb", old["official_pdb_output_sha256"])):
        assert sha256(directory / filename) == expected_sha
    ca = atoms.loc[atoms.atom == "CA"].copy()
    ca["amino_acid"] = ca.residue.map(AA)
    with np.load(output, allow_pickle=False) as saved:
        assert saved["prediction"].shape == (len(ca),) and np.isfinite(saved["prediction"]).all()
        for key in ("chain", "author_residue_id", "amino_acid"):
            assert np.array_equal(saved[key].astype(str), ca[key].astype(str).to_numpy())
        assert saved["residue_index"].shape == saved["valid"].shape == saved["target"].shape == (len(ca),)
        assert saved["residue_index"].dtype.kind in "iu" and saved["valid"].dtype.kind == "b"
    return {**old, "cache_reused_without_forward": True, "cache_algorithm_lifecycle_version": version,
            "validation_file_sha256": sha256(success), "cache_frozen_truth_revalidated_this_invocation":False,
            "final_scorer_frozen_truth_revalidation_required":True,
            "cache_verification_scope":"Input/source/all recorded output SHA plus frozen-PDB-to-NPZ chain/authorID/amino-acid identity and array shapes. Target/valid/residue_index frozen-truth equality is deferred to the independent final scorer."}


def load_missing_truth(records, missing_ids, pd):
    if not missing_ids:
        return {}
    missing_ids = set(missing_ids)
    reference_paths = {records[sample]["reference_table_path"] for sample in missing_ids}
    assert len(reference_paths) == 1
    found = []
    for chunk in pd.read_csv(next(iter(reference_paths)), keep_default_na=False, low_memory=False, chunksize=50000,
                             usecols=TRUTH_COLUMNS, dtype={"chain":str,"author_residue_id":str,"amino_acid":str}):
        selected = chunk.loc[chunk.sample_id.astype(str).isin(missing_ids)]
        if len(selected):
            found.append(selected.copy())
    assert found
    groups = dict(tuple(pd.concat(found,ignore_index=True).groupby("sample_id",sort=False)))
    assert set(groups) == missing_ids
    for sample,truth in groups.items():
        assert len(truth) == int(records[sample]["n_residues"])
    return groups


def competing_jobs(psutil):
    blocked = []
    # Windows venv Python may retain a launcher parent with the exact same
    # script/argument vector. It is this invocation, not a second runner.
    # Exclude only actual Python ancestors whose args match byte-for-byte.
    current = psutil.Process(os.getpid())
    own_pids = {current.pid}
    current_arguments = current.cmdline()[1:]
    for ancestor in current.parents():
        try:
            if ancestor.name().lower() in {"python.exe", "pythonw.exe", "python", "pythonw"} and ancestor.cmdline()[1:] == current_arguments:
                own_pids.add(ancestor.pid)
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    for process in psutil.process_iter(["pid", "name", "cmdline"]):
        try:
            name = (process.info["name"] or "").lower()
            command = " ".join(process.info["cmdline"] or []).replace("\\", "/").lower()
            signatures = ("/fastprotflex/run_local.py", "/fastprotflex/run_local_lifetime_v2.py",
                          "/fastprotflex/memory_lifetime_v2/validate_candidate.py", "/fastprotflex/memory_lifetime_v2/probe_extreme_resources.py")
            if process.pid not in own_pids and (name in {"r.exe", "rscript.exe", "rterm.exe", "rgui.exe"} or ("python" in name and any(signature in command for signature in signatures))):
                blocked.append({"pid": process.pid, "name": name})
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            pass
    return blocked


def memory_gate(n_atoms, sample_id, calibration, maximum_wait, psutil):
    if n_atoms > calibration["maximum_measured_atoms"]:
        return {"status": "blocked_unvalidated_atom_count_requires_root_resource_review", "sample_id": sample_id,
                "input_protein_heavy_atoms": n_atoms, "maximum_measured_atoms": calibration["maximum_measured_atoms"],
                "resource_gate_version": calibration["version"], "reason": "No extrapolation or automatic gate reduction for unvalidated extremes."}
    requirements = resource_requirements(calibration,n_atoms)
    started = time.perf_counter()
    while True:
        available = psutil.virtual_memory().available/1024**2
        commit = windows_commit_status()
        jobs = competing_jobs(psutil)
        ready = available >= requirements["required_available_RAM_MiB"] and commit["available_system_commit_MiB"] is not None and commit["available_system_commit_MiB"] >= requirements["required_available_commit_MiB"] and not jobs
        status = {"status": "ready" if ready else "waiting_for_RAM_commit_or_existing_R", "sample_id": sample_id,
                  "input_protein_heavy_atoms": n_atoms, **requirements, **commit,
                  "estimated_R_process_tree_rss_mib":requirements["estimated_R_process_tree_RSS_MiB"],
                  "required_available_RAM_mib":requirements["required_available_RAM_MiB"], "available_RAM_before_launch_mib":available,
                  "memory_wait_seconds": time.perf_counter()-started, "blocking_jobs": jobs,
                  "resource_gate_version": calibration["version"], "resource_gate_calibration_sha256": calibration["calibration_sha256"]}
        if ready:
            return status
        if status["memory_wait_seconds"] >= maximum_wait:
            status["status"] = "blocked_finite_memory_wait_timeout"
            return status
        write_json(HERE / "waiting_for_memory.lifetime_v2.json", status)
        print(json.dumps(status), flush=True)
        time.sleep(min(30, maximum_wait-status["memory_wait_seconds"]))


def launch_missing(record, truth, atoms, provenance, calibration, args, np, pd, psutil):
    directory = HERE / "records" / record["sample_id"]
    directory.mkdir(parents=True, exist_ok=True)
    gate = memory_gate(len(atoms), record["sample_id"], calibration, args.memory_wait_seconds, psutil)
    if gate["status"] != "ready":
        write_json(HERE / "waiting_for_memory.lifetime_v2.json", gate)
        return {"status": gate["status"], "sample_id": record["sample_id"], "gate": gate}
    for source, filename, expected_sha in ((UPSTREAM / "predict.r", "predict.r", provenance["official_predict_sha256"]),
                             (CANDIDATE_ROOT / "candidate/FUNCTION_GDV.r", "FUNCTION_GDV.r", provenance["candidate_function_sha256"]),
                             (Path(record["structure_path"]), "protein.pdb", record["structure_sha256"]),
                             (HERE / "capture_full_precision.R", "capture_full_precision.R", provenance["capture_wrapper_sha256"])):
        assert sha256(source) == expected_sha
        shutil.copy2(source, directory / filename)
        assert sha256(directory / filename) == expected_sha
    environment = os.environ.copy()
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "RCPP_PARALLEL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        environment[name] = "1"
    environment["R_LIBS_USER"] = ";".join(R_LIBRARIES)
    started_utc = datetime.now(timezone.utc).isoformat()
    started = time.perf_counter()
    peak_rss = peak_threads = 0
    cpu_by_pid = {}
    affinity_cpu = psutil.cpu_count()-1
    with (directory / "stdout.log").open("w", encoding="utf-8") as stdout, (directory / "stderr.log").open("w", encoding="utf-8") as stderr:
        process = subprocess.Popen([str(RSCRIPT), "--vanilla", "capture_full_precision.R"], cwd=directory, env=environment, stdout=stdout, stderr=stderr)
        monitor = psutil.Process(process.pid)
        while process.poll() is None:
            try:
                rss = threads = 0
                for child in [monitor]+monitor.children(recursive=True):
                    try:
                        child.cpu_affinity([affinity_cpu])
                        rss += child.memory_info().rss
                        threads += child.num_threads()
                        usage = child.cpu_times()
                        cpu_by_pid[child.pid] = max(cpu_by_pid.get(child.pid, 0), usage.user+usage.system)
                    except (psutil.NoSuchProcess, psutil.AccessDenied):
                        pass
                peak_rss, peak_threads = max(peak_rss, rss), max(peak_threads, threads)
            except psutil.NoSuchProcess:
                pass
            time.sleep(.1)
        code = process.wait()
    performance = {"sample_id": record["sample_id"], "dataset": record["dataset"], "mode": "full",
                   "input_sha256": record["structure_sha256"], "started_utc": started_utc,
                   "wall_seconds": time.perf_counter()-started, "return_code": code,
                   "peak_process_tree_rss_bytes": peak_rss, "peak_process_tree_rss_mib": peak_rss/1024**2,
                   "observed_max_process_cpu_seconds": max(cpu_by_pid.values(), default=0),
                   "observed_process_tree_CPU_seconds": sum(cpu_by_pid.values()), "observed_process_ids": sorted(cpu_by_pid),
                   "cpu_affinity_logical_cpu": affinity_cpu, "thread_library_caps": 1,
                   "peak_os_threads_including_runtime_threads": peak_threads,
                   "input_protein_heavy_atoms": len(atoms), "n_residues": int(record["n_residues"]), "n_chains": int(record["n_chains"]),
                   "official_branch": "GDVallatonce" if len(atoms)<=10000 else "GDVperpartes_new2", "official_atom_threshold": 10000,
                   "algorithm_lifecycle_version": "lifetime_v2", "official_algorithm_source_sha256": provenance["official_function_sha256"],
                   "candidate_function_sha256": provenance["candidate_function_sha256"],
                   "official_algorithm_source_path": str(UPSTREAM / "FUNCTION_GDV.r"),
                   "candidate_function_path": str(CANDIDATE_ROOT / "candidate/FUNCTION_GDV.r"),
                   "execution_entrypoint_path": str(Path(__file__).resolve()),
                   "execution_entrypoint_sha256": sha256(__file__),
                   "lifecycle_validation_audit_path": calibration["validation_audit_path"],
                   "lifecycle_validation_audit_sha256": calibration["validation_audit_sha256"],
                   "resource_gate_version": calibration["version"], "resource_gate_calibration_sha256": calibration["calibration_sha256"],
                   "resource_policy_version":RESOURCE_POLICY_VERSION,
                   **{key:calibration[key] for key in ("resource_extension_audit_path", "resource_extension_audit_sha256", "resource_extension_probe_path", "resource_extension_probe_sha256") if key in calibration},
                   **{key: value for key, value in gate.items() if key not in {"status", "sample_id"}}}
    write_json(directory / "execution_performance.json", performance)
    if code:
        performance["status"] = "failed_model_execution"
        write_json(directory / "validation.json", performance)
        return performance
    raw = pd.read_csv(directory / "atom_predictions_full_precision.csv", keep_default_na=False)
    raw["chain"] = raw.chain.astype(str).replace("", "_")
    raw["author_residue_id"] = [author_id(a,b) for a,b in zip(raw["resno"],raw["insert"])]
    assert len(raw) == len(atoms)
    for original_key, raw_key in (("atom","elety"),("residue","resid"),("chain","chain"),("author_residue_id","author_residue_id")):
        assert atoms[original_key].astype(str).tolist() == raw[raw_key].astype(str).tolist()
    assert np.max(np.abs(atoms[["x","y","z"]].to_numpy()-raw[["x","y","z"]].to_numpy())) < 1e-9
    assert np.isfinite(raw.prediction.to_numpy(dtype=float)).all()
    ca = raw[raw.elety=="CA"].copy()
    ca["amino_acid"] = ca.resid.map(AA)
    assert identities(ca) == identities(truth)
    prediction = ca.prediction.to_numpy(dtype=np.float64)
    temporary = directory / "prediction.lifetime_v2.tmp.npz"
    np.savez_compressed(temporary, prediction=prediction, target=truth.target_rmsf.to_numpy(dtype=np.float64), valid=truth.valid_rmsf.to_numpy(dtype=bool),
                        chain=truth.chain.astype(str).to_numpy(dtype=str), author_residue_id=truth.author_residue_id.astype(str).to_numpy(dtype=str),
                        amino_acid=truth.amino_acid.astype(str).to_numpy(dtype=str), residue_index=truth.residue_index.to_numpy(dtype=np.int64))
    temporary.replace(directory / "prediction.npz")
    performance.update(status="complete", ca_positions_verified=len(prediction), all_atom_identities_and_coordinates_verified=True,
                       ca_identities_equal_frozen_truth_in_original_order=True, all_ca_predictions_finite=True,
                       output_quantity="normalized_relative_flexibility_CA_atom_score", absolute_RMSF_Angstrom=False, Hchi_supported=False,
                       reference_used_only_for_identity_and_evaluation=True, prediction_sha256=sha256(directory / "prediction.npz"),
                       full_precision_atoms_sha256=sha256(directory / "atom_predictions_full_precision.csv"),
                       rds_exact_output_sha256=sha256(directory / "atom_predictions_exact.rds"), official_pdb_output_sha256=sha256(directory / "proteinRMSF.pdb"))
    write_json(directory / "validation.json", performance)
    return performance


def archive_protocol_and_write(provenance, calibration, schedule):
    original_path = HERE / "run_protocol.json"
    prior_bytes = original_path.read_bytes()
    prior = json.loads(prior_bytes.decode("utf-8-sig"))
    if prior.get("algorithm_lifecycle_version") == "lifetime_v2":
        original_archive = Path(prior["official_v1_protocol_archive_path"])
        assert sha256(original_archive) == prior["official_v1_protocol_archive_sha256"]
    else:
        assert prior["upstream_sha256"]["FUNCTION_GDV.r"] == provenance["official_function_sha256"]
        original_archive = HERE / ("run_protocol.official_v1.archived."+hashlib.sha256(prior_bytes).hexdigest()[:16]+".json")
        if original_archive.exists():
            assert original_archive.read_bytes() == prior_bytes
        else:
            original_archive.write_bytes(prior_bytes)
    protocol = {
        "model": "FastProtFlex_official_GDV_linear_coefficients", "algorithm_lifecycle_version": "lifetime_v2",
        "Rscript": str(RSCRIPT), "R_libraries": R_LIBRARIES, "upstream_directory": str(UPSTREAM),
        "upstream_sha256": {name: sha256(UPSTREAM / name) for name in ("predict.r", "FUNCTION_GDV.r", "LICENSE", "README.md")},
        "official_v1_protocol_archive_path": str(original_archive), "official_v1_protocol_archive_sha256": sha256(original_archive),
        "official_algorithm_source_sha256": provenance["official_function_sha256"], "candidate_function_sha256": provenance["candidate_function_sha256"],
        "official_algorithm_source_path": str(UPSTREAM / "FUNCTION_GDV.r"),
        "candidate_function_path": str(CANDIDATE_ROOT / "candidate/FUNCTION_GDV.r"),
        "official_predict_sha256": provenance["official_predict_sha256"], "capture_wrapper_sha256": provenance["capture_wrapper_sha256"],
        "execution_entrypoint_path": str(Path(__file__).resolve()), "execution_entrypoint_sha256": sha256(__file__),
        "wrapper_path": str(Path(__file__).resolve()), "wrapper_sha256": sha256(__file__),
        "capture_wrapper_path": str(HERE / "capture_full_precision.R"),
        "lifecycle_validation_audit_path": calibration["validation_audit_path"],
        "resume_entrypoint_sha256": sha256(__file__), "lifecycle_validation_audit_sha256": calibration["validation_audit_sha256"],
        "resource_gate_version": calibration["version"], "resource_calibration": calibration,
        "resource_policy_version":RESOURCE_POLICY_VERSION, "policy_change_reason":POLICY_CHANGE_REASON,
        "truth_loading_policy":"Stream seven default-parser columns; retain only missing/incomplete record truth. Completed caches use input/source/output SHA plus PDB/NPZ identity. Final independent scorer must revalidate target/valid/residue_index against all frozen truth.",
        **{key:calibration[key] for key in ("resource_extension_audit_path", "resource_extension_audit_sha256", "resource_extension_probe_path", "resource_extension_probe_sha256") if key in calibration},
        "full_schedule_sha256": sha256(HERE / "full_schedule.csv"), "scheduled_records": len(schedule),
        "input_cohort_manifest_sha256": sha256(COHORT / "records.csv"), "new_case_source": "validated lifecycle_v2 only; no numerical algorithm changes",
        "cache_policy": "Reuse and hash/identity verify complete official_v1 and validated lifetime_v2 caches, without rewriting their predictions or validation metadata. Only missing/incomplete records receive a new forward pass.",
        "input_context": "complete_structure", "logical_cpus": 1, "library_threads": 1,
        "complete_structure_input": True, "independent_chain_inference": False,
        "branch_threshold_atoms": 10000, "regression_coefficients_unchanged": True, "no_model_fitting": True,
        "output_units": "normalized relative score; not Angstrom", "Hchi_supported": False,
        "extremes": "No missing structure above actual calibration max atom count can run. Root must review a separate resource-validation protocol; there is no automatic override.",
    }
    write_json(original_path, protocol)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--run", action="store_true")
    parser.add_argument("--validation-audit", type=Path)
    parser.add_argument("--resource-extension-audit", type=Path)
    parser.add_argument("--memory-wait-seconds", type=float, default=600)
    args = parser.parse_args()
    assert args.memory_wait_seconds >= 0
    provenance, schedule, records = inspect_sources()
    calibration = load_actual_calibration(args.validation_audit, provenance) if args.validation_audit else None
    if args.resource_extension_audit:
        assert calibration is not None, "--resource-extension-audit never replaces the actual three-source numerical-equivalence audit"
        calibration = load_resource_extension(args.resource_extension_audit, calibration, provenance, records)
    if not args.run:
        print(json.dumps({"status": "static_inspection_only_no_R_no_protocol_changes", "schedule_records": len(schedule),
                          "schedule_sha256": sha256(HERE / "full_schedule.csv"), "candidate_sha256": provenance["candidate_function_sha256"],
                          "actual_validation_audit_supplied": calibration is not None,
                          "calibration": calibration, "run_requirement": "Actual three-source exact-output audit required; original Fast/R stopped; no unmeasured extreme launch."}, indent=2))
        return 0
    assert calibration is not None, "--run requires --validation-audit from an actual successful three-source R validation"
    import numpy as np
    import pandas as pd
    import psutil
    assert RSCRIPT.is_file()
    assert not competing_jobs(psutil), "Stop/review original Fast runner and existing R processes before using this separate resume entry point"
    lock = HERE / "lifetime_v2.runner.lock"
    descriptor = os.open(lock, os.O_WRONLY|os.O_CREAT|os.O_EXCL)
    with os.fdopen(descriptor,"w") as handle:
        handle.write(json.dumps({"pid": os.getpid(), "script": str(Path(__file__).resolve()), "created_utc": datetime.now(timezone.utc).isoformat()}))
    try:
        # Preserve the original runner's copied-target parser. Scientific scores
        # remain the evaluator's authoritative frozen round-trip truth values.
        missing_ids = []
        for scheduled in schedule:
            validation = HERE / "records" / scheduled["sample_id"] / "validation.json"
            if not validation.exists() or json_read(validation).get("status") != "complete":
                missing_ids.append(scheduled["sample_id"])
        groups = load_missing_truth(records,missing_ids,pd)
        print(json.dumps({"truth_loading_policy":"stream_only_missing_records", "missing_records":len(missing_ids),
                          "truth_rows_retained":sum(len(frame) for frame in groups.values()),
                          "cached_target_valid_residue_index_revalidation":"independent_final_scorer_all_frozen_truth"}),flush=True)
        archive_protocol_and_write(provenance, calibration, schedule)
        results = []
        blocked = None
        for position, scheduled in enumerate(schedule, 1):
            record = records[scheduled["sample_id"]]
            cached = None
            try:
                atoms = validate_input(record, groups.get(record["sample_id"]), pd)
                assert len(atoms) == int(scheduled["protein_heavy_atoms"])
                cached = cached_prediction(record, atoms, provenance, np)
                result = cached if cached else launch_missing(record, groups[record["sample_id"]], atoms, provenance, calibration, args, np, pd, psutil)
            except Exception as error:
                # A corrupt cache is a failure, never permission to recompute it.
                # Preserve existing completed predictions and their metadata.
                result = {"sample_id": record["sample_id"], "dataset": record["dataset"],
                          "status": "failed_validation_or_execution", "error": repr(error),
                          "algorithm_lifecycle_version": "lifetime_v2", "cache_recomputed_after_failure": False}
            results.append(result)
            if not cached or position % 100 == 0 or position == len(schedule) or result["status"] != "complete":
                pd.DataFrame(results).to_csv(HERE / "full_results.csv", index=False)
            print(json.dumps({"progress": f"{position}/{len(schedule)}", "sample_id": record["sample_id"], "status": result["status"],
                              "cached": bool(cached), "lifecycle_version": result.get("cache_algorithm_lifecycle_version", result.get("algorithm_lifecycle_version"))}), flush=True)
            if result["status"].startswith("blocked"):
                blocked = result
                break
            if result["status"] != "complete":
                break
        complete = len(results)==3063 and all(row["status"]=="complete" for row in results)
        audit = {"mode":"full", "algorithm_lifecycle_version":"lifetime_v2", "records_requested":3063,
                 "records_attempted_or_cached":len(results), "records_complete":sum(row["status"]=="complete" for row in results),
                 "status":"complete" if complete else "blocked_resource_review_or_wait" if blocked else "some_failures",
                 "blocked":blocked, "lifecycle_validation_audit_sha256":calibration["validation_audit_sha256"],
                 "resource_gate_version":calibration["version"],
                 "resource_policy_version":RESOURCE_POLICY_VERSION,
                 "truth_records_loaded":len(groups), "truth_rows_retained":sum(len(frame) for frame in groups.values()),
                 "cached_frozen_truth_revalidated_this_invocation":False,
                 "final_scorer_frozen_truth_revalidation_required":True,
                 **{key:calibration[key] for key in ("resource_extension_audit_path", "resource_extension_audit_sha256", "resource_extension_probe_path", "resource_extension_probe_sha256") if key in calibration}, "results":results}
        write_json(HERE / "full_audit.lifetime_v2.json", audit)
        write_json(HERE / "full_audit.json", audit)
        return 0 if complete else 3 if blocked else 2
    finally:
        lock.unlink(missing_ok=True)


if __name__ == "__main__":
    raise SystemExit(main())
