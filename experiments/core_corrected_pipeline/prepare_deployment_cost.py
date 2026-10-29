"""Prepare an append-only, label-free deployment plan from frozen real inputs.

Preparation is not a benchmark result. Missing final weights are recorded,
never substituted with older weights. Local and remote paths are explicit.
"""
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
import uuid

HERE = Path(__file__).resolve().parent
FREEZE_SHA = "5b994bb59684e3a33c97f4a1801ee9a606e5dd856db9cfe23bf18515697707d3"
BF2_SHA = "1af6e13a0e0161b2b92fa6fb9cddd2dee8f3e71c284f5059c77d7a30b9c527b9"
PRIMARY_METHODS = ["full_single", "full_five", "intrachain_single", "static0", "static5", "backflip2"]
STRICT_METHODS = ["strict_full_single", "strict_intrachain_single", "strict_static", "backflip2"]


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def add_file(files, path, expected=None):
    path = Path(path).resolve()
    actual = sha(path)
    if expected is not None and actual != expected:
        raise ValueError("Frozen input/code/model hash mismatch: " + str(path))
    files[str(path)] = actual
    return str(path)


def build_plan(args):
    primary, external, static, bf2, source, strict = [p.resolve() for p in
        [args.primary_root, args.external_root, args.static_root, args.bf2_root,
         args.bf2_source, args.strict_root]]
    files = {}
    add_file(files, primary / "training_freeze.json", FREEZE_SHA)
    freeze = load(primary / "training_freeze.json")
    for name, expected in freeze["execution_hashes"].items():
        add_file(files, primary.joinpath(*name.replace("\\", "/").split("/")), expected)
    ext_freeze = load(external / "external_payload_freeze.json")
    add_file(files, external / "external_payload_freeze.json")
    add_file(files, external / "prespecified_protocol.json", ext_freeze["external_protocol_sha256"])
    add_file(files, external / "ready_test.jsonl")
    rows = [json.loads(x) for x in (external / "ready_test.jsonl").read_text(encoding="utf-8").splitlines()]
    if len(rows) != 32 or len({r["sample_id"] for r in rows}) != 32:
        raise ValueError("Frozen external cohort must have exactly 32 unique records")
    sys.path.insert(0, str(primary / "src"))
    import numpy as np
    from flexppi.data.structures import parse_pdb
    ext_hashes = {k.replace("\\", "/"): v for k, v in ext_freeze["data_sha256"].items()}
    bf_freeze = load(bf2 / "input_freeze.json")
    add_file(files, bf2 / "input_freeze.json")
    add_file(files, bf2 / "input_manifest.csv", bf_freeze["input_manifest_sha256"])
    add_file(files, bf2 / "residue_mapping.csv", bf_freeze["mapping_sha256"])
    chain_rows = [r for r in csv.DictReader((bf2 / "input_manifest.csv").open(encoding="utf-8"))
                  if r["role"] == "external"]
    mapping = {}
    for r in csv.DictReader((bf2 / "residue_mapping.csv").open(encoding="utf-8")):
        if r["role"] == "external":
            mapping.setdefault(r["input_id"], []).append(r)
    if len(chain_rows) != 95 or len(mapping) != 95:
        raise ValueError("All 95 frozen external chains must be present")
    records = []
    identities = {}
    for r in sorted(rows, key=lambda x: x["sample_id"]):
        # Resolve by frozen filename, not the platform-specific manifest path.
        pdb = external / "structures" / (r["sample_id"] + ".pdb")
        expected = ext_hashes["structures/" + pdb.name]
        add_file(files, pdb, expected)
        if expected != r["metadata"]["input_sha256"]:
            raise ValueError("Ready manifest and external freeze disagree")
        s = parse_pdb(str(pdb))
        n = s.num_residues
        if n != r["metadata"]["residue_count"] or n > 1024:
            raise ValueError("Prespecified input residue count changed")
        target = [c in r["target_chains"] for c in s.chain_ids]
        identity = dict(amino_acid=list(s.sequence), chain=s.chain_ids,
            author_residue_id=s.author_residue_ids, geometry_mask=s.frame_mask.tolist(),
            target_mask=target, chains=[])
        observed = [c for c in chain_rows if c["sample_id"] == r["sample_id"]]
        if {c["chain"] for c in observed} != set(s.unique_chains()):
            raise ValueError("Isolated and complex chain inventories disagree")
        visits = np.zeros(n, dtype=int)
        for c in observed:
            path = bf2 / c["pdb_file"]
            add_file(files, path, c["input_sha256"])
            m = sorted(mapping[c["input_id"]], key=lambda x: int(x["input_index"]))
            ids = s.indices_for_chain(c["chain"])
            isolated = parse_pdb(str(path))
            if [int(x["input_index"]) for x in m] != list(range(len(ids))):
                raise ValueError("Chain input index is not complete and ordered")
            if [int(x["record_residue_index"]) for x in m] != ids.tolist():
                raise ValueError("Complex-to-chain residue map changed")
            if [x["amino_acid"] for x in m] != [s.sequence[i] for i in ids]:
                raise ValueError("Mapped amino-acid identities disagree")
            if [x["author_residue_id"] for x in m] != [s.author_residue_ids[i] for i in ids]:
                raise ValueError("Mapped author residue identities disagree")
            if isolated.sequence != s.chain_sequence(c["chain"]) or int(c["n_residues"]) != len(ids):
                raise ValueError("Isolated and complex sequences disagree")
            if not np.array_equal(isolated.backbone_mask, s.backbone_mask[ids]) or not np.array_equal(isolated.backbone, s.backbone[ids]):
                raise ValueError("Isolated and complex backbone inputs disagree")
            visits[ids] += 1
            identity["chains"].append(dict(input_id=c["input_id"], chain=c["chain"],
                structure_path=str(path.resolve()), record_indices=ids.tolist()))
        if not (visits == 1).all():
            raise ValueError("Frozen isolated inputs must cover every complex residue once")
        identities[r["sample_id"]] = identity
        # No label_path, trajectory, target values, or feature caches are exposed.
        records.append(dict(sample_id=r["sample_id"], structure_path=str(pdb.resolve()),
            target_chains=r["target_chains"], metadata={"residue_count": n},
            input_sha256=expected, chain_count=len(observed)))
    checkpoint = source / "models/backflip-2.1/backflip-2.1.ckpt"
    add_file(files, checkpoint, BF2_SHA)
    add_file(files, source / "models/backflip-2.1/config.yaml")
    for p in source.rglob("*.py"):
        if "__pycache__" not in p.parts:
            add_file(files, p)
    for name in ["static_feature_definitions.py", "deployment_cost_protocol_zh.md",
                 "deployment_cost_backflip_adapter.py", "deployment_cost_neural_static_adapter.py",
                 "deployment_cost_worker.py", "run_deployment_cost.py", "prepare_deployment_cost.py"]:
        add_file(files, HERE / name)
    for path in args.dependency_provenance:
        add_file(files, path)
    schema = add_file(files, static / "training_schema.json")
    static_schema = load(schema)
    add_file(files, HERE / "static_feature_definitions.py", static_schema["feature_source_sha256"])
    matched = load(static / "matched_external_membership_audit.json")
    add_file(files, static / "matched_external_membership_audit.json")
    if matched["training_freeze_sha256"] != FREEZE_SHA or matched["new_fits"] != 0:
        raise ValueError("Matched static training provenance changed")
    add_file(files, static / "external_schema.json", matched["external_schema_sha256"])
    members = {m["fold"]: m for m in matched["members"]}
    if set(members) != set(range(5)):
        raise ValueError("Exactly five matched static training members are required")
    for fold, member in members.items():
        train_file = primary / f"splits/fold_{fold}_train.jsonl"
        test_file = primary / f"splits/fold_{fold}_test.jsonl"
        add_file(files, train_file, member["train_manifest_sha256"])
        add_file(files, test_file, member["test_manifest_sha256"])
        declared = [json.loads(x)["sample_id"] for x in train_file.read_text(encoding="utf-8").splitlines()]
        if len(declared) != member["train_records"] or set(declared) != set(member["training_members"]):
            raise ValueError("Matched static and neural training-member sets disagree")
    methods = {}
    missing = []
    def neural(method, variant, folds, is_strict=False):
        configs, weights = [], []
        for fold in folds:
            config = (strict / f"configs/strict_{variant}.yaml" if is_strict
                      else primary / f"configs/{variant}_fold{fold}.yaml")
            configs.append(add_file(files, config))
            member = "strict_" + variant if is_strict else variant
            weight = args.weights_root.resolve() / f"{member}/fold_{fold}/epoch_040.pt"
            weights.append(str(weight))
            if weight.exists():
                add_file(files, weight)
            else:
                missing.append(str(weight))
        methods[method] = dict(kind="neural", checkpoint_paths=weights, config_paths=configs,
            protocol="strict_single" if is_strict else "primary",
            output_metrics=["rmsf", "chi_entropy"], software_environment="neural")
    neural("full_single", "full", [0])
    neural("full_five", "full", list(range(5)))
    neural("intrachain_single", "intrachain", [0])
    for name, folds in [("static0", [0]), ("static5", list(range(5)))]:
        paths = []
        for fold in folds:
            paths.append({metric: add_file(files, static / f"static_fold{fold}_{metric}.joblib",
                members[fold]["static_checkpoint_sha256"][metric]) for metric in ["rmsf", "chi_entropy"]})
        methods[name] = dict(kind="static", static_paths=paths, schema_path=schema,
            output_metrics=["rmsf", "chi_entropy"], protocol="primary", software_environment="neural")
    if args.include_strict:
        strict_freeze = load(strict / "strict_external_freeze.json")
        if strict_freeze["primary_freeze_sha256"] != FREEZE_SHA or strict_freeze["epoch"] != 40:
            raise ValueError("Strict frozen training contract changed")
        add_file(files, strict / "strict_external_freeze.json")
        for name, expected in strict_freeze["execution_hashes"].items():
            add_file(files, strict.joinpath(*name.replace("\\", "/").split("/")), expected)
        strict_audit = load(static / "strict_run_audit.json")
        add_file(files, static / "strict_run_audit.json")
        if (strict_audit["status"] != "complete" or strict_audit["training_records"] != 2431
            or strict_audit["external_records"] != 32 or strict_audit["no_test_tuning"] is not True
            or strict_audit["strict_freeze_sha256"] != sha(strict / "strict_external_freeze.json")
            or strict_audit["feature_schema_sha256"] != sha(schema)):
            raise ValueError("Strict static control does not match the frozen 2431-record protocol")
        neural("strict_full_single", "full", [0], True)
        neural("strict_intrachain_single", "intrachain", [0], True)
        methods["strict_static"] = dict(kind="static", static_paths=[{metric:
            add_file(files, static / f"strict_static_{metric}.joblib") for metric in ["rmsf", "chi_entropy"]}],
            schema_path=schema, output_metrics=["rmsf", "chi_entropy"], protocol="strict_single",
            software_environment="neural")
    methods["backflip2"] = dict(kind="backflip2", source_root=str(source),
        checkpoint_path=str(checkpoint), output_metrics=["rmsf"],
        protocol="published_comparator", software_environment="backflip2")
    return dict(status="ready_for_benchmark" if not missing else "preparation_only_missing_final_weights",
        created_utc=datetime.now(timezone.utc).isoformat(), new_fits=0,
        accuracy_evaluation=False, final_cost_results=False, family_components=1,
        cpu_threads=2, interop_threads=1, repetitions=3, batch_size=1,
        primary_methods=PRIMARY_METHODS, strict_methods=STRICT_METHODS if args.include_strict else [],
        order_rule="Rotate the fixed method list left by round index (0,1,2); serial processes",
        worker_rule="One fresh process per method and round; imports, model load and fixed warmup separately timed",
        serialization_scope="Prediction arrays and residue identities in compressed NPZ; OS-buffered close included, durable fsync excluded",
        primary_root=str(primary), methods=methods, records=records, identities=identities,
        files_sha256=files, missing_final_weights=sorted(set(missing)),
        python_executables={"neural": str(args.neural_python), "backflip2": str(args.bf2_python)},
        dependency_roots=[str(p.resolve()) for p in args.dependency_root],
        completion_gates=[str(primary / "queue_state.json"), str(primary / "strict_queue_state.json"),
            str(primary / "external_evaluation/evaluation_audit.json"),
            str(primary / "strict_external_evaluation/evaluation_audit.json")],
        training_freeze_sha256=FREEZE_SHA,
        scope="Prepared deployment input -> parse/features -> CPU forward -> conversion/identity validation -> NPZ serialization; hashes/ground truth/accuracy/initialization excluded",
        backflip2_identity_contract="Official outputs have no original author IDs. Frozen isolated-chain sequence/backbone hashes and ordered maps determine residue identities; output lengths and exactly-once complex coverage are checked for every timed call.",
        limitation="32 related complexes constitute one family component. Shared-host load and software environments are disclosed. No independent-family CI or p value.")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--primary-root", type=Path, default=HERE / "training_payload")
    parser.add_argument("--external-root", type=Path, default=HERE / "external_payload/external")
    parser.add_argument("--static-root", type=Path, default=HERE / "static_corrected")
    parser.add_argument("--strict-root", type=Path, default=HERE / "strict_external_payload")
    parser.add_argument("--weights-root", type=Path, default=HERE / "trained_checkpoints")
    parser.add_argument("--bf2-root", type=Path, default=HERE / "backflip2_corrected_inputs")
    parser.add_argument("--bf2-source", type=Path, default=HERE.parent / "backflip2/source")
    parser.add_argument("--neural-python", type=Path, default=Path(sys.executable))
    parser.add_argument("--bf2-python", type=Path, default=Path(sys.executable))
    parser.add_argument("--include-strict", action="store_true")
    parser.add_argument("--dependency-root", type=Path, action="append", default=[])
    parser.add_argument("--dependency-provenance", type=Path, action="append", default=[])
    parser.add_argument("--output-root", type=Path, default=HERE / "deployment_cost/preparations")
    args = parser.parse_args()
    plan = build_plan(args)
    tag = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ") + "_" + uuid.uuid4().hex[:8]
    out = args.output_root.resolve() / tag
    out.mkdir(parents=True, exist_ok=False)
    path = out / "plan.json"
    path.write_text(json.dumps(plan, indent=2, ensure_ascii=True), encoding="utf-8")
    print(json.dumps({"plan": str(path), "status": plan["status"], "records": len(plan["records"]),
        "chains": sum(r["chain_count"] for r in plan["records"]),
        "missing_final_weights": len(plan["missing_final_weights"]),
        "new_fits": 0, "benchmark_executed": False}, indent=2))


if __name__ == "__main__":
    main()
