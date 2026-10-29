"""Audit official chain predictions and score corrected, identity-matched RMSF.

No calibration, rescaling, fitted threshold or error-based sample filtering is
performed. Ordinary and purged static fits are reported as separate protocols.
The external cohort has one family component; this script produces no CI or p.
"""
from pathlib import Path
import json
import hashlib
import numpy as np
import pandas as pd
from run_corrected_static import HERE, OUT as STATIC, metric_values, parse_pdb, target_mask

INPUT = HERE / "backflip2_corrected_inputs"
PRED = HERE / "backflip2_corrected_completed/predictions"
OUT = HERE / "backflip2_corrected_scoring"


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def main():
    OUT.mkdir(exist_ok=True)
    freeze = json.loads((INPUT / "input_freeze.json").read_text())
    metadata = json.loads((PRED / "inference_metadata.json").read_text())
    assert sha(INPUT / "input_manifest.csv") == freeze["input_manifest_sha256"] == metadata["manifest_sha256"]
    assert sha(INPUT / "residue_mapping.csv") == freeze["mapping_sha256"]
    assert metadata["checkpoint_sha256"] == "1af6e13a0e0161b2b92fa6fb9cddd2dee8f3e71c284f5059c77d7a30b9c527b9"
    assert metadata["source_commit"] == "887c6a73b216e3fe2ff9f551ae6e2c983c21b292"
    assert metadata["device"] == "cpu" and metadata["covariance_unit"] == "nm2"
    manifest = pd.read_csv(INPUT / "input_manifest.csv", keep_default_na=False)
    mapping = pd.read_csv(INPUT / "residue_mapping.csv", keep_default_na=False,
                          dtype={"author_residue_id": str, "chain": str})
    frozen_records = json.loads((INPUT / "record_manifest.json").read_text())
    assert len(manifest) == 5068 and len(frozen_records) == 3095
    assert manifest.input_id.is_unique
    groups = dict(tuple(mapping.groupby("input_id", sort=False)))
    status = pd.DataFrame([json.loads(line) for line in (PRED / "status.jsonl").read_text().splitlines()])
    assert len(status) == 5068 and status.input_id.is_unique
    assert set(status.input_id) == set(manifest.input_id)
    assert (status.status == "complete").all(), "Disclose inference failures before making a complete-cohort comparison"
    assert len(list(PRED.glob("*.npz"))) == 5068
    predictions, output_hashes, conversion_errors = {}, [], []
    for entry in manifest.to_dict("records"):
        key = entry["input_id"]
        assert sha(INPUT / entry["pdb_file"]) == entry["input_sha256"]
        path = PRED / (key + ".npz")
        with np.load(path) as z:
            rmsf = z["rmsf_angstrom"].copy()
            covariance = z["covariance_nm2"]
            assert rmsf.shape == (entry["n_residues"],)
            assert covariance.shape == (entry["n_residues"], 3, 3)
            assert np.isfinite(rmsf).all() and np.isfinite(covariance).all()
            trace = np.trace(covariance, axis1=-2, axis2=-1)
            assert (trace >= 0).all()
            converted = 10 * np.sqrt(trace)
            assert np.allclose(rmsf, converted, atol=2e-6, rtol=2e-6)
            conversion_errors.append(float(np.max(np.abs(rmsf - converted))))
        chain_map = groups[key].sort_values("input_index")
        assert np.array_equal(chain_map.input_index, np.arange(entry["n_residues"]))
        assert (chain_map.sample_id == entry["sample_id"]).all()
        assert (chain_map.role == entry["role"]).all() and (chain_map.chain == entry["chain"]).all()
        predictions[key] = rmsf
        output_hashes.append(dict(input_id=key, sha256=sha(path), n_residues=entry["n_residues"]))
    record_lookup = {(row["role"], row["record"]["sample_id"]): row for row in frozen_records}
    overlap = pd.read_csv(HERE / "sequence_audit/record_overlap_audit.csv", keep_default_na=False)
    per_record, residue_tables, coverage = [], [], []
    for role in ["training", "external"]:
        schema = json.loads((STATIC / (role + "_schema.json")).read_text())
        with np.load(STATIC / (role + "_features.npz")) as z:
            target, valid = z["target"], z["valid"]
        bf = np.full(len(target), np.nan, dtype=np.float64)
        identity = []
        for info in schema["records"]:
            frozen = record_lookup[(role, info["sample_id"])]
            record = frozen["record"]
            assert sha(record["structure_path"]) == info["input_sha256"] == frozen["input_sha256"]
            assert sha(record["label_path"]) == info["label_sha256"] == frozen["label_sha256"]
            s = parse_pdb(record["structure_path"])
            n = s.num_residues
            assert n == info["stop"] - info["start"]
            with np.load(record["label_path"]) as z:
                recomputed = z["mask"][:, :2] & s.frame_mask[:, None] & target_mask(s, record["target_chains"])[:, None]
                assert np.array_equal(recomputed, valid[info["start"]:info["stop"]])
                assert np.array_equal(z["values"][:, :2], target[info["start"]:info["stop"]], equal_nan=True)
            assembled = np.full(n, np.nan, dtype=np.float64)
            visits = np.zeros(n, dtype=np.int32)
            entries = manifest[(manifest.role == role) & (manifest.sample_id == info["sample_id"])]
            assert set(entries.chain) == set(s.unique_chains())
            for entry in entries.to_dict("records"):
                chain_map = groups[entry["input_id"]].sort_values("input_index")
                ids = chain_map.record_residue_index.to_numpy(int)
                assert np.array_equal(ids, s.indices_for_chain(entry["chain"]))
                assert chain_map.amino_acid.tolist() == [s.sequence[i] for i in ids]
                assert chain_map.author_residue_id.tolist() == [s.author_residue_ids[i] for i in ids]
                assembled[ids] = predictions[entry["input_id"]]
                visits[ids] += 1
            assert (visits == 1).all() and np.isfinite(assembled).all()
            bf[info["start"]:info["stop"]] = assembled
            identity.append(pd.DataFrame(dict(sample_id=info["sample_id"], dataset=info["dataset"],
                biological_system=info["biological_system"], fold=info["fold"], residue_index=np.arange(n),
                amino_acid=list(s.sequence), chain=s.chain_ids, author_residue_id=s.author_residue_ids)))
            coverage.append(dict(role=role, sample_id=info["sample_id"], chains=len(entries), residues=n,
                valid_rmsf=int(recomputed[:, 0].sum()), valid_chi_entropy=int(recomputed[:, 1].sum()),
                backflip2_complete=True))
        assert np.isfinite(bf).all()
        np.savez_compressed(OUT / (role + "_aligned.npz"), prediction=bf, target=target[:, 0], valid=valid[:, 0])
        suffix = "oof" if role == "training" else "external"
        variants = {"backflip2_official": bf}
        for name in ["static_sequence_complex", "aa_only"]:
            with np.load(STATIC / (name + "_" + suffix + ".npz")) as z:
                variants[name] = z["prediction"][:, 0]
        protocols = {"ordinary_" + suffix: variants}
        if role == "external":
            strict = {"backflip2_official": bf}
            for name in ["strict_static_sequence_complex", "strict_aa_only"]:
                with np.load(STATIC / (name + "_external.npz")) as z:
                    strict[name] = z["prediction"][:, 0]
            protocols["purged_external"] = strict
        source = pd.concat(identity, ignore_index=True)
        assert len(source) == len(target)
        source["role"] = role
        source["target_rmsf"] = target[:, 0]
        source["valid_rmsf"] = valid[:, 0]
        for protocol, predictions_for_protocol in protocols.items():
            for name, p in predictions_for_protocol.items():
                assert np.isfinite(p[valid[:, 0]]).all()
                source[name] = p
                for info in schema["records"]:
                    indices = np.arange(info["start"], info["stop"])
                    indices = indices[valid[indices, 0]]
                    assert len(indices) > 0
                    per_record.append(dict(role=role, protocol=protocol, variant=name, metric="rmsf",
                        sample_id=info["sample_id"], dataset=info["dataset"], fold=info["fold"],
                        biological_system=info["biological_system"], **metric_values(target[indices, 0], p[indices])))
        source = source.merge(overlap[["sample_id", "role", "novelty", "component_id"]],
                              on=["sample_id", "role"], how="left", validate="many_to_one")
        assert source.component_id.notna().all()
        source.to_csv(OUT / (role + "_per_residue.csv.gz"), index=False, compression="gzip")
        residue_tables.append(source)
        print(role, len(schema["records"]), int(valid[:, 0].sum()), flush=True)
    per = pd.DataFrame(per_record).merge(overlap[["sample_id", "role", "novelty", "component_id"]],
                                        on=["sample_id", "role"], how="left", validate="many_to_one")
    per.to_csv(OUT / "per_record_metrics.csv", index=False)
    numeric = ["mae", "rmse", "bias", "spearman", "pearson"]
    # Both estimands are retained: record macro matches the existing exporter;
    # biology macro gives each source-specific biological system one vote.
    rows = []
    for keys, group in per.groupby(["role", "protocol", "variant", "dataset"]):
        for estimand in ["record_macro", "biological_system_macro"]:
            units = group if estimand == "record_macro" else group.groupby("biological_system")[numeric].mean()
            rows.append(dict(zip(["role", "protocol", "variant", "dataset"], keys), estimand=estimand,
                records=len(group), biological_systems=group.biological_system.nunique(),
                family_components=group.component_id.nunique(), residues=int(group["count"].sum()),
                **{m: float(units[m].mean()) for m in numeric},
                spearman_records_defined=int(group.spearman.notna().sum())))
    summary = pd.DataFrame(rows)
    summary.to_csv(OUT / "source_macro_metrics.csv", index=False)
    per.groupby(["protocol", "variant", "dataset", "fold"])[numeric].mean().to_csv(OUT / "source_fold_macro_metrics.csv")
    pd.DataFrame(coverage).to_csv(OUT / "coverage_audit.csv", index=False)
    chain_times = manifest.merge(status[["input_id", "seconds", "min_rmsf", "max_rmsf"]], on="input_id", validate="one_to_one")
    chain_times.to_csv(OUT / "cpu_chain_cost.csv", index=False)
    costs = []
    for role, group in chain_times.groupby("role"):
        costs.append(dict(role=role, chains=len(group), residues=int(group.n_residues.sum()),
            timed_seconds_sum=float(group.seconds.sum()), median_seconds=float(group.seconds.median()),
            p05_seconds=float(group.seconds.quantile(.05)), p95_seconds=float(group.seconds.quantile(.95)),
            max_seconds=float(group.seconds.max())))
    pd.DataFrame(costs).to_csv(OUT / "cpu_cost_summary.csv", index=False)
    audit = dict(status="complete", records=3095, chains=5068, inference_failures=0,
        rmsf_valid_training=961321, rmsf_valid_external=18748,
        full_identity_and_input_hash_verification=True, maximum_covariance_conversion_difference=max(conversion_errors),
        no_fitted_calibration=True, no_error_based_exclusions=True,
        external_family_components=int(per.loc[per.role == "external", "component_id"].nunique()),
        external_uncertainty="One connected family component; no independent-family CI or p value is identifiable.",
        ordinary_training_members=3063, purged_training_members=2431,
        training_static_protocol="Held-out fold only; BackFlip2 uses its released external weights, not retrained folds.",
        external_static_protocol="Ordinary controls train on all3063; purged controls train on2431. Separate comparisons.",
        backflip2_pretraining_overlap="Released comparator training overlap and MD protocol differences remain limitations.",
        cost_scope="One completed pass, CPU2threads; per-chain input hash, official predict_from_pdb, conversion and NPZ serialization; excludes model initialization and input preparation. Concurrent GPU training and host workload preclude an unqualified speed ratio.",
        official_inference_metadata=metadata, output_hashes=output_hashes,
        mapping_sha256=sha(INPUT / "residue_mapping.csv"), source_sha256=sha(__file__),
        training_freeze_sha256=sha(HERE / "training_payload/training_freeze.json"),
        strict_freeze_sha256=sha(HERE / "strict_external_payload/strict_external_freeze.json"),
        external_feature_schema_sha256=sha(STATIC / "external_schema.json"),
        search_audit_sha256=sha(HERE / "sequence_audit/record_overlap_audit.csv"))
    assert audit["external_family_components"] == 1
    assert sum(r["valid_rmsf"] for r in coverage if r["role"] == "training") == audit["rmsf_valid_training"]
    assert sum(r["valid_rmsf"] for r in coverage if r["role"] == "external") == audit["rmsf_valid_external"]
    (OUT / "scoring_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
