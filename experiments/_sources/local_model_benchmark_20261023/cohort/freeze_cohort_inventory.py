"""Read-only validation of the existing cohort; write new manifests only here.

No inference, fitting, calibration, sequence clustering, network access, or
changes to the original structures, predictions, labels, or frozen hashes.
"""
from pathlib import Path
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import csv
import hashlib
import json
import sys

import numpy as np
import pandas as pd

sys.stdout.reconfigure(encoding="utf-8")
ROOT = Path(__file__).resolve().parent
EXPERIMENTS = Path("E:/Project/FlexPPI/experiments")
BASE = EXPERIMENTS / "_sources/baselines"
COHORT = BASE / "cohort"
CORE = EXPERIMENTS / "core_corrected_pipeline"
FAMILY = EXPERIMENTS / "_sources/family_folds"
ANNOTATION = EXPERIMENTS / "_sources/functional_annotations_complete"
assert ROOT == (EXPERIMENTS / "_sources/local_model_benchmark_20261023/cohort").resolve()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_csv(name, rows):
    assert rows
    with (ROOT / name).open("w", encoding="utf-8-sig", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path, **kwargs):
    return pd.read_csv(path, keep_default_na=False, float_precision="round_trip", **kwargs)


def main():
    freeze = read_json(COHORT / "cohort_freeze.json")
    inputs = {COHORT / "cohort_freeze.json", Path(__file__)}
    for name, expected in freeze["files"].items():
        path = COHORT / name
        assert sha(path) == expected, (name, "current local freeze mismatch")
        inputs.add(path)
    structures = read_csv(COHORT / "structure_manifest.tsv", sep="\t")
    chains = read_csv(COHORT / "chain_manifest.csv")
    truth = read_csv(COHORT / "truth_residues.csv.gz", dtype={"chain": str, "author_residue_id": str})
    index = read_json(COHORT / "record_index.json")
    schema_path = CORE / "static_corrected/training_schema.json"
    schema = read_json(schema_path)["records"]
    family_path = FAMILY / "record_family_folds.csv"
    family = read_csv(family_path).set_index("sample_id")
    inputs.update([schema_path, family_path, FAMILY / "training_freeze.json"])
    assert len(structures) == len(index) == len(schema) == len(family) == 3063
    assert structures.sample_id.is_unique and len(chains) == 4973
    assert len(truth) == 961321 and truth.sample_id.nunique() == 3063
    assert int(truth.valid_rmsf.sum()) == 961321
    assert int(truth.valid_chi_entropy.sum()) == 824988
    assert not truth.duplicated(["sample_id", "chain", "author_residue_id"]).any()
    assert family.component_id.nunique() == 1662
    sample_set = set(structures.sample_id)
    assert sample_set == set(family.index) == {x["sample_id"] for x in index} == {x["sample_id"] for x in schema}
    metadata_path = CORE / "neural_exports/evaluation/prediction_manifest.jsonl"
    inputs.add(metadata_path)
    metadata = {}
    for line in metadata_path.open(encoding="utf-8"):
        row = json.loads(line)
        if row["variant"] == "full":
            assert row["sample_id"] not in metadata
            metadata[row["sample_id"]] = row
    assert set(metadata) == sample_set
    annotations = read_csv(ANNOTATION / "record_annotations.csv", encoding="utf-8-sig").set_index("sample_id")
    inputs.add(ANNOTATION / "record_annotations.csv")
    source_entries = {}
    for path in sorted(ANNOTATION.glob("batch_*.json")):
        inputs.add(path)
        for entry in read_json(path)["response"]["data"]["entries"]:
            if entry:
                source_entries[entry["rcsb_id"].upper()] = entry

    static_paths = {method: CORE / f"static_corrected/{name}_oof.npz" for method, name in [("Static", "static_sequence_complex"), ("AA", "aa_only")]}
    static = {}
    for method, path in static_paths.items():
        inputs.add(path)
        with np.load(path, allow_pickle=False) as saved:
            static[method] = saved["prediction"].copy()
        assert static[method].shape == (961321, 2)
    backflip_path = CORE / "backflip2_corrected_scoring/training_aligned.npz"
    inputs.add(backflip_path)
    with np.load(backflip_path, allow_pickle=False) as saved:
        backflip = {key: saved[key].copy() for key in saved.files}
    assert backflip["prediction"].shape == (961321,)
    assert backflip["valid"].all() and np.isfinite(backflip["prediction"]).all()

    truth_by = dict(tuple(truth.groupby("sample_id", sort=False)))
    chain_by = dict(tuple(chains.groupby("sample_id", sort=False)))
    structure_by = structures.set_index("sample_id")
    index_by = {x["sample_id"]: x for x in index}
    records, chain_rows, predictions, candidate_rows = [], [], [], []
    hist = Counter()
    for info in schema:
        sid = info["sample_id"]
        current = truth_by[sid].sort_values("residue_index")
        source = structure_by.loc[sid]
        old_index = index_by[sid]
        chain_group = chain_by[sid]
        model_record = metadata[sid]["record"]
        model_meta = model_record.get("metadata", {})
        fam = family.loc[sid]
        start, stop = int(info["start"]), int(info["stop"])
        assert start == int(old_index["start"]) and stop == int(old_index["stop"])
        assert stop - start == len(current) == int(source.n_residues)
        assert len(chain_group) == int(source.n_chains)
        assert np.array_equal(current.residue_index, np.arange(len(current)))
        assert int(source.fold) == int(info["fold"]) == int(fam.original_fold)
        assert np.array_equal(backflip["target"][start:stop], current.target_rmsf.to_numpy(dtype=np.float32), equal_nan=True)
        structure_path = COHORT / "structures" / Path(source.container_pdb_path).name
        assert structure_path.is_file()
        inputs.add(structure_path)
        n_chains = len(chain_group)
        group = "single_chain" if n_chains == 1 else "tcr_pmhc_itrace" if source.cohort == "itrace" else "other_multichain"
        is_multi = n_chains > 1
        is_tcr = source.cohort == "itrace"
        assert not is_tcr or (is_multi and model_record.get("source_group") == "itrace_complex")
        ppi_annotation = model_meta.get("ppi_classification", "unknown")
        r = dict(sample_id=sid, dataset=str(source.cohort), biological_system=str(source.biological_system),
                 original_fold=int(source.fold), family_fold=int(fam.fold), family_component_id=str(fam.component_id),
                 biological_cluster_id=model_meta.get("cluster_id", ""), source_group=model_record.get("source_group", ""),
                 main_group_candidate=group, group_status="assumption_pending_optional_user_preference",
                 n_chains=n_chains, n_residues=len(current), n_valid_rmsf=int(current.valid_rmsf.sum()),
                 n_valid_hchi=int(current.valid_chi_entropy.sum()), is_multichain=is_multi, is_tcr_pmhc_itrace=is_tcr,
                 ppi_classification_original=ppi_annotation,
                 explicit_high_confidence_ppi_or_tcr=ppi_annotation == "high_confidence_ppi" or is_tcr,
                 source_designated_complex=model_record.get("source_group", "").endswith("_complex"),
                 observed_distinct_chain_sequences=int(chain_group.sequence.nunique()),
                 observed_all_chain_sequences_identical=chain_group.sequence.nunique() == 1,
                 structure_path=str(structure_path), structure_sha256_expected=str(source.input_sha256),
                 original_structure_path=model_record["structure_path"], original_input_sha256=str(source.original_input_sha256),
                 original_label_path=model_record["label_path"], original_label_sha256=str(source.label_sha256),
                 reference_table_path=str(COHORT / "truth_residues.csv.gz"), reference_sample_filter=sid,
                 aligned_array_start=start, aligned_array_stop=stop,
                 sequence_structure_sha256=hashlib.sha256("|".join(chain_group.sequence).encode()).hexdigest())
        hist[group] += 1
        for c in chain_group.to_dict("records"):
            local = COHORT / "chains" / Path(c["pdb_file"]).name
            assert local.is_file()
            inputs.add(local)
            chain_rows.append(dict(sample_id=sid, input_id=c["input_id"], chain=c["chain"], n_residues=int(c["n_residues"]),
                                   sequence=c["sequence"], chain_sequence_sha256=hashlib.sha256(c["sequence"].encode()).hexdigest(),
                                   structure_path=str(local), structure_sha256_expected=c["input_sha256"],
                                   original_structure_path=c["pdb_file"], main_group_candidate=group))
        # Diagnostic candidates only. Canonical input chains can be renamed;
        # author-label coincidence is not sufficient to freeze bio-oligomer status.
        annotation = annotations.loc[sid]
        pdb_ids = str(annotation.pdb_ids).split(";")
        entry = source_entries.get(pdb_ids[0]) if len(pdb_ids) == 1 else None
        entity_map = defaultdict(set)
        entities = (entry or {}).get("polymer_entities") or []
        uniprots = set()
        for entity in entities:
            ids = entity.get("rcsb_polymer_entity_container_identifiers") or {}
            uniprots.update(ids.get("uniprot_ids") or [])
            for auth_chain in ids.get("auth_asym_ids") or []:
                entity_map[auth_chain].add(entity["rcsb_id"])
        matched = all(len(entity_map[c]) == 1 for c in chain_group.chain)
        matched_entities = sorted({next(iter(entity_map[c])) for c in chain_group.chain}) if matched else []
        candidate = "single_chain" if not is_multi else "unresolved_multichain"
        if is_multi and matched:
            candidate = "candidate_one_entity_multichain" if len(matched_entities) == 1 else "candidate_multi_entity_complex"
        elif is_multi and len(entities) == 1:
            candidate = "candidate_one_source_entity_multichain"
        candidate_rows.append(dict(sample_id=sid, pdb_ids=annotation.pdb_ids, annotation_status=annotation.annotation_status,
                                   input_chain_ids=";".join(chain_group.chain), source_entity_count=len(entities),
                                   source_uniprot_ids=";".join(sorted(uniprots)),
                                   input_chain_labels_match_all_source_auth_labels=matched,
                                   matched_source_entity_ids=";".join(matched_entities),
                                   candidate_biological_architecture=candidate,
                                   classification_status="diagnostic_only_not_adopted_input_chain_renaming_not_resolved",
                                   pdb_classification=annotation.pdb_classification, protein_descriptions=annotation.protein_descriptions))
        for method, dirname in [("Full", "full"), ("Intra", "intrachain")]:
            path = CORE / "neural_exports/evaluation" / dirname / (sid + ".npz")
            assert path.is_file()
            inputs.add(path)
            with np.load(path, allow_pickle=False) as saved:
                values = saved["mean"]
                valid = saved["valid"]
                target = saved["target"]
                assert values.shape == (len(current), 2)
                assert saved["amino_acid"].astype(str).tolist() == current.amino_acid.tolist()
                assert saved["chain"].astype(str).tolist() == current.chain.tolist()
                assert saved["author_residue_id"].astype(str).tolist() == current.author_residue_id.tolist()
                for j, task in enumerate(["rmsf", "chi_entropy"]):
                    mask = current["valid_" + task].to_numpy(bool)
                    assert np.array_equal(valid[:, j], mask)
                    assert np.array_equal(target[mask, j], current.loc[mask, "target_" + task].to_numpy(dtype=target.dtype))
                    assert np.isfinite(values[mask, j]).all()
            r[dirname + "_prediction_path"] = str(path)
            predictions.append(dict(sample_id=sid, method=method, task_support="rmsf;chi_entropy", local_raw_available=True,
                                    path=str(path), array_key="mean", row_start=0, row_stop=len(current),
                                    rmsf_column=0, hchi_column=1, n_valid_rmsf=int(current.valid_rmsf.sum()),
                                    n_valid_hchi=int(current.valid_chi_entropy.sum()), coverage_basis="this_run_array_identity_and_finite_checks",
                                    input_context="complete_structure" if method == "Full" else "within_chain_encoder"))
        for method in ["Static", "AA"]:
            values = static[method][start:stop]
            for j, task in enumerate(["rmsf", "chi_entropy"]):
                assert np.isfinite(values[current["valid_" + task].to_numpy(bool), j]).all()
            predictions.append(dict(sample_id=sid, method=method, task_support="rmsf;chi_entropy", local_raw_available=True,
                                    path=str(static_paths[method]), array_key="prediction", row_start=start, row_stop=stop,
                                    rmsf_column=0, hchi_column=1, n_valid_rmsf=int(current.valid_rmsf.sum()),
                                    n_valid_hchi=int(current.valid_chi_entropy.sum()), coverage_basis="this_run_array_finite_checks_and_frozen_training_schema",
                                    input_context="complete_structure_feature_summaries" if method == "Static" else "training_fold_amino_acid_type_mean"))
        predictions.append(dict(sample_id=sid, method="BackFlip2_2.1", task_support="rmsf", local_raw_available=True,
                                path=str(backflip_path), array_key="prediction", row_start=start, row_stop=stop,
                                rmsf_column="", hchi_column="", n_valid_rmsf=int(current.valid_rmsf.sum()),
                                n_valid_hchi=0, coverage_basis="this_run_finite_and_target_checks_frozen_training_schema_identity",
                                input_context="independent_chain_predictions_reassembled_per_record"))
        records.append(r)
    assert hist == Counter(single_chain=2193, other_multichain=626, tcr_pmhc_itrace=244)
    print("Validated records, identities and Full/Intra/Static/AA/BackFlip arrays:", dict(hist), flush=True)

    # The locally retained full-residue ProDy file is authoritative for ANM/GNM.
    prody_path = BASE / "prody_common_evaluation/common_residue_predictions.csv.gz"
    inputs.add(prody_path)
    prody = read_csv(prody_path, dtype={"chain": str, "author_residue_id": str})
    keys = ["sample_id", "residue_index", "chain", "author_residue_id"]
    assert len(prody) == 961321 and not prody.duplicated(keys).any()
    joined = truth[keys + ["target_rmsf"]].merge(prody[keys + ["target_rmsf", "anm", "gnm"]], on=keys, how="left", suffixes=("_truth", "_saved"), validate="one_to_one")
    assert np.array_equal(joined.target_rmsf_truth, joined.target_rmsf_saved)
    assert np.isfinite(joined[["anm", "gnm"]].to_numpy()).all()
    for row in records:
        for method in ["ANM", "GNM"]:
            predictions.append(dict(sample_id=row["sample_id"], method=method, task_support="rmsf_relative", local_raw_available=True,
                                    path=str(prody_path), array_key=method.lower(), row_start="", row_stop="",
                                    rmsf_column="", hchi_column="", n_valid_rmsf=row["n_valid_rmsf"], n_valid_hchi=0,
                                    coverage_basis="this_run_exact_residue_key_join_target_equality_and_finite_checks",
                                    input_context="complete_structure_CA_network"))

    external_rows = []
    for method, sub in [("PLANET_full", "planet_full"), ("PLANET_sequence", "planet_seq"), ("RMSF_net_PDB", "rmsf_net")]:
        audit_path = BASE / sub / "verification_audit.json"
        meta_path = BASE / sub / "metadata.json"
        inputs.update([audit_path, meta_path])
        audit, meta = read_json(audit_path), read_json(meta_path)
        raw_count = len(list((BASE / sub).glob("*.npz")))
        assert audit["status"] == "complete" and audit["verified_records"] == 3063 and audit["verified_residues"] == 961321
        assert audit["failure_count"] == 0
        external_rows.append(dict(method=method, historical_status="complete", historically_verified_records=3063,
                                  historically_verified_positions=961321, expected_raw_items=int(audit["expected_items"]),
                                  observed_local_raw_npz=raw_count, local_per_residue_reusable=raw_count == int(audit["expected_items"]),
                                  local_raw_directory=str(BASE / sub), audit_path=str(audit_path), metadata_path=str(meta_path),
                                  per_record_metrics_path=str(BASE / "additional_original_oof/per_record_metrics.csv"),
                                  unit="Angstrom" if sub.startswith("planet") else "relative_score",
                                  input_context="independent_chain" if sub.startswith("planet") else "complete_structure_simulated_density",
                                  official_training_overlap="unresolved_not_audited", weights_sha256=meta["weights_sha256"],
                                  scientific_scope="SHP_is_not_Hchi" if sub.startswith("planet") else "PDB_only_simulated_density_adapter_not_absolute_RMSF"))
    for sub in ["additional_original_oof", "additional_family_oof"]:
        inputs.update(p for p in (BASE / sub).iterdir() if p.is_file())
    per_record = read_csv(BASE / "additional_original_oof/per_record_metrics.csv")
    assert len(per_record) == 3063 * 7
    assert all(set(q.sample_id) == sample_set for _, q in per_record.groupby("method"))
    assert np.all(per_record.common_fraction == 1)

    # Hash every local structure, isolated chain and prediction actually read.
    all_paths = sorted(inputs, key=str)
    def describe(path):
        return dict(path=str(path), bytes=path.stat().st_size, sha256=sha(path))
    file_rows = []
    with ThreadPoolExecutor(max_workers=6) as pool:
        for i, row in enumerate(pool.map(describe, all_paths), 1):
            file_rows.append(row)
            if i % 2000 == 0:
                print("Hashed", i, "of", len(all_paths), "existing files", flush=True)
    hashes = {x["path"]: x["sha256"] for x in file_rows}
    for row in records:
        assert hashes[row["structure_path"]] == row["structure_sha256_expected"]
        row["structure_sha256"] = hashes[row["structure_path"]]
        row["full_prediction_sha256"] = hashes[row["full_prediction_path"]]
        row["intrachain_prediction_sha256"] = hashes[row["intrachain_prediction_path"]]
    for row in chain_rows:
        assert hashes[row["structure_path"]] == row["structure_sha256_expected"]
        row["structure_sha256"] = hashes[row["structure_path"]]
    for row in predictions:
        row["file_sha256"] = hashes[row["path"]]
    write_csv("records.csv", records)
    write_csv("chains.csv", chain_rows)
    write_csv("existing_prediction_manifest.csv", predictions)
    write_csv("external_official_result_inventory.csv", external_rows)
    write_csv("annotation_architecture_candidates.csv", candidate_rows)
    write_csv("source_file_hashes.csv", file_rows)
    bootstrap = [{key: r[key] for key in ["sample_id", "dataset", "biological_system", "original_fold", "family_fold", "family_component_id", "biological_cluster_id", "main_group_candidate"]} for r in records]
    write_csv("bootstrap_membership.csv", bootstrap)
    groups = []
    selections = {"all_records": lambda r: True, "single_chain": lambda r: r["main_group_candidate"] == "single_chain",
                  "other_multichain": lambda r: r["main_group_candidate"] == "other_multichain",
                  "tcr_pmhc_itrace": lambda r: r["main_group_candidate"] == "tcr_pmhc_itrace",
                  "all_multichain_overlapping_subset": lambda r: r["is_multichain"],
                  "explicit_high_confidence_ppi_or_tcr_overlapping_subset": lambda r: r["explicit_high_confidence_ppi_or_tcr"]}
    for group, select in selections.items():
        for source in ["all", "atlas", "dynarepo", "itrace"]:
            items = [r for r in records if select(r) and (source == "all" or r["dataset"] == source)]
            groups.append(dict(group=group, dataset=source, records=len(items), chains=sum(x["n_chains"] for x in items),
                               rmsf_positions=sum(x["n_valid_rmsf"] for x in items), hchi_positions=sum(x["n_valid_hchi"] for x in items),
                               families=len({x["family_component_id"] for x in items}), biological_systems=len({x["biological_system"] for x in items})))
    write_csv("group_support.csv", groups)
    protocol = dict(status="assumption_pending_optional_user_preference", primary_groups=["single_chain", "other_multichain", "tcr_pmhc_itrace"],
                    primary_partition_exhaustive_and_disjoint=True, all_multichain_is_union=["other_multichain", "tcr_pmhc_itrace"],
                    other_multichain_definition="All observed multi-chain inputs outside ITrace, including biologically supported complexes; not a non-complex class.",
                    biological_architecture="No new homology threshold. Exact sequence differences do not establish heteromer status; fragmented/mutated homomers can differ. Cached source author-chain/entity matching is diagnostic only until chain renaming is resolved.",
                    explicit_complex_subset="Original high_confidence_ppi annotation OR ITrace TCR-pMHC source membership; overlaps all_multichain.",
                    statistical_unit="record; same residue support and same record identities for every paired method comparison",
                    legacy_uncertainty="Existing manuscript AP uses source-stratified paired record bootstrap. Its intervals are not independent-family intervals.",
                    family_robustness="Retain family_component_id for paired family-cluster resampling as a separate sensitivity estimand; do not merge original and family-trained Full/Intra predictions.",
                    itrace_family_limitation="All 244 ITrace records belong to one connected family and family fold0. Independent-family uncertainty within this subset is not identifiable.",
                    official_training_overlap="Published comparator weights were not retrained on FlexPPI folds. Training overlap, particularly PLANET ATLAS/mdCATH, is unresolved; neither record/family bootstrap nor local coverage verification establishes absence of training overlap.",
                    legacy_artifacts="E:/Project/FlexPPI/artifacts/evaluation/additional_baselines is a distinct old 276+244 cohort using BackFlip1.0 and older targets/estimands; excluded from current primary reuse pending explicit identity/label/input-hash matching.")
    (ROOT / "grouping_protocol.json").write_text(json.dumps(protocol, ensure_ascii=False, indent=2), encoding="utf-8")
    itrace = [r for r in records if r["is_tcr_pmhc_itrace"]]
    assert len({r["family_component_id"] for r in itrace}) == 1 and {r["family_fold"] for r in itrace} == {0}
    audit = dict(status="complete_local_inventory_and_identity_checks_grouping_pending_optional_preference",
                 created_utc=datetime.now(timezone.utc).isoformat(), records=3063, chains=4973, rmsf_support=961321, hchi_support=824988,
                 primary_group_records=dict(hist), family_components=1662, itrace_family_components=1,
                 original_source_hashes_unchanged=True, local_pdb_files_hash_verified=3063+4973,
                 full_intrachain_npz_verified=6126, all_original_position_identity_and_valid_masks_preserved=True,
                 reusable_local_per_residue_methods=["Full", "Intra", "Static", "AA", "BackFlip2_2.1", "ANM", "GNM"],
                 official_external_raw_missing=[x["method"] for x in external_rows if not x["local_per_residue_reusable"]],
                 declared_counts_and_current_hashes_verified_from=str(COHORT / "cohort_freeze.json"),
                 parent_freeze_sha256=sha(COHORT / "cohort_freeze.json"), reference_table_sha256=sha(COHORT / "truth_residues.csv.gz"),
                 no_model_inference_started=True, no_network_access=True, original_files_modified=False,
                 manifest_sha256={p.name: sha(p) for p in ROOT.glob("*.csv")}, grouping_protocol_sha256=sha(ROOT / "grouping_protocol.json"))
    (ROOT / "inventory_audit.json").write_text(json.dumps(audit, ensure_ascii=False, indent=2), encoding="utf-8")
    (ROOT / "README.md").write_text("""# Local benchmark cohort inventory\n\nRead records.csv for complete-structure inputs and chains.csv for independent-chain inputs. All paths point to existing local immutable files; their current content hashes were verified without rewriting original paths or original manifests. reference_table_path in records.csv resolves the common residue truth and masks. aligned_array_start/stop use the verified original static training-schema order for BackFlip/Static/AA NPZ slicing. Existing raw Full/Intra NPZ files additionally contain residue identities, masks and references.\n\nThe candidate disjoint partition is single_chain / other_multichain / tcr_pmhc_itrace, pending optional user preference. other_multichain includes genuine complexes. all_multichain is an overlapping union, not a fourth disjoint group. annotation_architecture_candidates.csv is diagnostic only: source PDB author chains may have been renamed in the actual MD input, and no arbitrary sequence-identity threshold was introduced.\n\nANM/GNM raw values are retained in their common-residue CSV; join on sample_id/residue_index/chain/author_residue_id. PLANET/RMSF-net historical completion certificates and per-record scores are local, but their 3063-cohort per-position NPZ predictions are not in the local copied baseline directories. Do not treat a certificate or a per-record score as raw per-position data. Their current inventory is explicit in external_official_result_inventory.csv.\n\nbootstrap_membership.csv preserves both original fold and family fold plus one stable family component for every sample. All methods must remain paired on sample_id. All ITrace records share one connected family, so that subgroup cannot yield independent-family confidence intervals. Published comparator training overlap is unresolved. No inference, fitting, normalization, calibration, new bootstrap or original-data mutation was performed.\n""", encoding="utf-8")
    print(json.dumps(audit, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
