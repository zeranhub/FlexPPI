"""Protein-equal analysis of the frozen five-fold FlexPPI predictions.

Identity parsing is deliberately separate from training parsing: protonated
histidines HSD/HSE/HSP are histidine for identity only. No label or model input
is changed. The complete parsed protein chain, including residues absent from
the frozen model cache, defines exact-sequence identity. This is not a homology
family definition. Every chain and every valid residue is retained.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
from flexppi.constants import THREE_TO_ONE

AA_NAMES = dict(THREE_TO_ONE, HSD="H", HSE="H", HSP="H")
PATH_ROOTS = [
    ("/ranze/repos/FlexPPI", str(ROOT)),
    ("/ranze/data/ATLAS_extracted", "E:/data/protein_md_datasets/ATLAS_extracted"),
    ("/ranze/data/dynarepo/dynarepo_non_nucleic_full", "D:/data/dynarepo_non_nucleic_full"),
    ("/ranze/data/ITrace", "E:/data/ITrace"),
]
METRICS = ["rmsf", "chi_entropy"]
MEASURES = ["mae_flexppi", "mae_aa", "mae_gain", "rho_flexppi", "rho_aa",
            "rho_residual", "rho_gain"]


def digest(value):
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def local_path(value):
    for old, new in PATH_ROOTS:
        if value.startswith(old):
            return Path(new + value[len(old):])
    return Path(value)


def parse_identity(path):
    """Return all known protein residues and a conservative whole-state hash.

    Only model one and accepted alternate locations are used, as in the model
    parser. State hashing includes non-water heavy atoms, including ligands.
    Chain labels and author residue numbers remain in the state hash, so this
    only collapses truly identical coordinate records, not equivalent folds.
    """
    chains = {}
    state_atoms = []
    seen_atoms = set()
    saw_model = False
    with Path(path).open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            kind = line[:6].strip()
            if kind == "MODEL":
                if saw_model:
                    break
                saw_model = True
                continue
            if kind == "ENDMDL":
                break
            if kind not in {"ATOM", "HETATM"} or len(line) < 54:
                continue
            if line[16:17] not in {" ", "A", "1"}:
                continue
            name = line[17:20].strip().upper()
            chain = line[21:22].strip() or "_"
            author_id = line[22:26].strip() + line[26:27].strip()
            atom = line[12:16].strip().upper()
            try:
                xyz = [float(line[a:b]) for a, b in [(30, 38), (38, 46), (46, 54)]]
            except ValueError:
                continue
            if name in AA_NAMES:
                residue = chains.setdefault(chain, {})
                if author_id in residue and residue[author_id]["aa"] != AA_NAMES[name]:
                    raise ValueError(f"Conflicting residue identity: {path} {chain} {author_id}")
                residue.setdefault(author_id, {"aa": AA_NAMES[name], "name": name})
            element = line[76:78].strip().upper() if len(line) >= 78 else ""
            atom_key = (chain, author_id, name, atom)
            if (name in {"HOH", "WAT", "SOL", "TIP", "TIP3"}
                    or element in {"H", "D"}
                    or (not element and atom.lstrip("0123456789").startswith(("H", "D")))
                    or atom_key in seen_atoms):
                continue
            seen_atoms.add(atom_key)
            canonical_name = AA_NAMES.get(name, name)
            state_atoms.append((chain, author_id, canonical_name, atom, *xyz))
    if not chains:
        raise ValueError(f"No recognized protein chains: {path}")
    sequences = {chain: "".join(r["aa"] for r in residues.values())
                 for chain, residues in chains.items()}
    coordinate_hash = digest(json.dumps(sorted(state_atoms), separators=(",", ":")))
    return chains, sequences, coordinate_hash


def build_chain_map(manifest, data):
    cache_keys = data.groupby(["sample_id", "chain"], sort=False)
    rows, audits = [], []
    for count, item in enumerate(manifest, 1):
        record = item["record"]
        path = local_path(record["structure_path"])
        chains, sequences, structure_hash = parse_identity(path)
        ids = {c: "protein_" + digest(seq)[:20] for c, seq in sequences.items()}
        composition_ids = sorted(ids.values())  # List preserves copy numbers.
        composition_key = "composition_" + digest("|".join(composition_ids))[:20]
        state_id = "state_" + digest(json.dumps([structure_hash, record.get("temperature")]))[:20]
        for chain, seq in sequences.items():
            key = (item["sample_id"], chain)
            subset = cache_keys.get_group(key) if key in cache_keys.indices else None
            cache_count = 0 if subset is None else len(subset)
            missing_histidines = sum(r["name"] in {"HSD", "HSE", "HSP"} for r in chains[chain].values())
            if subset is not None:
                for residue in subset[["author_residue_id", "amino_acid"]].itertuples(index=False):
                    identity = chains[chain].get(str(residue.author_residue_id))
                    if identity is None or identity["aa"] != residue.amino_acid:
                        raise ValueError(f"Cache/identity mismatch {key} {residue}")
            rows.append(dict(
                sample_id=item["sample_id"], chain=chain, protein_id=ids[chain],
                state_id=state_id, structure_state_hash=structure_hash,
                composition_key=composition_key, composition_protein_ids="|".join(composition_ids),
                chain_length=len(seq), full_chain_length=len(seq), sequence=seq,
                short_chain=len(seq) < 20, cached_residue_count=cache_count,
                identity_only_residue_count=len(seq) - cache_count,
                normalized_hsd_hse_hsp_count=missing_histidines,
                in_prediction_cache=subset is not None, source=item["dataset"], fold=int(item["fold"]),
                biological_system=record.get("metadata", {}).get("biological_system", ""),
                cluster_id=record.get("metadata", {}).get("cluster_id", ""),
                pdb_id=record.get("metadata", {}).get("pdb_id", ""),
                temperature=record.get("temperature"), window_ns=record.get("window_ns"),
                source_group=record.get("source_group", ""),
                label_protocol=record.get("metadata", {}).get("label_protocol", "unspecified_in_manifest"),
                target_chains="|".join(record.get("target_chains", [])),
                structure_path=str(path), structure_file_sha256=sha256_file(path),
                all_protein_chain_count=len(sequences),
            ))
        audits.append(dict(sample_id=item["sample_id"], state_id=state_id,
                           composition_key=composition_key, parsed_chains=len(chains),
                           parsed_residues=sum(map(len, sequences.values())),
                           cached_residues=int((data.sample_id == item["sample_id"]).sum())))
        if count % 500 == 0:
            print(f"Parsed identity {count}/{len(manifest)}", flush=True)
    chain_map = pd.DataFrame(rows)
    if chain_map.duplicated(["sample_id", "chain"]).any():
        raise ValueError("Duplicate chain-map keys")
    joined = data[["sample_id", "chain"]].drop_duplicates().merge(
        chain_map[["sample_id", "chain"]], how="left", indicator=True)
    if not joined._merge.eq("both").all():
        raise ValueError("Some cached chains are absent from identity map")
    return chain_map, pd.DataFrame(audits)


def safe_rho(prediction, target):
    p, y = np.asarray(prediction, dtype=float), np.asarray(target, dtype=float)
    if len(y) < 3:
        return np.nan, "fewer_than_3_valid_residues"
    if np.ptp(y) == 0:
        return np.nan, "constant_target"
    if np.ptp(p) == 0:
        return np.nan, "constant_prediction"
    return float(spearmanr(p, y).statistic), "defined"


def chain_metric_table(data, chain_map):
    """Compute metrics on identical task-specific masks; reusable after joins."""
    rows = []
    lookup = chain_map.set_index(["sample_id", "chain"])
    for (sample, chain), group in data.groupby(["sample_id", "chain"], sort=False):
        common = lookup.loc[(sample, chain)].to_dict()
        common.update(sample_id=sample, chain=chain)
        for metric in METRICS:
            part = group.loc[group["valid_" + metric]]
            y, p, aa = [part[prefix + metric].to_numpy(float)
                        for prefix in ["target_", "pred_", "aa_prior_"]]
            if not np.isfinite([y, p, aa]).all():
                raise ValueError(f"Non-finite valid observations: {sample}/{chain}/{metric}")
            rho_p, reason_p = safe_rho(p, y)
            rho_aa, reason_aa = safe_rho(aa, y)
            rho_res, reason_res = safe_rho(p - aa, y - aa)
            mae_p = float(np.abs(p - y).mean()) if len(y) else np.nan
            mae_aa = float(np.abs(aa - y).mean()) if len(y) else np.nan
            rows.append(dict(**common, metric=metric, valid_count=len(y),
                             mae_flexppi=mae_p, mae_aa=mae_aa, mae_gain=mae_aa - mae_p,
                             rho_flexppi=rho_p, rho_aa=rho_aa, rho_residual=rho_res,
                             rho_gain=rho_p - rho_aa, rho_flexppi_status=reason_p,
                             rho_aa_status=reason_aa, rho_residual_status=reason_res,
                             target_mean=float(y.mean()) if len(y) else np.nan))
    return pd.DataFrame(rows)


def aggregate_hierarchy(chain_table):
    # All arithmetic averages are explicit: copies -> records -> states -> protein.
    keys = ["sample_id", "protein_id", "state_id", "composition_key", "metric"]
    copied = chain_table.groupby(keys, sort=False, as_index=False)[MEASURES + ["valid_count"]].mean()
    copy_counts = chain_table.groupby(keys, sort=False).size().rename("n_chain_copies").reset_index()
    copied = copied.merge(copy_counts, on=keys, validate="one_to_one")
    state_keys = ["protein_id", "state_id", "metric"]
    states = copied.groupby(state_keys, sort=False, as_index=False)[MEASURES + ["valid_count"]].mean()
    states = states.merge(copied.groupby(state_keys, sort=False).sample_id.nunique().rename("n_records").reset_index(), on=state_keys)
    states = states.merge(copied.groupby(state_keys, sort=False).n_chain_copies.sum().rename("n_chain_observations").reset_index(), on=state_keys)
    proteins = states.groupby(["protein_id", "metric"], sort=False, as_index=False)[MEASURES + ["valid_count"]].mean()
    proteins = proteins.rename(columns={"valid_count": "mean_valid_residues_per_chain"})
    counts = states.groupby(["protein_id", "metric"], sort=False).agg(
        n_states=("state_id", "size"), n_records=("n_records", "sum"),
        n_chain_observations=("n_chain_observations", "sum"),
        n_states_with_mae=("mae_flexppi", "count"),
        n_states_with_rho_flexppi=("rho_flexppi", "count"),
        n_states_with_rho_aa=("rho_aa", "count"),
        n_states_with_rho_residual=("rho_residual", "count"),
    ).reset_index()
    proteins = proteins.merge(counts, on=["protein_id", "metric"], validate="one_to_one")
    metadata = chain_table.groupby("protein_id", sort=False).agg(
        full_chain_length=("full_chain_length", "first"), sequence=("sequence", "first"),
        sources=("source", lambda x: "|".join(sorted(set(x)))),
        folds=("fold", lambda x: "|".join(map(str, sorted(set(x))))),
        n_compositions=("composition_key", "nunique"),
        biological_systems=("biological_system", lambda x: "|".join(sorted(set(x)))),
        cluster_ids=("cluster_id", lambda x: "|".join(sorted(set(x)))),
    ).reset_index()
    proteins = proteins.merge(metadata, on="protein_id", validate="many_to_one")
    proteins["short_chain"] = proteins.full_chain_length < 20
    proteins["relative_mae_reduction"] = proteins.mae_gain / proteins.mae_aa.replace(0, np.nan)
    return copied, states, proteins


def bootstrap_mean(values, seed=20260919, repeats=5000):
    values = np.asarray(values, float)
    values = values[np.isfinite(values)]
    if not len(values):
        return np.nan, np.nan, np.nan
    rng = np.random.default_rng(seed)
    means = np.empty(repeats)
    for start in range(0, repeats, 100):
        n = min(100, repeats - start)
        means[start:start + n] = values[rng.integers(0, len(values), size=(n, len(values)))].mean(axis=1)
    return float(values.mean()), *map(float, np.quantile(means, [.025, .975]))


def summary_table(proteins, repeats):
    rows = []
    for metric in METRICS:
        table = proteins.loc[proteins.metric == metric]
        for stratum, group in [("all_chains", table), ("protein_chains_ge20", table.loc[~table.short_chain]),
                               ("short_peptides_lt20", table.loc[table.short_chain])]:
            for measure in MEASURES + ["relative_mae_reduction"]:
                vals = group[measure].to_numpy(float)
                mean, low, high = bootstrap_mean(vals, repeats=repeats)
                finite = vals[np.isfinite(vals)]
                rows.append(dict(metric=metric, stratum=stratum, measure=measure,
                                 n_proteins_total=len(group), n_proteins_defined=len(finite),
                                 n_proteins_undefined=len(group) - len(finite),
                                 protein_equal_mean=mean, ci_low=low, ci_high=high,
                                 median=float(np.median(finite)) if len(finite) else np.nan,
                                 positive_fraction=float((finite > 0).mean()) if len(finite) else np.nan))
    return pd.DataFrame(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--residues", type=Path, default=ROOT / "artifacts/evaluation_atlas_hchi_v1/analysis/oof_per_residue.csv.gz")
    parser.add_argument("--manifest", type=Path, default=ROOT / "artifacts/evaluation_atlas_hchi_v1/oof/prediction_manifest.jsonl")
    parser.add_argument("--output", type=Path, default=ROOT / "output/information_sources_20260919/frozen_analysis")
    parser.add_argument("--bootstrap-repeats", type=int, default=5000)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    manifest = [json.loads(line) for line in args.manifest.read_text().splitlines()]
    data = pd.read_csv(args.residues, dtype={"sample_id": str, "chain": str,
                                            "author_residue_id": str, "amino_acid": str}, keep_default_na=False)
    for metric in METRICS:
        col = "valid_" + metric
        if data[col].dtype != bool:
            data[col] = data[col].astype(str).map({"True": True, "False": False})
            if data[col].isna().any():
                raise ValueError(f"Invalid boolean validity column: {col}")
    if data.duplicated(["sample_id", "chain", "author_residue_id"]).any():
        raise ValueError("Duplicated cached residues")
    if set(data.sample_id) != {x["sample_id"] for x in manifest}:
        raise ValueError("Residue table and manifest sample inventory differs")
    chain_map, identity_audit = build_chain_map(manifest, data)
    chain_map.to_csv(args.output / "chain_map.csv", index=False)
    identity_audit.to_csv(args.output / "identity_audit.csv", index=False)
    print(f"Identity complete: {len(chain_map)} chains, {chain_map.protein_id.nunique()} identities", flush=True)
    chain_table = chain_metric_table(data, chain_map)
    copied, states, proteins = aggregate_hierarchy(chain_table)
    summary = summary_table(proteins, args.bootstrap_repeats)
    for name, frame in [("chain_metrics", chain_table), ("record_protein_metrics", copied),
                        ("state_metrics", states), ("protein_metrics", proteins),
                        ("summary_statistics", summary)]:
        frame.to_csv(args.output / (name + ".csv"), index=False)
    rho_status = []
    for metric in METRICS:
        for measure in ["rho_flexppi", "rho_aa", "rho_residual"]:
            for status, count in chain_table.loc[chain_table.metric == metric, measure + "_status"].value_counts().items():
                rho_status.append(dict(metric=metric, measure=measure, status=status, n_chain_records=int(count)))
    pd.DataFrame(rho_status).to_csv(args.output / "rho_eligibility.csv", index=False)
    audit = dict(
        status="complete", input_residue_rows=len(data), input_records=len(manifest),
        parsed_chain_records=len(chain_map), cached_chain_records=int(chain_map.in_prediction_cache.sum()),
        unique_protein_sequences=int(chain_map.protein_id.nunique()),
        unique_short_peptide_sequences=int(chain_map.loc[chain_map.short_chain].protein_id.nunique()),
        unique_states=int(chain_map.state_id.nunique()),
        unique_compositions=int(chain_map.composition_key.nunique()),
        valid_residues={m: int(data["valid_" + m].sum()) for m in METRICS},
        summed_chain_valid_residues={m: int(chain_table.loc[chain_table.metric == m, "valid_count"].sum()) for m in METRICS},
        identity_only_residue_count=int(chain_map.identity_only_residue_count.sum()),
        parsed_histidine_variant_count=int(chain_map.normalized_hsd_hse_hsp_count.sum()),
        bootstrap=dict(unit="exact-sequence protein identity", repeats=args.bootstrap_repeats, seed=20260919,
                       scope="conditional descriptive CI; homologous proteins may remain dependent"),
        hierarchy="chain metrics -> same-sequence copies within record -> records sharing exact coordinate state and temperature -> states per protein -> proteins equal weight",
        interpretation="Full-model gain over a training-fold amino-acid mean combines all context and capacity effects. It does not isolate structural information.",
        identity_definition="SHA256 of complete first-model parsed chain sequence; HSD/HSE/HSP normalized to H for identity only; short chains retained; not homology families.",
        state_definition="SHA256 of first-model accepted-altloc non-water heavy-atom names, chain/author identifiers, residue identities and coordinates, then temperature; ligands included. Conservative exact-coordinate identity, not conformational clustering.",
        composition_definition="Sorted list of sequence identities for every parsed protein chain, with duplicate copies and short peptides retained; non-protein ligands belong to state_id, not protein composition_key.",
        rho_definition="Spearman within each chain/task; undefined for fewer than 3 valid residues or constant input. Residual rho correlates (prediction - training-fold AA prior) and (target - same prior). Hierarchical means exclude undefined rho only; eligibility counts are exported.",
        frozen_label_coverage="No label regeneration. Existing ILE/HSD and other masks remain exactly as in frozen OOF input.",
        file_sha256={str(args.residues): sha256_file(args.residues), str(args.manifest): sha256_file(args.manifest),
                     str(Path(__file__)): sha256_file(__file__)},
    )
    if audit["valid_residues"] != audit["summed_chain_valid_residues"]:
        raise ValueError("Valid residue coverage was not preserved")
    (args.output / "run_summary.json").write_text(json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(audit, indent=2, ensure_ascii=False), flush=True)
    print(summary.loc[(summary.stratum == "protein_chains_ge20") & summary.measure.isin(["mae_gain", "rho_residual"])].to_string(index=False))


if __name__ == "__main__":
    main()
