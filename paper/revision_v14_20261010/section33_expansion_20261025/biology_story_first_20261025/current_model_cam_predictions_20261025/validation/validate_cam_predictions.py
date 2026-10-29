"""Independent descriptive check of CaM predictions against deposited NMR S2.

No model is imported or run.  Hchi, Ca-alpha RMSF and NMR order parameters
are different observables; this script evaluates only inverse ordering as a
declared hypothesis.  A residue, not a methyl probe or a fitted error bar,
is the counting unit.  All output is restricted to this validation folder.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata


HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent.parent
PARTNERS = ["CaMKK", "eNOS", "nNOS", "smMLCK"]
KEY_SITES = [18, 39, 105, 109, 144]
SEED = 20261025
TINY_EXP = 1e-12
TINY_PRED = 1e-8  # numerical equality only, not a biological threshold


def sha256(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")


def scalar(value):
    if value is None:
        return None
    return float(value) if np.isfinite(value) else None


def correlation(a, b, rank=False):
    a, b = np.asarray(a, float), np.asarray(b, float)
    good = np.isfinite(a) & np.isfinite(b)
    a, b = a[good], b[good]
    if len(a) < 3 or np.ptp(a) < 1e-14 or np.ptp(b) < 1e-14:
        return None
    if rank:
        a, b = rankdata(a), rankdata(b)
    return scalar(np.corrcoef(a, b)[0, 1])


def effect_summary(frame, exp_column, pred_column, bootstrap):
    """Input frame has one row per residue; positive experimental change
    means less order (-delta S2), positive predicted change means more motion.
    Residue bootstrap is descriptive and does not create biological replicates.
    """
    good = np.isfinite(frame[exp_column]) & np.isfinite(frame[pred_column])
    values = frame.loc[good, [exp_column, pred_column]].to_numpy(float)
    x, y = values[:, 0], values[:, 1]
    nonzero = (np.abs(x) > TINY_EXP) & (np.abs(y) > TINY_PRED)
    concordant = (x[nonzero] * y[nonzero]) > 0
    result = {
        "n_residues": len(values),
        "spearman_inverse_order": correlation(x, y, rank=True),
        "pearson_inverse_order": correlation(x, y),
        "n_nonzero_directions": int(nonzero.sum()),
        "n_concordant_directions": int(concordant.sum()),
        "direction_fraction": scalar(concordant.mean()) if len(concordant) else None,
        "n_predicted_ties": int((np.abs(y) <= TINY_PRED).sum()),
        "n_experimental_ties": int((np.abs(x) <= TINY_EXP).sum()),
        "spearman_residue_bootstrap_95_percentile_interval": None,
        "direction_residue_bootstrap_95_percentile_interval": None,
    }
    if bootstrap and len(values) >= 3:
        rng = np.random.default_rng(SEED)
        rhos, fractions = [], []
        for _ in range(bootstrap):
            ii = rng.integers(0, len(values), len(values))
            xx, yy = x[ii], y[ii]
            rho = correlation(xx, yy, rank=True)
            if rho is not None:
                rhos.append(rho)
            nz = (np.abs(xx) > TINY_EXP) & (np.abs(yy) > TINY_PRED)
            if nz.any():
                fractions.append(float(np.mean(xx[nz] * yy[nz] > 0)))
        if rhos:
            result["spearman_residue_bootstrap_95_percentile_interval"] = np.quantile(rhos, [.025, .975]).tolist()
        if fractions:
            result["direction_residue_bootstrap_95_percentile_interval"] = np.quantile(fractions, [.025, .975]).tolist()
    return result


def prepare_experimental(output):
    methyl_path = SOURCE / "common_44_unambiguous_methyl_probes_wide.csv"
    nh_path = SOURCE / "common_87_NH_probes_wide.csv"
    methyl = pd.read_csv(methyl_path)
    nh = pd.read_csv(nh_path)
    required = ["seq_id", "amino_acid", "atom"] + [f"{p}_{suffix}" for p in PARTNERS for suffix in ["s2", "fit_error"]]
    for table, expected_n, label in [(methyl, 44, "methyl"), (nh, 87, "NH")]:
        assert set(required) <= set(table), (label, "missing columns")
        assert len(table) == expected_n, (label, len(table))
        assert not table.duplicated(["seq_id", "atom"]).any(), (label, "duplicate probe")
        assert table.groupby("seq_id").amino_acid.nunique().max() == 1
    long_rows, aggregate_rows = [], []
    for seq_id, group in methyl.groupby("seq_id", sort=True):
        aa = group.amino_acid.iloc[0]
        pair_delta = group.smMLCK_s2 - group.nNOS_s2
        ambiguous = bool((pair_delta > TINY_EXP).any() and (pair_delta < -TINY_EXP).any())
        for partner in PARTNERS:
            s2 = group[f"{partner}_s2"]
            long_rows.extend({
                "partner": partner, "seq_id": int(seq_id), "amino_acid": aa,
                "atom": row.atom, "s2": float(row[f"{partner}_s2"]),
                "fit_error": float(row[f"{partner}_fit_error"]),
            } for _, row in group.iterrows())
            aggregate_rows.append({
                "partner": partner, "seq_id": int(seq_id), "amino_acid": aa,
                "methyl_s2_mean": float(s2.mean()),
                "methyl_s2_median": float(s2.median()),
                "methyl_s2_min": float(s2.min()), "methyl_s2_max": float(s2.max()),
                "n_methyl_probes": len(group), "methyl_atoms": ";".join(group.atom),
                "methyl_s2_per_probe": ";".join(f"{atom}:{value:.6g}" for atom, value in zip(group.atom, s2)),
                "methyl_fit_error_per_probe": ";".join(f"{atom}:{value:.6g}" for atom, value in zip(group.atom, group[f"{partner}_fit_error"])),
                "nNOS_to_smMLCK_opposite_probe_directions": ambiguous,
                "has_chi_target": aa not in {"A", "G"},
            })
    aggregate = pd.DataFrame(aggregate_rows)
    assert len(aggregate) == 34 * 4
    nh_rows = []
    for _, row in nh.iterrows():
        for partner in PARTNERS:
            nh_rows.append({"partner": partner, "seq_id": int(row.seq_id),
                            "amino_acid": row.amino_acid,
                            "NH_s2": float(row[f"{partner}_s2"]),
                            "NH_fit_error": float(row[f"{partner}_fit_error"])})
    nh_long = pd.DataFrame(nh_rows)
    aggregate.to_csv(output / "experimental_methyl_residue_aggregate.csv", index=False)
    pd.DataFrame(long_rows).to_csv(output / "experimental_methyl_probes.csv", index=False)
    nh_long.to_csv(output / "experimental_NH_residues.csv", index=False)
    methyl.loc[methyl.seq_id.isin(aggregate.loc[aggregate.nNOS_to_smMLCK_opposite_probe_directions, "seq_id"])].to_csv(
        output / "experimental_opposite_methyl_directions_nNOS_smMLCK.csv", index=False)
    audit = {
        "methyl_probe_count": len(methyl), "methyl_residue_count": int(methyl.seq_id.nunique()),
        "methyl_residue_count_with_chi_target": int(aggregate.query("partner == 'nNOS' and has_chi_target").seq_id.nunique()),
        "NH_probe_count": len(nh),
        "opposite_methyl_direction_residues_nNOS_smMLCK": sorted(map(int, aggregate.loc[aggregate.nNOS_to_smMLCK_opposite_probe_directions, "seq_id"].unique())),
        "source_sha256": {str(methyl_path): sha256(methyl_path), str(nh_path): sha256(nh_path)},
        "key_sites_specified_before_predictions": KEY_SITES,
        "aggregation_rule": "Arithmetic mean S2 over the deposited common methyl probes of each residue; one residue gets one weight. Mean and median agree for one or two probes, so they are not independent sensitivity checks.",
        "probe_errors": "Deposited fitting errors retained separately. No independent-probe or replicate assumption and no propagated aggregate experimental uncertainty.",
    }
    write_json(output / "experimental_schema_audit.json", audit)
    return aggregate, nh_long, audit


def load_predictions(path):
    pred = pd.read_csv(path)
    required = {"variant", "partner", "input_mode", "seq_id", "amino_acid", "pred_rmsf", "pred_hchi"}
    missing = required - set(pred)
    if missing:
        raise ValueError(f"{path}: missing columns {missing}")
    if "fold" not in pred:
        pred["fold"] = "aggregate"
    pred["fold"] = pred.fold.astype(str)
    pred["seq_id"] = pd.to_numeric(pred.seq_id, errors="raise").astype(int)
    pred = pred[pred.partner.isin(PARTNERS)].copy()
    assert len(pred), "No recognised CaM partner rows"
    key = ["variant", "input_mode", "fold", "partner", "seq_id"]
    if pred.duplicated(key).any():
        raise ValueError(f"{path}: duplicate state/residue predictions; aggregate any model repeats explicitly before validation")
    for col in ["pred_rmsf", "pred_hchi"]:
        pred[col] = pd.to_numeric(pred[col], errors="coerce")
    # Do not impute absent residues, missing predictions or Ala/Gly Hchi.
    pred.loc[pred.amino_acid.isin(["A", "G"]), "pred_hchi"] = np.nan
    return pred


def group_name(values):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", "__".join(map(str, values)))


def pair_table(merged, left, right, exp_col, pred_col):
    left_table = merged.loc[merged.partner == left].copy()
    right_table = merged.loc[merged.partner == right].copy()
    result = left_table.merge(right_table, on=["seq_id", "amino_acid"], suffixes=(f"_{left}", f"_{right}"), validate="one_to_one")
    result["delta_experimental_s2"] = result[f"{exp_col}_{right}"] - result[f"{exp_col}_{left}"]
    result["negative_delta_experimental_s2"] = -result.delta_experimental_s2
    result["delta_prediction"] = result[f"{pred_col}_{right}"] - result[f"{pred_col}_{left}"]
    result["inverse_direction_concordant"] = result.negative_delta_experimental_s2 * result.delta_prediction > 0
    result["predicted_numerical_tie"] = result.delta_prediction.abs() <= TINY_PRED
    return result


def centered_summary(table, exp_column, pred_column):
    good = table[np.isfinite(table[exp_column]) & np.isfinite(table[pred_column])].copy()
    counts = good.groupby("seq_id").partner.nunique()
    good = good[good.seq_id.isin(counts[counts == 4].index)].copy()
    if not len(good):
        return {"n_residues": 0}, good
    good["negative_experimental_centered"] = -(good[exp_column] - good.groupby("seq_id")[exp_column].transform("mean"))
    good["prediction_centered"] = good[pred_column] - good.groupby("seq_id")[pred_column].transform("mean")
    return {
        "n_residues": int(good.seq_id.nunique()),
        "n_state_values": len(good),
        "pooled_residue_centered_spearman_descriptive": correlation(good.negative_experimental_centered, good.prediction_centered, rank=True),
        "pooled_residue_centered_pearson_descriptive": correlation(good.negative_experimental_centered, good.prediction_centered),
        "inference_limit": "Four partner values of one residue are dependent; no state-level p value is calculated. Residue centering removes constant per-residue and amino-acid baselines.",
    }, good


def compare_input(path, methyl, nh, output, bootstrap):
    pred = load_predictions(path)
    # Keep exported paths below the common Windows 260-character limit.
    output_stem = {"predictions_fivefold_mean": "core_mean", "predictions_by_fold": "core_folds"}.get(path.stem, path.stem)
    input_dir = output / output_stem
    input_dir.mkdir(parents=True, exist_ok=True)
    # Check sequence before any statistics, not after selective filtering.
    exp_sequence = pd.concat([methyl[["seq_id", "amino_acid"]], nh[["seq_id", "amino_acid"]]]).drop_duplicates()
    check = pred.merge(exp_sequence, on="seq_id", suffixes=("_model", "_experiment"))
    bad = check[check.amino_acid_model != check.amino_acid_experiment]
    if len(bad):
        bad.to_csv(input_dir / "ERROR_amino_acid_mismatch.csv", index=False)
        raise ValueError(f"Sequence mismatch in {path}; see ERROR_amino_acid_mismatch.csv")
    all_summary, key_output, state_descriptive = [], [], []
    pair_rows_by_group = {}
    for identifiers, group in pred.groupby(["variant", "input_mode", "fold"], sort=True):
        name = group_name(identifiers)
        group_dir = input_dir / name
        group_dir.mkdir(exist_ok=True)
        side = methyl.merge(group, on=["partner", "seq_id", "amino_acid"], how="left", validate="one_to_one")
        backbone = nh.merge(group, on=["partner", "seq_id", "amino_acid"], how="left", validate="one_to_one")
        side.to_csv(group_dir / "all_methyl_residues_matched.csv", index=False)
        backbone.to_csv(group_dir / "all_NH_residues_matched.csv", index=False)
        summary = {"variant": identifiers[0], "input_mode": identifiers[1], "fold": identifiers[2], "groups_are_not_independent_replicates": True}
        side_pair = pair_table(side, "nNOS", "smMLCK", "methyl_s2_mean", "pred_hchi")
        side_pair["opposite_methyl_probe_directions"] = side_pair.nNOS_to_smMLCK_opposite_probe_directions_nNOS
        side_pair["has_chi_target"] = side_pair.has_chi_target_nNOS
        eligible = side_pair[side_pair.has_chi_target].copy()
        unambiguous = eligible[~eligible.opposite_methyl_probe_directions].copy()
        summary["nNOS_to_smMLCK_methyl_mean_all_chi_residues"] = effect_summary(eligible, "negative_delta_experimental_s2", "delta_prediction", bootstrap)
        summary["nNOS_to_smMLCK_methyl_mean_probe_direction_unambiguous"] = effect_summary(unambiguous, "negative_delta_experimental_s2", "delta_prediction", bootstrap)
        summary["nNOS_to_smMLCK_opposite_probe_direction_residues"] = sorted(map(int, eligible.loc[eligible.opposite_methyl_probe_directions, "seq_id"]))
        side_pair.to_csv(group_dir / "nNOS_to_smMLCK_methyl_residue_deltas.csv", index=False)
        nh_pair = pair_table(backbone, "nNOS", "smMLCK", "NH_s2", "pred_rmsf")
        summary["nNOS_to_smMLCK_NH_against_RMSF_exploratory"] = effect_summary(nh_pair, "negative_delta_experimental_s2", "delta_prediction", bootstrap)
        nh_pair.to_csv(group_dir / "nNOS_to_smMLCK_NH_residue_deltas.csv", index=False)
        summary["four_state_methyl_residue_centered"], centered_side = centered_summary(side[side.has_chi_target], "methyl_s2_mean", "pred_hchi")
        summary["four_state_NH_residue_centered_exploratory"], centered_nh = centered_summary(backbone, "NH_s2", "pred_rmsf")
        centered_side.to_csv(group_dir / "four_state_methyl_residue_centered.csv", index=False)
        centered_nh.to_csv(group_dir / "four_state_NH_residue_centered.csv", index=False)
        pairs = []
        for left, right in itertools.combinations(PARTNERS, 2):
            pair = pair_table(side[side.has_chi_target], left, right, "methyl_s2_mean", "pred_hchi")
            pair["left_partner"], pair["right_partner"] = left, right
            pairs.append(pair)
        pd.concat(pairs, ignore_index=True).to_csv(group_dir / "all_six_partner_pair_methyl_residue_deltas.csv", index=False)
        # Raw state correlations are secondary; residue identity can dominate them.
        for partner in PARTNERS:
            state = side[(side.partner == partner) & side.has_chi_target]
            item = {"variant": identifiers[0], "input_mode": identifiers[1], "fold": identifiers[2], "partner": partner,
                    "n_residues": int((np.isfinite(state.methyl_s2_mean) & np.isfinite(state.pred_hchi)).sum()),
                    "raw_state_spearman_minus_S2_vs_Hchi": correlation(-state.methyl_s2_mean, state.pred_hchi, rank=True),
                    "limitation": "Different amino-acid and residue identities can drive raw-state correlation; this is not evidence of partner-response prediction."}
            state_descriptive.append(item)
        key_side = side_pair[side_pair.seq_id.isin(KEY_SITES)].copy()
        key_nh = nh_pair[nh_pair.seq_id.isin(KEY_SITES)].copy()
        key = key_side[["seq_id", "amino_acid", "methyl_s2_mean_nNOS", "methyl_s2_mean_smMLCK", "n_methyl_probes_nNOS", "methyl_atoms_nNOS", "methyl_s2_per_probe_nNOS", "methyl_s2_per_probe_smMLCK", "opposite_methyl_probe_directions", "pred_hchi_nNOS", "pred_hchi_smMLCK", "delta_experimental_s2", "delta_prediction", "inverse_direction_concordant", "predicted_numerical_tie"]].rename(columns={"delta_experimental_s2": "delta_methyl_S2", "delta_prediction": "delta_Hchi"})
        key = key.merge(key_nh[["seq_id", "NH_s2_nNOS", "NH_s2_smMLCK", "NH_fit_error_nNOS", "NH_fit_error_smMLCK", "pred_rmsf_nNOS", "pred_rmsf_smMLCK", "delta_experimental_s2", "delta_prediction"]].rename(columns={"delta_experimental_s2": "delta_NH_S2", "delta_prediction": "delta_RMSF"}), on="seq_id", validate="one_to_one")
        for column, item in zip(["variant", "input_mode", "fold"], identifiers):
            key[column] = item
        key["interpretation_limit"] = "S2 and Hchi/RMSF are different observables; directional concordance alone does not establish calibrated dynamics or entropy prediction. M109 is a small methyl-S2 change comparator."
        key.to_csv(group_dir / "predeclared_five_key_sites.csv", index=False)
        key_output.append(key)
        pair_rows_by_group[name] = eligible
        write_json(group_dir / "metrics.json", summary)
        all_summary.append(summary)
    pd.concat(key_output, ignore_index=True).to_csv(input_dir / "all_groups_predeclared_five_key_sites.csv", index=False)
    pd.DataFrame(state_descriptive).to_csv(input_dir / "raw_state_correlations_secondary.csv", index=False)
    # Paired comparison uses the identical measured residues for both groups.
    contrasts = []
    for left_name, right_name in itertools.combinations(pair_rows_by_group, 2):
        l, r = pair_rows_by_group[left_name], pair_rows_by_group[right_name]
        joint = l[["seq_id", "negative_delta_experimental_s2", "delta_prediction", "opposite_methyl_probe_directions"]].merge(
            r[["seq_id", "delta_prediction"]], on="seq_id", suffixes=("_left", "_right"), validate="one_to_one")
        joint = joint[~joint.opposite_methyl_probe_directions]
        valid = np.isfinite(joint.negative_delta_experimental_s2) & np.isfinite(joint.delta_prediction_left) & np.isfinite(joint.delta_prediction_right)
        joint = joint[valid]
        a = effect_summary(joint, "negative_delta_experimental_s2", "delta_prediction_left", 0)
        b = effect_summary(joint, "negative_delta_experimental_s2", "delta_prediction_right", 0)
        contrasts.append({"left": left_name, "right": right_name, "n_shared_residues": len(joint),
                          "left_spearman": a["spearman_inverse_order"], "right_spearman": b["spearman_inverse_order"],
                          "left_direction_fraction": a["direction_fraction"], "right_direction_fraction": b["direction_fraction"],
                          "comparison_limit": "Shared-residue descriptive comparison. Input mode/variant/fold are not experimental replicates."})
    write_json(input_dir / "all_metrics.json", all_summary)
    write_json(input_dir / "shared_residue_group_comparisons.json", contrasts)
    write_json(input_dir / "input_audit.json", {"path": str(path), "sha256": sha256(path), "n_rows": len(pred), "group_count": len(all_summary)})
    return {"input": str(path), "output": str(input_dir), "groups": all_summary}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--predictions", action="append", type=Path, default=[])
    parser.add_argument("--schema-check", action="store_true", help="Prepare/audit the fixed experimental tables without predictions")
    parser.add_argument("--bootstrap", type=int, default=0)
    parser.add_argument("--out", type=Path, default=HERE / "results")
    args = parser.parse_args()
    output = args.out.resolve()
    if not output.is_relative_to(HERE):
        raise ValueError("Output must stay inside this validation directory")
    output.mkdir(parents=True, exist_ok=True)
    methyl, nh, audit = prepare_experimental(output)
    reports = []
    if not args.schema_check and not args.predictions:
        parser.error("Supply --predictions or --schema-check")
    if len({path.stem for path in args.predictions}) != len(args.predictions):
        raise ValueError("Prediction filenames must have different stems to avoid overwriting output")
    for path in args.predictions:
        reports.append(compare_input(path.resolve(), methyl, nh, output, args.bootstrap))
    write_json(output / "validation_manifest.json", {
        "script_sha256": sha256(__file__), "seed": SEED, "bootstrap_draws": args.bootstrap,
        "experimental_audit": audit, "input_count": len(reports),
        "outputs": [{"input": r["input"], "output": r["output"]} for r in reports],
        "predeclared_primary_comparison": "Family-held-out model only; nNOS to smMLCK; one row per chi-bearing residue; inverse sign and Spearman of delta Hchi versus -delta mean methyl S2. Opposite methyl probe directions are explicitly excluded from directional primary and retained in all-residue sensitivity.",
        "physical_limit": "NMR directional order S2 is not torsion entropy Hchi or positional RMSF. No cross-observable numeric error, fitted S2 prediction, thermodynamic entropy claim, or formal NH-stability claim is computed.",
        "baseline_limit": "Residue-constant/AA-only baseline yields zero partner deltas and no defined Spearman; intrachain and partner_deleted controls are assessed on common residues. Five-fold homolog-exposed predictions are sensitivity analyses, not family-held-out confirmation.",
        "statistical_limit": "One protein and four targets. Residue bootstrap intervals are descriptive, not evidence of independent biological replication, and no per-methyl or per-state significance tests are performed.",
    })
    print(json.dumps({"output": str(output), "audit": audit, "prediction_inputs_processed": len(reports)}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
