"""Independent, exhaustive numerical case audit. New outputs only.

No model-fit or trajectory recomputation is performed here. The preexisting
out-of-fold predictions and raw-trajectory verification tables are joined by
explicit residue keys. All possible qualifying pairs are retained; neither
truth differences nor model errors determine membership in the primary set.
"""
import os
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'
os.environ['OPENBLAS_NUM_THREADS'] = '1'
from pathlib import Path
import json
import hashlib
from itertools import combinations
import numpy as np
import pandas as pd

OUT = Path(__file__).resolve().parent
ROOT = OUT.parent
SAMPLES = ["A0064", "A00C9"]
METHODS = ["full", "intra", "aa_sasa_contacts", "static47_sasa"]


def rho(a, b):
    a, b = np.asarray(a, float), np.asarray(b, float)
    ok = np.isfinite(a) & np.isfinite(b)
    return float(np.corrcoef(pd.Series(a[ok]).rank(method="average"),pd.Series(b[ok]).rank(method="average"))[0,1]) if ok.sum() >= 3 and np.ptp(a[ok]) > 1e-10 and np.ptp(b[ok]) > 1e-10 else np.nan


def sha(p):
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1048576), b""):
            h.update(chunk)
    return h.hexdigest()


def direction_score(delta):
    a = np.asarray(delta, float)
    return np.where(a > 1e-7, 1., np.where(a < -1e-7, 0., .5))


def all_pairs(frame):
    output = []
    for (sid, chain, aa), g in frame.groupby(["sample_id", "chain", "amino_acid"], sort=True):
        if len(g) < 4 or aa == "P":
            continue
        g = g.copy()
        g["sasa_rank"] = g.sidechain_sasa_complex_A2.rank(method="average", pct=True)
        # distance=0 is the end itself; >=3 excludes exactly three residues.
        g = g[g.distance_to_chain_end >= 3]
        for a, b in combinations(g.to_dict("records"), 2):
            if abs(a["sidechain_sasa_complex_A2"] - b["sidechain_sasa_complex_A2"]) > 10:
                continue
            if abs(a["md_rmsf"] - b["md_rmsf"]) > .1:
                continue
            low, high = (a, b) if a["md_hchi"] <= b["md_hchi"] else (b, a)
            row = dict(sample_id=sid, chain=chain, amino_acid=aa,
                       aa_count_before_end_exclusion=len(frame[(frame.sample_id == sid) & (frame.chain == chain) & (frame.amino_acid == aa)]),
                       low_residue_index=low["residue_index"], high_residue_index=high["residue_index"],
                       low_residue=aa + str(low["author_residue_id"]), high_residue=aa + str(high["author_residue_id"]),
                       low_hchi=low["md_hchi"], high_hchi=high["md_hchi"], true_gap=high["md_hchi"]-low["md_hchi"],
                       low_sasa_A2=low["sidechain_sasa_complex_A2"], high_sasa_A2=high["sidechain_sasa_complex_A2"],
                       sasa_abs_gap_A2=abs(high["sidechain_sasa_complex_A2"]-low["sidechain_sasa_complex_A2"]),
                       sasa_rank_abs_gap=abs(high["sasa_rank"]-low["sasa_rank"]),
                       low_rmsf_A=low["md_rmsf"], high_rmsf_A=high["md_rmsf"],
                       rmsf_abs_gap_A=abs(high["md_rmsf"]-low["md_rmsf"]),
                       low_interface=low["interface_ca8"], high_interface=high["interface_ca8"],
                       minimum_end_distance=min(high["distance_to_chain_end"], low["distance_to_chain_end"]))
            for method in METHODS:
                lo, hi = low[f"pred_{method}_hchi"], high[f"pred_{method}_hchi"]
                row[f"{method}_low"] = lo
                row[f"{method}_high"] = hi
                row[f"{method}_gap"] = hi-lo
                row[f"{method}_gap_abs_error"] = abs((hi-lo)-row["true_gap"])
                row[f"{method}_site_mae"] = (abs(lo-row["low_hchi"])+abs(hi-row["high_hchi"]))/2
            output.append(row)
    return pd.DataFrame(output)


def physical_pairs(pairs):
    result = []
    for sid, p in pairs.groupby("sample_id"):
        location = ROOT / "verified_case_trajectories" / sid
        audit = json.loads((location/"audit.json").read_text())
        assert audit["status"] == "passed" and audit["duration_ns"] == 500 and audit["n_original_windows"] == 5
        s = pd.read_csv(location/"all_sidechains_occupancy.csv").set_index("residue_index")
        w = pd.read_csv(location/"all_sidechains_per_window.csv")
        assert not w.duplicated(["residue_index", "block"]).any()
        table = {k: w.pivot(index="residue_index", columns="block", values=k) for k in
                 ["hchi", "fine_entropy_mean_per_chi", "width_mean_per_chi", "shifted60_joint_entropy"]}
        for row in p.to_dict("records"):
            lo, hi = int(row["low_residue_index"]), int(row["high_residue_index"])
            assert abs(s.loc[lo, "recomputed_hchi"]-row["low_hchi"]) < 1e-5
            assert abs(s.loc[hi, "recomputed_hchi"]-row["high_hchi"]) < 1e-5
            row["n_chi"] = int(s.loc[lo, "n_chi"])
            assert s.loc[lo, "n_chi"] == s.loc[hi, "n_chi"]
            for metric, t in table.items():
                gaps = t.loc[hi].to_numpy()-t.loc[lo].to_numpy()
                row[f"{metric}_block_positive_count"] = int((gaps > 1e-7).sum())
                row[f"{metric}_block_tie_count"] = int((np.abs(gaps) <= 1e-7).sum())
                row[f"{metric}_block_gap_mean"] = float(gaps.mean())
                row[f"{metric}_block_gap_min"] = float(gaps.min())
                row[f"{metric}_block_gap_max"] = float(gaps.max())
                for b, gap in enumerate(gaps, 1):
                    row[f"{metric}_block_{b}_gap"] = gap
            for metric in ["fine_entropy_mean_per_chi", "width_mean_per_chi"]:
                row[f"{metric}_whole500ns_gap"] = s.loc[hi, metric]-s.loc[lo, metric]
            result.append(row)
    return pd.DataFrame(result)


def pair_summary(pairs):
    rows=[]
    for (sid, chain), part in pairs.groupby(["sample_id", "chain"]):
        for ends, p in [("exclude_3", part), ("previous_distance_gt3", part[part.minimum_end_distance > 3])]:
            for rank, p2 in [("no_sasa_rank_caliper", p), ("sasa_rank_gap_le025", p[p.sasa_rank_abs_gap <= .25])]:
                for name, sub in [("all", p2), ("at_least_one_interface", p2[p2.low_interface | p2.high_interface]),
                                  ("both_interface", p2[p2.low_interface & p2.high_interface]),
                                  ("neither_interface", p2[~p2.low_interface & ~p2.high_interface])]:
                    for gap, q in [("all_nonzero", sub[sub.true_gap > 1e-7]), ("gap_ge_0p4", sub[sub.true_gap >= .4])]:
                        if q.empty:
                            continue
                        nres = len(set(q.low_residue_index) | set(q.high_residue_index))
                        for method in METHODS:
                            rows.append(dict(sample_id=sid, chain=chain, terminal_rule=ends, sasa_rank_rule=rank,
                                             interface_subset=name, truth_gap_subset=gap, method=method, n_pairs=len(q),
                                             n_distinct_residues=nres, n_amino_acid_types=q.amino_acid.nunique(),
                                             direction_accuracy=direction_score(q[f"{method}_gap"]).mean(),
                                             pair_gap_mae=q[f"{method}_gap_abs_error"].mean(),
                                             pair_weighted_site_mae=q[f"{method}_site_mae"].mean(),
                                             true_gap_mean=q.true_gap.mean(),
                                             hchi_same_order_all5_fraction=(q.hchi_block_positive_count == 5).mean(),
                                             hchi_same_order_ge4_fraction=(q.hchi_block_positive_count >= 4).mean(),
                                             fine_entropy_whole500_same_order_fraction=(q.fine_entropy_mean_per_chi_whole500ns_gap > 0).mean(),
                                             width_whole500_same_order_fraction=(q.width_mean_per_chi_whole500ns_gap > 0).mean(),
                                             fine_entropy_same_order_ge4_fraction=(q.fine_entropy_mean_per_chi_block_positive_count >= 4).mean(),
                                             width_same_order_ge4_fraction=(q.width_mean_per_chi_block_positive_count >= 4).mean(),
                                             shifted60_mean_same_order_fraction=(q.shifted60_joint_entropy_block_gap_mean > 0).mean(),
                                             shifted60_same_order_ge4_fraction=(q.shifted60_joint_entropy_block_positive_count >= 4).mean()))
    return pd.DataFrame(rows)


def residue_summary(frame, pairs):
    rows=[]
    for (sid, chain), chainpart in frame.groupby(["sample_id", "chain"]):
        pp=pairs[(pairs.sample_id == sid) & (pairs.chain == chain)]
        ids = set(pp.low_residue_index) | set(pp.high_residue_index)
        for subset, p in [("all_valid_nonpro", chainpart), ("primary_matched_unique_residues", chainpart[chainpart.residue_index.isin(ids)])]:
            for intf, q in [("all", p), ("interface", p[p.interface_ca8]), ("noninterface", p[~p.interface_ca8])]:
                if not len(q):
                    continue
                goodtypes = set(q.groupby("amino_acid").size().loc[lambda s:s >= 4].index)
                t = q[q.amino_acid.isin(goodtypes)].copy()
                for method in METHODS:
                    pred = f"pred_{method}_hchi"
                    typed_rho = np.nan
                    if len(t):
                        y = t.md_hchi-t.groupby("amino_acid").md_hchi.transform("mean")
                        v = t[pred]-t.groupby("amino_acid")[pred].transform("mean")
                        typed_rho = rho(y, v)
                    rows.append(dict(sample_id=sid, chain=chain, subset=subset, interface_subset=intf, method=method,
                                     n_residues=len(q), n_residues_typed=len(t), n_typed_groups=len(goodtypes),
                                     mae=np.abs(q.md_hchi-q[pred]).mean(), spearman_rho=rho(q.md_hchi, q[pred]),
                                     aa_centered_spearman_rho=typed_rho))
    return pd.DataFrame(rows)


def aa_pair_summary(pairs):
    rows=[]
    for (sid,chain,aa),p in pairs.groupby(["sample_id","chain","amino_acid"]):
        p=p[p.true_gap>1e-7]
        if p.empty:
            continue
        for method in METHODS:
            rows.append(dict(sample_id=sid,chain=chain,amino_acid=aa,method=method,n_pairs=len(p),
                             n_distinct_residues=len(set(p.low_residue_index)|set(p.high_residue_index)),
                             direction_accuracy=direction_score(p[f"{method}_gap"]).mean(),
                             pair_gap_mae=p[f"{method}_gap_abs_error"].mean(),
                             n_gap_ge_0p4=int((p.true_gap>=.4).sum())))
    return pd.DataFrame(rows)


def leave_one_out_audit(pairs):
    # A sensitivity to a single highly reused residue or amino-acid group;
    # this is not cross-validation and does not refit predictions.
    rows=[]
    for (sid,chain),p in pairs.groupby(["sample_id","chain"]):
        p=p[p.true_gap>1e-7]
        residue_ids=sorted(set(p.low_residue_index)|set(p.high_residue_index))
        for unit, values in [("residue",residue_ids),("amino_acid",sorted(p.amino_acid.unique()))]:
            for value in values:
                if unit=="residue":
                    q=p[(p.low_residue_index!=value)&(p.high_residue_index!=value)]
                else:
                    q=p[p.amino_acid!=value]
                if q.empty:
                    continue
                for other in METHODS[1:]:
                    rows.append(dict(sample_id=sid,chain=chain,removed_unit=unit,removed_value=value,
                                     comparator=other,n_pairs_remaining=len(q),
                                     direction_accuracy_gain_full=direction_score(q.full_gap).mean()-direction_score(q[f"{other}_gap"]).mean(),
                                     pair_gap_mae_reduction_full=q[f"{other}_gap_abs_error"].mean()-q.full_gap_abs_error.mean(),
                                     pair_weighted_site_mae_reduction_full=q[f"{other}_site_mae"].mean()-q.full_site_mae.mean()))
    return pd.DataFrame(rows)


def main():
    print('Imports ready; reading residue feature subset',flush=True)
    featurefile=ROOT/"feature_baseline/features_all.parquet"
    predfile=ROOT/"feature_baseline/oof_predictions_all.parquet"
    cols=["sample_id","residue_index","chain","author_residue_id","amino_acid","valid_hchi","valid_rmsf",
          "sidechain_complete","md_hchi","md_rmsf","interface_ca8","sidechain_sasa_complex_A2","distance_to_chain_end"]
    frame=pd.read_parquet(featurefile,columns=cols,filters=[("sample_id","in",SAMPLES)])
    print(f'Features loaded: {len(frame)} rows',flush=True)
    pcols=["sample_id","residue_index"]+[f"pred_{m}_hchi" for m in METHODS]
    pred=pd.read_parquet(predfile,columns=pcols,filters=[("sample_id","in",SAMPLES)])
    print(f'Predictions loaded: {len(pred)} rows',flush=True)
    frame=frame.merge(pred,on=["sample_id","residue_index"],validate="one_to_one")
    frame=frame[frame.valid_hchi & frame.valid_rmsf & frame.sidechain_complete & (frame.amino_acid != "P")].copy()
    frame.to_csv(OUT/"all_valid_nonpro_residues.csv",index=False)
    pairs=physical_pairs(all_pairs(frame))
    assert not pairs.duplicated(["sample_id","chain","low_residue_index","high_residue_index"]).any()
    pairs.to_csv(OUT/"all_qualifying_pairs_with_physical_verification.csv",index=False)
    pairs[pairs.true_gap>=.4].to_csv(OUT/"all_large_gap_pairs_no_model_selection.csv",index=False)
    pair_summary(pairs).to_csv(OUT/"pair_summary.csv",index=False)
    residue_summary(frame,pairs).to_csv(OUT/"unique_residue_summary.csv",index=False)
    aa_pair_summary(pairs).to_csv(OUT/"per_amino_acid_pair_summary.csv",index=False)
    leave_one_out_audit(pairs).to_csv(OUT/"leave_one_residue_or_amino_acid_out.csv",index=False)
    contract=dict(samples=SAMPLES,primary_criteria="Same sample, same chain and amino acid, >=4 valid complete nonPro sidechains of that type before terminal exclusion; exactly first/last 3 excluded (distance>=3); absolute sidechain complex SASA difference <=10A2 and true C-alpha RMSF difference <=0.10A; all unordered pairs, no Hchi-gap or prediction selection.",
                  sensitivities="Previous end-distance>3 and within-AA SASA-percentile difference<=0.25 are both reported separately; interface masks and truth-gap>=0.4 are post-stratifications of the entire retained pair set.",
                  scope="No amino-acid substitutions, sequence designs, biological experiments or new model training. Proline excluded because the frozen Hchi label is constant zero and gives no ordering information. Directional ties in model count 0.5; true numerical ties excluded (threshold 1e-7).",
                  caution="Pairs share residues and are not independent observations. Five 100ns windows are segments of one 500ns trajectory, not independent replicates. No pairwise p-values or independence-based confidence intervals are used. Physical fine metrics are per-chi averages with four 5degree histogram origins. Shifted60 120degree joint bins are an alternative coarse definition, not an experimental truth.",
                  n_valid_nonpro_residues=len(frame), n_pairs=len(pairs),
                  verified_input_audits={sid:json.loads((ROOT/"verified_case_trajectories"/sid/"audit.json").read_text()) for sid in SAMPLES},
                  script_sha256=sha(Path(__file__)))
    (OUT/"audit_contract.json").write_text(json.dumps(contract,indent=2),encoding="utf-8")
    print(json.dumps({"n_pairs":len(pairs),"by_chain":pairs.groupby(["sample_id","chain"]).size().to_dict().__repr__(),"out":str(OUT)},indent=2))


if __name__ == "__main__":
    main()
