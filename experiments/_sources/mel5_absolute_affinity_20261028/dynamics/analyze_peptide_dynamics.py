"""Recompute MEL5 four-state MD labels and audit the complete mutation response.

Uses the archived family-holdout predictions and the exact label implementation
already validated for the v125 manuscript. All original inputs are read-only.
Replicates are independent simulations; no frame is treated as an independent
experimental replicate, and no pairing between replica numbers is assumed.
"""
from pathlib import Path
import importlib.util
import hashlib
import json
import os
import sys

for thread_var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ[thread_var] = "1"

OUT = Path(__file__).resolve().parent
EXP = OUT.parent.parent
ARCHIVE = EXP / "tcr_expanded_screen_v41_20261017"
FIG = Path(r"E:\Project\FlexPPI\paper\revision_v14_20261010\section33_expansion_20261025\model_first_story_search_20261025\revision_v124_figure4_compact_palette")
sys.path.insert(0, str(EXP / "label_sensitivity/runtime"))
import numpy as np
import pandas as pd

spec = importlib.util.spec_from_file_location("archived_label_recompute", ARCHIVE / "recompute_allprotein_pair_labels.py")
old = importlib.util.module_from_spec(spec)
spec.loader.exec_module(old)
old.STAGE = Path(os.environ["TEMP"]) / "flexppi_mel5_peptide_dynamics_20261028"
old.STAGE.mkdir(parents=True, exist_ok=True)

STATES = ("3hg1", "4jff", "4jfd", "4jfe")
CONTRASTS = (("3hg1", "4jff", "TCR maturation"), ("4jff", "4jfd", "Peptide G4A"), ("4jff", "4jfe", "Peptide I7A"))
LABEL_DIR = OUT / "recomputed_labels"
LABEL_DIR.mkdir(parents=True, exist_ok=True)
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rank_corr(x, y):
    x, y = np.asarray(x), np.asarray(y)
    if len(x) < 3 or np.ptp(x) == 0 or np.ptp(y) == 0:
        return np.nan
    return float(pd.Series(x).rank().corr(pd.Series(y).rank()))


def bool_value(value):
    return value if isinstance(value, (bool, np.bool_)) else str(value).lower() == "true"


def safe_json(data):
    if isinstance(data, dict):
        return {str(k): safe_json(v) for k, v in data.items()}
    if isinstance(data, (list, tuple)):
        return [safe_json(x) for x in data]
    if isinstance(data, (np.integer,)):
        return int(data)
    if isinstance(data, (np.floating, float)):
        return float(data) if np.isfinite(data) else None
    if isinstance(data, np.bool_):
        return bool(data)
    return data


def main():
    states, all_frames, provenance, cache_checks, all_labels = {}, [], [], [], []
    for pdb in STATES:
        label_csv = LABEL_DIR / f"{pdb}_all_residue_per_replica_labels.csv"
        label_npz = LABEL_DIR / f"{pdb}_all_residue_replica_labels.npz"
        audit_file = LABEL_DIR / f"{pdb}_audit.json"
        if label_csv.exists() and label_npz.exists() and audit_file.exists():
            a = json.loads(audit_file.read_text(encoding="utf-8"))
            ppath = old.SOURCE / f"family_predictions/full/itrace_{pdb}.npz"
            assert sha(ppath) == a["prediction_sha256"]
            for source in a["provenance"]:
                assert sha(source["topology"]) == source["topology_sha256"]
                assert sha(source["trajectory"]) == source["trajectory_sha256"]
            with np.load(ppath, allow_pickle=False) as z:
                pred = {k: z[k].copy() for k in z.files}
            with np.load(label_npz, allow_pickle=False) as z:
                means, variance = z["mean"].copy(), z["sample_variance"].copy()
            frame = pd.read_csv(label_csv, dtype={"chain": str, "author_residue_id": str})
            prov, frames, checks = a["provenance"], a["frames"], a["cache_checks"]
        else:
            frame, pred, means, variance, prov, frames, checks = old.complex_labels(pdb, LABEL_DIR)
            a = dict(provenance=prov, frames=frames, cache_checks=checks,
                     prediction_sha256=sha(old.SOURCE / f"family_predictions/full/itrace_{pdb}.npz"))
            audit_file.write_text(json.dumps(a, indent=2), encoding="utf-8")
        with np.load(label_npz, allow_pickle=False) as z:
            values = z["values"].copy()
        states[pdb] = dict(pred=pred, mean=means, variance=variance, values=values,
                           lookup={(str(c), str(i)): j for j, (c, i) in enumerate(zip(pred["chain"], pred["author_residue_id"]))})
        provenance.extend(prov)
        all_frames.extend(frames)
        cache_checks.extend(checks)
        all_labels.append(frame)
    pd.concat(all_labels).to_csv(OUT / "all_four_states_per_residue_per_replica_labels.csv", index=False)

    cdr_src = pd.read_csv(FIG / "source_CDR_three_contrasts.csv", dtype={"chain_before": str, "chain_after": str, "id_before": str, "id_after": str})
    iface_src = pd.read_csv(FIG / "source_interface_three_contrasts.csv", dtype={"chain_before": str, "chain_after": str, "id_before": str, "id_after": str})
    pep_src = pd.read_csv(FIG / "source_peptide_four_states.csv", dtype={"chain": str, "id": str})
    assert len(cdr_src) == 135 and len(iface_src) == 90 and len(pep_src) == 40

    common_tcr = {}
    maturation = cdr_src[cdr_src.pdb_before.eq("3hg1") & cdr_src.pdb_after.eq("4jff")]
    for r in maturation.itertuples():
        common_tcr[(str(r.chain_after), str(r.id_after))] = str(r.id_before)
    assert len(common_tcr) == 45

    def entry(a, b, scope, chain_a, id_a, chain_b, id_b, region, role):
        sa, sb = states[a], states[b]
        ia, ib = sa["lookup"][(chain_a, id_a)], sb["lookup"][(chain_b, id_b)]
        aa_a, aa_b = str(sa["pred"]["amino_acid"][ia]), str(sb["pred"]["amino_acid"][ib])
        row = dict(pdb_before=a, pdb_after=b, contrast=f"{b}-{a}", set=scope,
                   role=role, region=region, chain_before=chain_a, id_before=id_a,
                   chain_after=chain_b, id_after=id_b, aa_before=aa_a, aa_after=aa_b,
                   same_amino_acid=aa_a == aa_b,
                   common_position=id_a if a == "3hg1" or role not in ("tcr_alpha", "tcr_beta") else common_tcr[(chain_a, id_a)])
        for j, metric in enumerate(("rmsf", "hchi")):
            eligible = bool(sa["pred"]["valid"][ia, j] and sb["pred"]["valid"][ib, j])
            if metric == "hchi":
                eligible = eligible and aa_a == aa_b
            row[metric + "_eligible"] = eligible
            for side, state, idx in (("before", sa, ia), ("after", sb, ib)):
                row[f"{metric}_model_{side}"] = float(state["pred"]["mean"][idx, j]) if eligible else np.nan
                row[f"{metric}_MD_{side}"] = float(state["mean"][idx, j]) if eligible else np.nan
                row[f"{metric}_MD_SD_{side}"] = float(np.sqrt(state["variance"][idx, j])) if eligible else np.nan
                for rep in range(3):
                    row[f"{metric}_MD_rep{rep+1}_{side}"] = float(state["values"][rep, idx, j]) if eligible else np.nan
            x = sb["values"][:, ib, j]
            y = sa["values"][:, ia, j]
            cross_delta = x[:, None] - y[None, :]
            row[f"{metric}_MD_delta"] = row[f"{metric}_MD_after"] - row[f"{metric}_MD_before"]
            row[f"{metric}_model_delta"] = row[f"{metric}_model_after"] - row[f"{metric}_model_before"]
            row[f"{metric}_delta_SE"] = float(np.sqrt((sa["variance"][ia, j] + sb["variance"][ib, j]) / 3)) if eligible else np.nan
            row[f"{metric}_same_direction"] = bool(row[f"{metric}_MD_delta"] * row[f"{metric}_model_delta"] > 0) if eligible else False
            row[f"{metric}_all9_cross_replica_deltas_positive"] = bool(np.all(cross_delta > 0)) if eligible else False
            row[f"{metric}_all9_cross_replica_deltas_negative"] = bool(np.all(cross_delta < 0)) if eligible else False
            row[f"{metric}_abs_delta_over_2SE"] = bool(abs(row[f"{metric}_MD_delta"]) > 2 * row[f"{metric}_delta_SE"]) if eligible else False
        return row

    records = []
    for source, scope in ((cdr_src, "CDR"), (iface_src, "pMHC_interface")):
        for r in source.itertuples():
            row = entry(r.pdb_before, r.pdb_after, scope, r.chain_before, r.id_before,
                        r.chain_after, r.id_after, r.region_before, r.role)
            # Published figure values are a separate frozen source to be verified.
            for metric in ("rmsf", "hchi"):
                if row[metric + "_eligible"]:
                    for side in ("before", "after"):
                        np.testing.assert_allclose(row[f"{metric}_model_{side}"], getattr(r, f"{metric}_full_{side}"), rtol=0, atol=1e-12)
                        np.testing.assert_allclose(row[f"{metric}_MD_{side}"], getattr(r, f"{metric}_md_{side}"), rtol=0, atol=3e-6)
            records.append(row)
    for a, b, _ in CONTRASTS:
        pa = pep_src[pep_src.pdb.eq(a)].sort_values("id", key=lambda s: s.astype(int))
        pb = pep_src[pep_src.pdb.eq(b)].sort_values("id", key=lambda s: s.astype(int))
        assert len(pa) == len(pb) == 10
        for ra, rb in zip(pa.itertuples(), pb.itertuples()):
            row = entry(a, b, "peptide", str(ra.chain), str(ra.id), str(rb.chain), str(rb.id), "peptide", "peptide")
            for metric in ("rmsf", "hchi"):
                if row[metric + "_eligible"]:
                    for side, r in (("before", ra), ("after", rb)):
                        np.testing.assert_allclose(row[f"{metric}_model_{side}"], getattr(r, f"{metric}_full"), rtol=0, atol=1e-12)
                        np.testing.assert_allclose(row[f"{metric}_MD_{side}"], getattr(r, f"{metric}_md"), rtol=0, atol=3e-6)
            records.append(row)
    frame = pd.DataFrame(records)
    assert len(frame) == 255
    frame.to_csv(OUT / "complete_three_contrasts_per_position.csv", index=False)
    mutation = frame[frame.pdb_before.eq("4jff")]
    mutation.to_csv(OUT / "peptide_mutations_complete_per_position.csv", index=False)

    region_stats, region_reps, absolute_checks = [], [], []
    for a, b, name in CONTRASTS:
        f = frame[frame.pdb_before.eq(a) & frame.pdb_after.eq(b)]
        scopes = {
            "all45_CDR": f.set.eq("CDR"),
            "all10_peptide": f.set.eq("peptide"),
            "unmutated_peptide": f.set.eq("peptide") & f.same_amino_acid,
            "all30_pMHC_interface": f.set.eq("pMHC_interface"),
            "HLA_interface": f.set.eq("pMHC_interface") & f.role.eq("mhc"),
        }
        for region in ("CDR1α", "CDR2α", "CDR3α", "CDR1β", "CDR2β", "CDR3β"):
            scopes[region] = f.set.eq("CDR") & f.region.eq(region)
        for scope, mask in scopes.items():
            for metric in ("rmsf", "hchi"):
                q = f[mask & f[metric + "_eligible"]]
                if q.empty:
                    continue
                x, y = q[metric + "_MD_delta"].to_numpy(), q[metric + "_model_delta"].to_numpy()
                stable = q[metric + "_abs_delta_over_2SE"].to_numpy(bool)
                nonoverlap = (q[metric + "_all9_cross_replica_deltas_positive"] | q[metric + "_all9_cross_replica_deltas_negative"]).to_numpy(bool)
                replica_means = {}
                for side, state in (("before", a), ("after", b)):
                    rv = np.array([q[f"{metric}_MD_rep{rep}_{side}"].mean() for rep in (1, 2, 3)])
                    replica_means[side] = rv
                    for rep, value in enumerate(rv, 1):
                        region_reps.append(dict(pdb=state, contrast=f"{b}-{a}", contrast_name=name,
                                                scope=scope, metric=metric, n_positions=len(q),
                                                replica=rep, value=value))
                    absolute_checks.append(dict(pdb=state, contrast=f"{b}-{a}", scope=scope, metric=metric,
                                                n=len(q), model_mean=q[f"{metric}_model_{side}"].mean(),
                                                MD_mean=q[f"{metric}_MD_{side}"].mean(),
                                                absolute_MAE=np.mean(abs(q[f"{metric}_model_{side}"] - q[f"{metric}_MD_{side}"])),
                                                absolute_Spearman=rank_corr(q[f"{metric}_MD_{side}"], q[f"{metric}_model_{side}"])))
                region_sem = float(np.sqrt(replica_means["before"].var(ddof=1) / 3 + replica_means["after"].var(ddof=1) / 3))
                region_cross = replica_means["after"][:, None] - replica_means["before"][None, :]
                region_stats.append(dict(contrast=f"{b}-{a}", contrast_name=name, scope=scope, metric=metric,
                                         n=len(q), n_unmutated=int(q.same_amino_acid.sum()),
                                         MD_before=q[metric + "_MD_before"].mean(), MD_after=q[metric + "_MD_after"].mean(),
                                         model_before=q[metric + "_model_before"].mean(), model_after=q[metric + "_model_after"].mean(),
                                         MD_mean_delta=x.mean(), model_mean_delta=y.mean(), region_delta_SE=region_sem,
                                         region_MD_delta_over_2SE=abs(x.mean()) > 2 * region_sem,
                                         region_all9_cross_replica_deltas_same_sign=bool(np.all(region_cross > 0) or np.all(region_cross < 0)),
                                         delta_Spearman=rank_corr(x, y), delta_MAE=np.mean(abs(x-y)),
                                         zero_delta_MAE=np.mean(abs(x)), n_direction_agree=int((x*y > 0).sum()),
                                         direction_fraction=float(np.mean(x*y > 0)),
                                         n_abs_delta_over_2SE=int(stable.sum()),
                                         direction_fraction_abs_delta_over_2SE=float(np.mean(x[stable]*y[stable] > 0)) if stable.any() else np.nan,
                                         n_all9_cross_replica_same_sign=int(nonoverlap.sum()),
                                         direction_fraction_all9_cross_replica_same_sign=float(np.mean(x[nonoverlap]*y[nonoverlap] > 0)) if nonoverlap.any() else np.nan))
    stats = pd.DataFrame(region_stats)
    stats.to_csv(OUT / "complete_predefined_region_change_metrics.csv", index=False)
    stats[stats.contrast.ne("4jff-3hg1")].to_csv(OUT / "peptide_mutations_predefined_region_change_metrics.csv", index=False)
    pd.DataFrame(region_reps).to_csv(OUT / "predefined_region_means_per_independent_replica.csv", index=False)
    pd.DataFrame(absolute_checks).to_csv(OUT / "absolute_state_validation_metrics.csv", index=False)
    mapping = []
    reference_cdr_keys = {(str(r.chain_before), str(r.id_before)) for r in maturation.itertuples()}
    raw_maps = {
        (str(row.pdb), int(row.index), int(row.replica)): int(row.raw_MD_resid)
        for labels in all_labels for row in labels.itertuples(index=False)
    }
    for pdb, state in states.items():
        p = state["pred"]
        for idx, (chain, resid, aa) in enumerate(zip(p["chain"], p["author_residue_id"], p["amino_acid"])):
            common = common_tcr.get((str(chain), str(resid)), "") if pdb != "3hg1" else (str(resid) if (str(chain), str(resid)) in reference_cdr_keys else "")
            mapping.append(dict(pdb=pdb, index=idx, chain=chain, model_author_residue_id=resid,
                                aa=aa, common_MEL5_CDR_position=common,
                                raw_MD_residue_id_replica1=raw_maps[(pdb, idx, 1)],
                                raw_MD_residue_id_replica2=raw_maps[(pdb, idx, 2)],
                                raw_MD_residue_id_replica3=raw_maps[(pdb, idx, 3)]))
    pd.DataFrame(mapping).to_csv(OUT / "full_state_numbering_map.csv", index=False)
    audit = dict(status="passed", new_MD=False, new_model_inference=False,
                 runtime=dict(python=sys.version, executable=sys.executable,
                              numpy=np.__version__, pandas=pd.__version__, MDAnalysis=old.mda.__version__),
                 archived_label_implementation=str(ARCHIVE / "recompute_allprotein_pair_labels.py"),
                 archived_label_implementation_sha256=sha(ARCHIVE / "recompute_allprotein_pair_labels.py"),
                 prediction_model="full-composition-chain-family-holdout fold0; archived predictions used by v125",
                 checkpoint_sha256s=sorted(set(x["checkpoint_sha256"] for x in provenance)),
                 data_scope="4 states x 3 independent MD replicas x every saved frame 0-200ns",
                 RMSF="Every protein chain independently C-alpha fitted to its own first frame, double precision; identical to original targets",
                 Hchi="Fixed 3 equal 120-degree bins for every chi; joint-state Shannon entropy; no-chi and amino-acid-changed comparisons excluded",
                 source_mutation_state_sequences={pdb: "".join(states[pdb]["pred"]["amino_acid"][states[pdb]["pred"]["chain"] == "C"]) for pdb in STATES},
                 all_prediction_targets_and_label_sample_variances_reproduced_to=3e-6,
                 all_v124_figure_model_values_reproduced_to=1e-12,
                 uncertainty="3 independent trajectories per state, unpaired. SE of difference from per-state replica variance. abs(delta)>2SE is a descriptive filter, not a significance test. All9 cross-replica deltas are the 3x3 cross differences, not 9 independent replicates.",
                 sources=[dict(path=str(FIG / name), sha256=sha(FIG / name)) for name in ("source_CDR_three_contrasts.csv", "source_interface_three_contrasts.csv", "source_peptide_four_states.csv")],
                 provenance=provenance, frames=all_frames, label_cache_checks=cache_checks,
                 script_sha256=sha(__file__))
    (OUT / "analysis_audit.json").write_text(json.dumps(safe_json(audit), ensure_ascii=False, indent=2), encoding="utf-8")
    brief = stats[stats.contrast.ne("4jff-3hg1") & (stats.metric.eq("rmsf") | stats.scope.eq("HLA_interface"))]
    print(brief.to_string(index=False), flush=True)
    print("OUTPUT", OUT, flush=True)


if __name__ == "__main__":
    main()
