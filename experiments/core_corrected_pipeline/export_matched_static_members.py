"""Reuse frozen OOF controls on external inputs, matching each neural train fold.

The all3063 fit remains a separately labelled reference. Primary neural external
predictions average five separately trained models; its exact training-member
control averages the same five OOF static members, not the all3063 static fit.
"""
import json
import numpy as np
import pandas as pd
import joblib
from threadpoolctl import threadpool_limits
from run_corrected_static import HERE, OUT, METRICS, sha, metric_values
from static_feature_definitions import weights


def main():
    schema = json.loads((OUT / "training_schema.json").read_text())
    ext_schema = json.loads((OUT / "external_schema.json").read_text())
    with np.load(OUT / "training_features.npz") as z:
        data = {k: z[k] for k in ["X", "target", "valid", "record_index", "fold"]}
    with np.load(OUT / "external_features.npz") as z:
        external = {k: z[k] for k in ["X", "target", "valid"]}
    name = schema["names"].index("aa_offset_+0")
    aa = data["X"][:, name].astype(int)
    ext_aa = external["X"][:, name].astype(int)
    predictions = {k: np.full((5,) + external["target"].shape, np.nan, np.float32)
                   for k in ["matched_static_sequence_complex", "matched_aa_only"]}
    member_audit, metrics = [], []
    with threadpool_limits(limits=3):
        for fold in range(5):
            train_path = HERE / f"training_payload/splits/fold_{fold}_train.jsonl"
            train_rows = [json.loads(line) for line in train_path.read_text().splitlines()]
            declared = {r["sample_id"] for r in train_rows}
            actual = {r["sample_id"] for r in schema["records"] if r["fold"] != fold}
            assert actual == declared and len(actual) == len(train_rows)
            test_path = HERE / f"training_payload/splits/fold_{fold}_test.jsonl"
            heldout = {json.loads(line)["sample_id"] for line in test_path.read_text().splitlines()}
            assert not actual.intersection(heldout)
            model_hashes = {}
            for m, metric in enumerate(METRICS):
                path = OUT / f"static_fold{fold}_{metric}.joblib"
                model = joblib.load(path)
                assert model.max_iter == 120 and not model.early_stopping and model.random_state == 7
                predictions["matched_static_sequence_complex"][fold, :, m] = model.predict(external["X"])
                train = (data["fold"] != fold) & data["valid"][:, m]
                balance = weights(schema["records"], data["record_index"], train)
                y = data["target"][:, m]
                centers = np.full(22, np.average(y[train], weights=balance))
                ids = np.flatnonzero(train)
                for code in np.unique(aa[train]):
                    pick = aa[ids] == code
                    centers[code] = np.average(y[ids[pick]], weights=balance[pick])
                predictions["matched_aa_only"][fold, :, m] = centers[ext_aa]
                model_hashes[metric] = sha(path)
            member_audit.append(dict(fold=fold, train_records=len(actual), heldout_records=len(heldout),
                train_manifest_sha256=sha(train_path), test_manifest_sha256=sha(test_path),
                static_checkpoint_sha256=model_hashes, training_members=sorted(actual)))
            print("Matched external static member", fold, len(actual), flush=True)
    for name, members in predictions.items():
        assert np.isfinite(members[:, external["valid"]]).all()
        # Preserve each member and the separately declared equal mean.
        ensemble = members.astype(np.float64).mean(axis=0)
        np.savez_compressed(OUT / (name + "_members_external.npz"), prediction=members)
        np.savez_compressed(OUT / (name + "_ensemble_external.npz"), prediction=ensemble)
        for member, p in [(str(i), members[i]) for i in range(5)] + [("five_fold_mean", ensemble)]:
            for r in ext_schema["records"]:
                for m, metric in enumerate(METRICS):
                    ids = np.arange(r["start"], r["stop"])
                    ids = ids[external["valid"][ids, m]]
                    metrics.append(dict(sample_id=r["sample_id"], dataset=r["dataset"],
                        biological_system=r["biological_system"], variant=name, member=member, metric=metric,
                        **metric_values(external["target"][ids, m], p[ids, m])))
    per = pd.DataFrame(metrics)
    per.to_csv(OUT / "matched_external_per_record_metrics.csv", index=False)
    numeric = ["mae", "rmse", "bias", "spearman", "pearson"]
    per.groupby(["variant", "member", "metric"])[numeric].mean().to_csv(OUT / "matched_external_macro_metrics.csv")
    audit = dict(status="complete", new_fits=0, reused_fits=10, training_freeze_sha256=sha(HERE / "training_payload/training_freeze.json"),
        external_schema_sha256=sha(OUT / "external_schema.json"), code_sha256=sha(__file__), members=member_audit,
        protocol="External predictions from the exact five OOF train folds, then equal five-member mean. All3063 static fits remain separate.",
        seed_comparison="Neural seeds7/11/23 share fold0 training members; compare each to static member0, not the five-fold static mean.",
        strict_protocol="Purged2431 static fits are separately preserved; these five-fold outputs must not replace them.",
        external_independence="32 systems in one shared family component; no independent-family CI.")
    (OUT / "matched_external_membership_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    print(per[per.member == "five_fold_mean"].groupby(["variant", "metric"])[numeric].mean().to_string(), flush=True)


if __name__ == "__main__":
    main()
