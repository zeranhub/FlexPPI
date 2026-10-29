"""Five-fold, information-restricted static probes; never a neural ablation.

All estimators use the same fixed boosting budget and frozen OOF label masks.
No source identifiers, simulation conditions or test labels are input features.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from pathlib import Path

os.environ.setdefault("OMP_NUM_THREADS", "4")
os.environ.setdefault("MKL_NUM_THREADS", "4")
import numpy as np
from scipy.spatial.distance import cdist

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
from flexppi.data.structures import parse_pdb
from export_v2_oof_predictions import load_cache

OUT = ROOT / "output/information_sources_20260919/static_probes"
MANIFEST = ROOT / "artifacts/evaluation_atlas_hchi_v1/oof/prediction_manifest.jsonl"
MAPPINGS = {
    "/ranze/repos/FlexPPI": str(ROOT),
    "/ranze/data/ATLAS_extracted": "E:/data/protein_md_datasets/ATLAS_extracted",
    "/ranze/data/dynarepo/dynarepo_non_nucleic_full": "D:/data/dynarepo_non_nucleic_full",
    "/ranze/data/ITrace": "E:/data/ITrace",
}
METRICS = ["rmsf", "chi_entropy"]


def digest(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")
    temporary.replace(path)


def local_path(path):
    for old, new in MAPPINGS.items():
        if path.startswith(old):
            return Path(new + path[len(old):])
    return Path(path)


def dihedral(a, b, c, d):
    b0, b1, b2 = a - b, c - b, d - c
    b1 = b1 / np.maximum(np.linalg.norm(b1, axis=-1, keepdims=True), 1e-8)
    v = b0 - (b0 * b1).sum(-1, keepdims=True) * b1
    w = b2 - (b2 * b1).sum(-1, keepdims=True) * b1
    return np.arctan2((np.cross(b1, v) * w).sum(-1), (v * w).sum(-1))


def features(s):
    n = s.num_residues
    names, values = [], []

    def add(name, v):
        names.append(name)
        values.append(np.asarray(v, np.float32))

    # Sequence features see the same parsed sequence as the frozen model.
    for delta in range(-4, 5):
        v = np.full(n, 21, np.float32)  # 21 is a chain-boundary token, not a residue.
        for chain in np.unique(s.chain_index):
            ids = np.flatnonzero(s.chain_index == chain)
            positions = np.arange(len(ids)) + delta
            ok = (positions >= 0) & (positions < len(ids))
            v[ids[ok]] = s.aatype[ids[positions[ok]]]
        add("aa_offset_%+d" % delta, v)
    seq_cat = names.copy()
    length = np.zeros(n)
    fraction = np.zeros(n)
    end_distance = np.zeros(n)
    composition = np.zeros((n, 20))
    for chain in np.unique(s.chain_index):
        ids = np.flatnonzero(s.chain_index == chain)
        length[ids] = len(ids)
        fraction[ids] = np.arange(len(ids)) / max(len(ids) - 1, 1)
        end_distance[ids] = np.minimum(np.arange(len(ids)), np.arange(len(ids))[::-1])
        composition[ids] = np.bincount(s.aatype[ids], minlength=21)[:20] / len(ids)
    add("log_chain_length", np.log1p(length))
    add("fraction_chain_position", fraction)
    add("log_distance_to_chain_end", np.log1p(end_distance))
    positional = names[len(seq_cat):].copy()
    for a in range(20):
        add("chain_aa_fraction_%02d" % a, composition[:, a])
    seq = names.copy()

    # All geometry comes from N/CA/C, never side-chain atoms or trajectories.
    ca = s.backbone[:, 1].astype(np.float64)
    dist = cdist(ca, ca)
    same = s.chain_index[:, None] == s.chain_index[None, :]
    sep = np.abs(s.residue_index[:, None] - s.residue_index[None, :])
    valid_pair = s.frame_mask[:, None] & s.frame_mask[None, :]
    before = len(names)
    for radius in (6, 8, 10, 12):
        add("intra_contacts_%dA" % radius, np.log1p(((dist < radius) & same & (sep > 2) & valid_pair).sum(1)))
    for radius in (8, 12):
        add("intra_longrange_contacts_%dA" % radius, np.log1p(((dist < radius) & same & (sep >= 12) & valid_pair).sum(1)))
    phi, psi = np.full(n, np.nan), np.full(n, np.nan)
    for chain in np.unique(s.chain_index):
        ids = np.flatnonzero(s.chain_index == chain)
        prev, curr = ids[:-1], ids[1:]
        # Do not connect residues across a missing segment in the input structure.
        bonded = (np.linalg.norm(s.backbone[prev, 2] - s.backbone[curr, 0], axis=1) < 2.0)
        bonded &= s.frame_mask[prev] & s.frame_mask[curr]
        prev, curr = prev[bonded], curr[bonded]
        phi[curr] = dihedral(s.backbone[prev, 2], s.backbone[curr, 0], s.backbone[curr, 1], s.backbone[curr, 2])
        psi[prev] = dihedral(s.backbone[prev, 0], s.backbone[prev, 1], s.backbone[prev, 2], s.backbone[curr, 0])
    for name, angles in (("phi", phi), ("psi", psi)):
        add(name + "_sin", np.sin(angles))
        add(name + "_cos", np.cos(angles))
    v, w = s.backbone[:, 0] - ca, s.backbone[:, 2] - ca
    cosine = (v * w).sum(1) / np.maximum(np.linalg.norm(v, axis=1) * np.linalg.norm(w, axis=1), 1e-8)
    cosine[~s.frame_mask] = np.nan
    add("N_CA_C_cosine", cosine)
    intra = names[before:].copy()
    before = len(names)
    for radius in (6, 8, 10, 12):
        add("partner_contacts_%dA" % radius, np.log1p(((dist < radius) & ~same & valid_pair).sum(1)))
    partner = names[before:].copy()
    sets = {
        "sequence": seq,
        "geometry": positional + intra + partner,
        "sequence_intra": seq + intra,
        "sequence_complex": seq + intra + partner,
    }
    return np.stack(values, axis=1), names, sets, seq_cat


def extract():
    OUT.mkdir(parents=True, exist_ok=True)
    entries = [json.loads(line) for line in MANIFEST.read_text().splitlines()]
    assert len(entries) == 3068
    fields = {k: [] for k in ("X", "target", "valid", "full", "fold", "record_index", "residue_index")}
    record_table, rows, start = [], 0, time.monotonic()
    for k, entry in enumerate(entries):
        arrays = load_cache(Path(entry["cache_path"]), entry["identity"])
        path = local_path(entry["record"]["structure_path"])
        if digest(path) != entry["identity"]["input_sha256"]:
            raise ValueError("Changed structure: " + str(path))
        s = parse_pdb(str(path))
        assert np.array_equal(arrays["aatype"], s.aatype)
        assert np.array_equal(arrays["chain"], s.chain_ids)
        assert np.array_equal(arrays["author_residue_id"], s.author_residue_ids)
        x, names, sets, cat = features(s)
        n = len(x)
        fields["X"].append(x)
        fields["target"].append(arrays["target"].astype(np.float32))
        fields["valid"].append(arrays["valid"])
        fields["full"].append(arrays["mean"].astype(np.float32))
        fields["fold"].append(np.full(n, entry["fold"], np.int8))
        fields["record_index"].append(np.full(n, k, np.int32))
        fields["residue_index"].append(np.arange(n, dtype=np.int32))
        record_table.append(dict(record_index=k, sample_id=entry["sample_id"], dataset=entry["dataset"],
                                 biological_system=entry["record"]["metadata"]["biological_system"],
                                 fold=entry["fold"], start=rows, stop=rows+n,
                                 structure_path=str(path), input_sha256=entry["identity"]["input_sha256"]))
        rows += n
        if (k + 1) % 250 == 0:
            print("Extracted %d/3068 structures; %.1fs" % (k + 1, time.monotonic() - start), flush=True)
    arrays = {k: np.concatenate(v) for k, v in fields.items()}
    assert arrays["valid"].sum(0).tolist() == [959347, 798385]
    np.savez_compressed(OUT / "features.npz", **arrays)
    write_json(OUT / "feature_schema.json", dict(names=names, feature_sets=sets, categorical=cat,
               records=record_table, rows=rows, valid_counts=arrays["valid"].sum(0).tolist(),
               manifest_sha256=digest(MANIFEST), script_sha256=digest(__file__),
               extraction_seconds=time.monotonic()-start))
    print("Features ready:", arrays["X"].shape, flush=True)


def weights(records, record_ids, train_mask):
    # Same expected source/system/record balance as the neural run; each metric
    # has uniform residues within a record. Unlike neural NLL, the loss is MSE.
    from collections import Counter
    counts = np.bincount(record_ids[train_mask], minlength=len(records))
    active = [r for r in records if counts[r["record_index"]] > 0]
    system_records = Counter((r["dataset"], r["biological_system"]) for r in active)
    source_systems = Counter(source for source, system in system_records)
    per_record = np.zeros(len(records))
    for r in active:
        idx = r["record_index"]
        per_record[idx] = 1 / (len(source_systems) * source_systems[r["dataset"]]
                               * system_records[(r["dataset"], r["biological_system"])] * counts[idx])
    w = per_record[record_ids[train_mask]]
    return w / w.mean()


def train(args):
    import joblib
    import sklearn
    from sklearn.ensemble import HistGradientBoostingRegressor
    from threadpoolctl import threadpool_limits
    data = np.load(OUT / "features.npz")
    schema = json.loads((OUT / "feature_schema.json").read_text())
    if digest(MANIFEST) != schema["manifest_sha256"]:
        raise ValueError("Manifest changed")
    names, records = schema["names"], schema["records"]
    config = dict(loss="squared_error", learning_rate=0.07, max_iter=args.iterations,
                  max_leaf_nodes=31, max_depth=None, min_samples_leaf=50,
                  l2_regularization=1.0, max_bins=63, early_stopping=False, random_state=7)
    run = dict(status="running", interpretation="Fixed-budget histogram boosting probes, NOT FlexPPI neural ablations",
               training_config=config, sklearn_version=sklearn.__version__, seed=7,
               threads=args.threads, source_features=False, metric_specific_frozen_masks=True,
               test_model_selection=False, feature_schema_sha256=digest(OUT / "feature_schema.json"),
               label_protocol="Frozen ATLAS-Hchi-v1 coverage; existing ILE/His limitations retained",
               weighting="Equal source, equal biological system in source, equal record in system, equal valid residue in record",
               jobs=[])
    jobs = OUT / "jobs"
    jobs.mkdir(exist_ok=True)
    selected = args.variants or list(schema["feature_sets"])
    allpred = {}
    start = time.monotonic()
    with threadpool_limits(limits=args.threads):
        for variant in selected:
            columns = [names.index(x) for x in schema["feature_sets"][variant]]
            categorical = [names[c] in schema["categorical"] for c in columns]
            x = data["X"][:, columns]
            pred = np.full(data["target"].shape, np.nan, np.float32)
            for fold in range(5):
                for metric_index, metric in enumerate(METRICS):
                    train_mask = (data["fold"] != fold) & data["valid"][:, metric_index]
                    test_mask = (data["fold"] == fold) & data["valid"][:, metric_index]
                    stem = "%s_fold%d_%s" % (variant, fold, metric)
                    npz = jobs / (stem + ".npz")
                    rowfile = jobs / (stem + ".json")
                    contract = dict(variant=variant, fold=fold, metric=metric, config=config,
                                    feature_schema_sha256=run["feature_schema_sha256"], columns=columns,
                                    train_count=int(train_mask.sum()), test_count=int(test_mask.sum()))
                    tick = time.monotonic()
                    if npz.exists() and rowfile.exists():
                        row = json.loads(rowfile.read_text())
                        assert row["contract"] == contract, "Incompatible cached job"
                        p = np.load(npz)["prediction"]
                    else:
                        model = HistGradientBoostingRegressor(**config, categorical_features=categorical)
                        model.fit(x[train_mask], data["target"][train_mask, metric_index],
                                  sample_weight=weights(records, data["record_index"], train_mask))
                        p = model.predict(x[test_mask]).astype(np.float32)
                        if not np.isfinite(p).all():
                            raise ValueError("Nonfinite predictions")
                        np.savez_compressed(npz, prediction=p, row_index=np.flatnonzero(test_mask))
                        joblib.dump(model, jobs / (stem + ".joblib"), compress=3)
                        row = dict(contract=contract, elapsed_seconds=time.monotonic()-tick,
                                   test_mae=float(np.mean(np.abs(p-data["target"][test_mask, metric_index]))),
                                   prediction_sha256=digest(npz))
                        write_json(rowfile, row)
                    pred[test_mask, metric_index] = p
                    run["jobs"].append(row)
                    run["elapsed_seconds"] = time.monotonic()-start
                    write_json(OUT / "run_summary.json", run)
                    print(stem, "n=", int(train_mask.sum()), "MAE=", row["test_mae"],
                          "seconds=", round(row["elapsed_seconds"], 1), flush=True)
            assert np.isfinite(pred[data["valid"]]).all()
            allpred[variant] = pred
            np.savez_compressed(OUT / (variant + "_oof.npz"), prediction=pred)
    run["status"] = "complete"
    run["variants"] = selected
    run["elapsed_seconds"] = time.monotonic()-start
    write_json(OUT / "run_summary.json", run)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("action", choices=["extract", "train"])
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--iterations", type=int, default=120)
    parser.add_argument("--variants", nargs="+")
    args = parser.parse_args()
    if args.action == "extract":
        extract()
    else:
        train(args)


if __name__ == "__main__":
    main()
