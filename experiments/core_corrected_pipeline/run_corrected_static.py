"""Corrected-input OOF and external controls; train-only fixed boosting budget."""
from pathlib import Path
import os
os.environ.setdefault("OMP_NUM_THREADS", "3")
os.environ.setdefault("MKL_NUM_THREADS", "3")
import sys
import json
import argparse
import hashlib
import time
import numpy as np
import pandas as pd
import joblib
from sklearn.ensemble import HistGradientBoostingRegressor
from scipy.stats import spearmanr, pearsonr
from threadpoolctl import threadpool_limits

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE/"training_payload/src"))
from flexppi.data.structures import parse_pdb, target_mask
from static_feature_definitions import features, weights
OUT = HERE/"static_corrected"
METRICS = ["rmsf", "chi_entropy"]

def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def metric_values(y, p):
    result = {"count":len(y), "mae":float(np.mean(abs(y-p))),
              "rmse":float(np.sqrt(np.mean((y-p)**2))), "bias":float(np.mean(p-y)),
              "spearman":np.nan, "pearson":np.nan}
    if len(y)>=4 and np.ptp(y)>0 and np.ptp(p)>0:
        result.update(spearman=float(spearmanr(y,p).statistic), pearson=float(pearsonr(y,p).statistic))
    return result

def corrected_records():
    original = HERE/"training_payload/splits/all_ready.jsonl"
    if not original.exists():
        raise RuntimeError("Training input freeze is required before static controls")
    rows=[json.loads(x) for x in original.read_text(encoding="utf-8").splitlines()]
    for r in rows:
        for key,directory in [("structure_path","structures"),("label_path","labels")]:
            r[key]=str(HERE/"training_payload/dataset"/directory/Path(r[key]).name)
    return rows

def extract_set(rows, stem):
    arrays, table, offset = {k:[] for k in ["X","target","valid","record_index","fold"]}, [], 0
    for index,r in enumerate(rows):
        structure=parse_pdb(r["structure_path"])
        x,names,sets,categorical=features(structure)
        assert x.shape==(structure.num_residues,47)
        with np.load(r["label_path"]) as z:
            y=z["values"][:,:2].copy()
            valid=z["mask"][:,:2] & structure.frame_mask[:,None] & target_mask(structure,r["target_chains"])[:,None]
        assert y.shape==valid.shape==(len(x),2)
        assert np.isfinite(y[valid]).all()
        n=len(x)
        arrays["X"].append(x)
        arrays["target"].append(y)
        arrays["valid"].append(valid)
        arrays["record_index"].append(np.full(n,index,np.int32))
        arrays["fold"].append(np.full(n,r["metadata"].get("fold",-1),np.int8))
        table.append(dict(record_index=index,sample_id=r["sample_id"],dataset=r["dataset"],
            biological_system=r["metadata"].get("biological_system",r["source_group"]),
            fold=r["metadata"].get("fold",-1),start=offset,stop=offset+n,
            input_sha256=sha(r["structure_path"]),label_sha256=sha(r["label_path"])))
        offset+=n
        if (index+1)%250==0: print("Features",stem,index+1,len(rows),flush=True)
    arrays={k:np.concatenate(v) for k,v in arrays.items()}
    np.savez_compressed(OUT/(stem+"_features.npz"),**arrays)
    schema=dict(records=table,names=names,sets=sets,categorical=categorical,
                feature_source_sha256=sha(HERE/"static_feature_definitions.py"),valid_counts=arrays["valid"].sum(0).tolist())
    (OUT/(stem+"_schema.json")).write_text(json.dumps(schema,indent=2),encoding="utf-8")
    return arrays,schema

def load_set(rows,stem):
    data=OUT/(stem+"_features.npz")
    metadata=OUT/(stem+"_schema.json")
    if not data.exists() or not metadata.exists(): return extract_set(rows,stem)
    schema=json.loads(metadata.read_text(encoding="utf-8"))
    assert [r["sample_id"] for r in rows]==[r["sample_id"] for r in schema["records"]]
    assert schema["feature_source_sha256"]==sha(HERE/"static_feature_definitions.py")
    for r,entry in zip(rows,schema["records"]):
        assert sha(r["structure_path"])==entry["input_sha256"]
        assert sha(r["label_path"])==entry["label_sha256"]
    return dict(np.load(data)),schema

def summarize(data,schema,predictions,role):
    rows=[]
    for variant,p in predictions.items():
        assert np.isfinite(p[data["valid"]]).all()
        for r in schema["records"]:
            for m,metric in enumerate(METRICS):
                ids=np.arange(r["start"],r["stop"])
                ids=ids[data["valid"][ids,m]]
                if len(ids): rows.append(dict(sample_id=r["sample_id"],dataset=r["dataset"],fold=r["fold"],
                    variant=variant,metric=metric,role=role,**metric_values(data["target"][ids,m],p[ids,m])))
    df=pd.DataFrame(rows)
    df.to_csv(OUT/(role+"_per_record_metrics.csv"),index=False)
    df.groupby(["variant","dataset","metric"])[["mae","rmse","bias","spearman","pearson"]].mean().to_csv(OUT/(role+"_macro_metrics.csv"))

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--external",action="store_true")
    args=ap.parse_args()
    OUT.mkdir(exist_ok=True)
    data,schema=load_set(corrected_records(),"training")
    x,y,valid=data["X"],data["target"],data["valid"]
    names=schema["names"]
    categorical=[name in schema["categorical"] for name in names]
    configuration=dict(loss="squared_error",learning_rate=.07,max_iter=120,max_leaf_nodes=31,max_depth=None,
        min_samples_leaf=50,l2_regularization=1.,max_bins=63,early_stopping=False,random_state=7)
    predictions={name:np.full(y.shape,np.nan,np.float32) for name in ["static_sequence_complex","aa_only"]}
    audit={"status":"running","configuration":configuration,"no_test_tuning":True,
        "training_freeze_sha256":sha(HERE/"training_payload/training_freeze.json"),"jobs":[]}
    with threadpool_limits(limits=3):
        for fold in range(5):
            for m,metric in enumerate(METRICS):
                train=(data["fold"]!=fold)&valid[:,m]
                test=(data["fold"]==fold)&valid[:,m]
                balance=weights(schema["records"],data["record_index"],train)
                path=OUT/f"static_fold{fold}_{metric}.joblib"
                if path.exists(): model=joblib.load(path)
                else:
                    model=HistGradientBoostingRegressor(**configuration,categorical_features=categorical)
                    model.fit(x[train],y[train,m],sample_weight=balance)
                    joblib.dump(model,path,compress=3)
                predictions["static_sequence_complex"][test,m]=model.predict(x[test])
                aa=x[:,names.index("aa_offset_+0")].astype(int)
                centers=np.full(22,np.average(y[train,m],weights=balance))
                train_indices=np.flatnonzero(train)
                for code in np.unique(aa[train]):
                    pick=aa[train_indices]==code
                    centers[code]=np.average(y[train_indices[pick],m],weights=balance[pick])
                predictions["aa_only"][test,m]=centers[aa[test]]
                audit["jobs"].append(dict(fold=fold,metric=metric,train_count=int(train.sum()),test_count=int(test.sum()),
                    checkpoint_sha256=sha(path)))
                (OUT/"run_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")
                print("Static OOF",fold,metric,flush=True)
        for name,p in predictions.items(): np.savez_compressed(OUT/(name+"_oof.npz"),prediction=p)
        summarize(data,schema,predictions,"oof")
        if args.external:
            extraction=json.loads((HERE/"external_thermopcd/extraction_summary.json").read_text(encoding="utf-8"))
            # Disclose all extraction failures; score every completed fixed-cohort input.
            external_rows=[json.loads(r) for r in (HERE/"external_thermopcd/ready_test.jsonl").read_text(encoding="utf-8").splitlines()]
            external,ext_schema=load_set(external_rows,"external")
            ext_predictions={name:np.full(external["target"].shape,np.nan,np.float32) for name in predictions}
            for m,metric in enumerate(METRICS):
                train=valid[:,m]
                balance=weights(schema["records"],data["record_index"],train)
                path=OUT/f"all_training_static_{metric}.joblib"
                if path.exists(): model=joblib.load(path)
                else:
                    model=HistGradientBoostingRegressor(**configuration,categorical_features=categorical)
                    model.fit(x[train],y[train,m],sample_weight=balance)
                    joblib.dump(model,path,compress=3)
                ext_predictions["static_sequence_complex"][:,m]=model.predict(external["X"])
                aa=x[:,names.index("aa_offset_+0")].astype(int)
                centers=np.full(22,np.average(y[train,m],weights=balance))
                train_indices=np.flatnonzero(train)
                for code in np.unique(aa[train]):
                    pick=aa[train_indices]==code
                    centers[code]=np.average(y[train_indices[pick],m],weights=balance[pick])
                ext_predictions["aa_only"][:,m]=centers[external["X"][:,names.index("aa_offset_+0")].astype(int)]
            for name,p in ext_predictions.items(): np.savez_compressed(OUT/(name+"_external.npz"),prediction=p)
            summarize(external,ext_schema,ext_predictions,"external")
            audit["external_extraction"]=extraction
            audit["external_valid_counts"]=external["valid"].sum(0).tolist()
        audit.update(status="complete",completed_utc=time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),code_sha256=sha(__file__))
        (OUT/"run_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")

if __name__=="__main__":main()
