"""Fixed-budget controls trained on the same 2,431 purged records as strict fits."""
from pathlib import Path
import json
import numpy as np
import joblib
from sklearn.ensemble import HistGradientBoostingRegressor
from threadpoolctl import threadpool_limits
import run_corrected_static as static
from static_feature_definitions import weights

def main():
    p=static.HERE
    allowed={json.loads(x)["sample_id"] for x in (p/"sequence_audit/purged_training_local.jsonl").read_text(encoding="utf-8").splitlines()}
    assert len(allowed)==2431
    data=dict(np.load(static.OUT/"training_features.npz"))
    schema=json.loads((static.OUT/"training_schema.json").read_text(encoding="utf-8"))
    external=dict(np.load(static.OUT/"external_features.npz"))
    ext_schema=json.loads((static.OUT/"external_schema.json").read_text(encoding="utf-8"))
    selected=[r["record_index"] for r in schema["records"] if r["sample_id"] in allowed]
    member=np.isin(data["record_index"],selected)
    names=schema["names"]
    x,y,valid=data["X"],data["target"],data["valid"]
    categorical=[name in schema["categorical"] for name in names]
    configuration=dict(loss="squared_error",learning_rate=.07,max_iter=120,max_leaf_nodes=31,max_depth=None,
        min_samples_leaf=50,l2_regularization=1.,max_bins=63,early_stopping=False,random_state=7)
    predictions={name:np.full(external["target"].shape,np.nan,np.float32) for name in ["strict_static_sequence_complex","strict_aa_only"]}
    with threadpool_limits(limits=3):
        for m,metric in enumerate(static.METRICS):
            train=member&valid[:,m]
            balance=weights(schema["records"],data["record_index"],train)
            model=HistGradientBoostingRegressor(**configuration,categorical_features=categorical)
            model.fit(x[train],y[train,m],sample_weight=balance)
            joblib.dump(model,static.OUT/f"strict_static_{metric}.joblib",compress=3)
            predictions["strict_static_sequence_complex"][:,m]=model.predict(external["X"])
            aa=x[:,names.index("aa_offset_+0")].astype(int)
            centers=np.full(22,np.average(y[train,m],weights=balance))
            ids=np.flatnonzero(train)
            for code in np.unique(aa[train]):
                pick=aa[ids]==code
                centers[code]=np.average(y[ids[pick],m],weights=balance[pick])
            predictions["strict_aa_only"][:,m]=centers[external["X"][:,names.index("aa_offset_+0")].astype(int)]
            print("Strict static completed",metric,int(train.sum()),flush=True)
    static.summarize(external,ext_schema,predictions,"strict_external")
    for name,prediction in predictions.items():np.savez_compressed(static.OUT/(name+"_external.npz"),prediction=prediction)
    audit={"status":"complete","training_records":2431,"external_records":32,"configuration":configuration,
        "strict_freeze_sha256":static.sha(p/"strict_external_payload/strict_external_freeze.json"),
        "feature_schema_sha256":static.sha(static.OUT/"training_schema.json"),"code_sha256":static.sha(__file__),
        "no_test_tuning":True,"independence":"One shared-family component in external data; no independent-family CI."}
    (static.OUT/"strict_run_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")

if __name__=="__main__":main()
