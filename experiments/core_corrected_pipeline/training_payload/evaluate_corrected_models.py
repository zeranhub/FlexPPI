"""Fixed-final-epoch OOF export; preserve per-record, source, fold and residue data."""
from pathlib import Path
import sys
import json
import time
import hashlib
import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr, pearsonr
ROOT = Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"src"))
from flexppi.config import ModelConfig
from flexppi.data.dataset import StateDataset, collate_states
from flexppi.data.structures import parse_pdb
from flexppi.model.information_ablation import build_information_ablation
from flexppi.trainer import move_to_device

def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(1024*1024),b""):h.update(b)
    return h.hexdigest()

def statistics(y,p):
    good=np.isfinite(y)&np.isfinite(p)
    y,p=y[good],p[good]
    result=dict(count=len(y),mae=float(np.mean(abs(y-p))),rmse=float(np.sqrt(np.mean((y-p)**2))),
                bias=float(np.mean(p-y)),spearman=np.nan,pearson=np.nan)
    if len(y)>=4 and np.ptp(y)>0 and np.ptp(p)>0:
        result.update(spearman=float(spearmanr(y,p).statistic),pearson=float(pearsonr(y,p).statistic))
    return result

def main():
    torch.set_num_threads(2)
    out=ROOT/"evaluation"
    out.mkdir(exist_ok=True)
    records=[]
    sources=[]
    manifest=[]
    for variant in ["full","intrachain"]:
        variant_dir=out/variant
        variant_dir.mkdir(exist_ok=True)
        for fold in range(5):
            path=ROOT/f"runs/{variant}/fold_{fold}/epoch_040.pt"
            checkpoint=torch.load(path,map_location="cpu",weights_only=False)
            assert checkpoint["epoch"]==40
            c=checkpoint["config"]
            assert c["model"]["information_ablation"]==variant and c["seed"]==7
            model=build_information_ablation(ModelConfig(**c["model"]),variant).to("cuda").eval()
            dataset=StateDataset(str(ROOT/f"splits/fold_{fold}_test.jsonl"),max_residues=1024)
            expected=[json.loads(x) for x in (ROOT/f"splits/fold_{fold}_test.jsonl").read_text().splitlines()]
            assert len(dataset)==len(expected)
            model.load_state_dict(checkpoint["model"],strict=True)
            checkpoint_hash=sha(path)
            for i,record in enumerate(dataset.records):
                ex=dataset[i]
                structure=parse_pdb(record.structure_path)
                batch=move_to_device(collate_states([ex]),torch.device("cuda"))
                with torch.inference_mode():
                    prediction=model.forward_state(batch)
                mean=prediction["mean"][0].cpu().numpy()
                variance=prediction["log_variance"][0].exp().cpu().numpy()
                target=ex["label_values"].numpy()
                valid=(ex["label_mask"] & ex["target_mask"][:,None] & prediction["mask"][0].cpu()[:,None]).numpy()
                assert np.isfinite(mean).all() and np.isfinite(variance).all()
                cache=variant_dir/(record.sample_id+".npz")
                np.savez_compressed(cache,mean=mean,model_variance=variance,target=target,valid=valid,
                    label_variance=ex["label_variances"].numpy(),amino_acid=np.asarray(list(structure.sequence)),
                    chain=np.asarray(structure.chain_ids),author_residue_id=np.asarray(structure.author_residue_ids),aatype=structure.aatype)
                manifest.append(dict(variant=variant,fold=fold,sample_id=record.sample_id,dataset=record.dataset,
                    cache_path=str(cache),checkpoint_sha256=checkpoint_hash,input_sha256=sha(record.structure_path),
                    label_sha256=sha(record.label_path),record=record.to_dict()))
                d={"sample_id":record.sample_id,"dataset":record.dataset,"fold":fold,"variant":variant,
                   "residue_index":np.arange(len(mean)),"amino_acid":list(structure.sequence),
                   "chain":structure.chain_ids,"author_residue_id":structure.author_residue_ids}
                for m,metric in enumerate(c["model"]["metrics"]):
                    if valid[:,m].any():
                        records.append(dict(sample_id=record.sample_id,dataset=record.dataset,fold=fold,variant=variant,
                            metric=metric,**statistics(target[valid[:,m],m],mean[valid[:,m],m])))
                    d.update({f"target_{metric}":target[:,m],f"pred_{metric}":mean[:,m],
                              f"valid_{metric}":valid[:,m],f"sd_{metric}":np.sqrt(variance[:,m])})
                sources.append(pd.DataFrame(d))
                if (i+1)%50==0:print(variant,fold,i+1,len(dataset),flush=True)
            del model
            torch.cuda.empty_cache()
    per=pd.DataFrame(records)
    per.to_csv(out/"per_record_metrics.csv",index=False)
    all_data=pd.concat(sources,ignore_index=True)
    all_data.to_csv(out/"oof_per_residue.csv.gz",index=False,compression="gzip")
    (out/"prediction_manifest.jsonl").write_text("\n".join(json.dumps(r,ensure_ascii=False) for r in manifest)+"\n",encoding="utf-8")
    numeric=["mae","rmse","bias","spearman","pearson"]
    folds=per.groupby(["variant","fold","dataset","metric"])[numeric].mean().reset_index()
    overall=per.groupby(["variant","fold","metric"])[numeric].mean().reset_index().assign(dataset="all")
    folds=pd.concat([folds,overall],ignore_index=True)
    folds.to_csv(out/"source_fold_macro_metrics.csv",index=False)
    summary=folds.groupby(["variant","dataset","metric"])[numeric].agg(["mean","std"])
    summary.columns=[f"{name}_{scope}" for name,scope in summary.columns]
    summary.to_csv(out/"five_fold_mean_sample_sd.csv")
    report={"status":"complete","variants":["full","intrachain"],"oof_trained_fits":10,"additional_seed_fits":4,"epoch":40,"seed":7,
            "records_per_variant":len(per[per.metric=="rmsf"])/2,
            "effective_residues_per_variant":all_data.groupby("variant")[["valid_rmsf","valid_chi_entropy"]].sum().to_dict(),
            "protocol":"Corresponding heldout fold only, no five-model OOF ensemble; corrected input corpus; fixed epoch40.",
            "input_freeze_sha256":sha(ROOT/"training_freeze.json"),"export_code_sha256":sha(__file__)}
    (out/"evaluation_audit.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    print(summary.to_string(),flush=True)

if __name__=="__main__":main()
