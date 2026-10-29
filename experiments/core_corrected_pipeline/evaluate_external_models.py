"""Prespecified source-out predictions from corrected fixed-epoch checkpoints."""
from pathlib import Path
import sys
import json
import time
import numpy as np
import pandas as pd
import torch

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT/"src"))
from flexppi.config import ModelConfig
from flexppi.data.dataset import StateDataset,collate_states
from flexppi.data.structures import parse_pdb
from flexppi.model.information_ablation import build_information_ablation
from flexppi.trainer import move_to_device
from evaluate_corrected_models import sha,statistics

def main():
    torch.set_num_threads(2)
    external=ROOT/"external"
    protocol=json.loads((external/"prespecified_protocol.json").read_text(encoding="utf-8"))
    extraction=json.loads((external/"extraction_summary.json").read_text(encoding="utf-8"))
    expected=[json.loads(x) for x in (external/"ready_test.jsonl").read_text(encoding="utf-8").splitlines()]
    assert len(expected)==extraction["completed"]
    assert extraction["prespecified"]==protocol["accepted_count"]
    dataset=StateDataset(str(external/"ready_test.jsonl"),max_residues=1024)
    assert len(dataset)==len(expected)>0
    out=ROOT/"external_evaluation"
    out.mkdir(exist_ok=True)
    manifests=[]
    per_record=[]
    per_residue=[]
    contract={"training_freeze_sha256":sha(ROOT/"training_freeze.json"),
        "prespecified_protocol_sha256":sha(external/"prespecified_protocol.json"),
        "external_payload_sha256":sha(external/"external_payload_freeze.json"),
        "export_code_sha256":sha(__file__)}
    for variant in ["full","intrachain"]:
        means=[[] for _ in range(len(dataset))]
        variances=[[] for _ in range(len(dataset))]
        valid_masks=[]
        for fold,seed in [(k,7) for k in range(5)]+[(0,11),(0,23)]:
            checkpoint_variant=variant if seed==7 else f"{variant}_seed{seed}"
            path=ROOT/f"runs/{checkpoint_variant}/fold_{fold}/epoch_040.pt"
            checkpoint=torch.load(path,map_location="cpu",weights_only=False)
            c=checkpoint["config"]
            assert checkpoint["epoch"]==40 and c["seed"]==seed
            assert c["model"]["information_ablation"]==variant
            model=build_information_ablation(ModelConfig(**c["model"]),variant).to("cuda").eval()
            model.load_state_dict(checkpoint["model"],strict=True)
            checkpoint_hash=sha(path)
            sub=out/checkpoint_variant/f"fold_{fold}"
            sub.mkdir(parents=True,exist_ok=True)
            for i,r in enumerate(dataset.records):
                ex=dataset[i]
                structure=parse_pdb(r.structure_path)
                batch=move_to_device(collate_states([ex]),torch.device("cuda"))
                torch.cuda.synchronize()
                started=time.perf_counter()
                with torch.inference_mode(): p=model.forward_state(batch)
                torch.cuda.synchronize()
                seconds=time.perf_counter()-started
                mean=p["mean"][0].cpu().numpy()
                variance=p["log_variance"][0].exp().cpu().numpy()
                valid=(ex["label_mask"] & ex["target_mask"][:,None] & p["mask"][0].cpu()[:,None]).numpy()
                y=ex["label_values"].numpy()
                assert np.isfinite(mean).all() and np.isfinite(variance).all()
                if fold==0 and seed==7: valid_masks.append(valid)
                else: assert np.array_equal(valid,valid_masks[i])
                np.savez_compressed(sub/(r.sample_id+".npz"),mean=mean,model_variance=variance,
                    target=y,valid=valid,label_variance=ex["label_variances"].numpy(),
                    amino_acid=np.asarray(list(structure.sequence)),chain=np.asarray(structure.chain_ids),
                    author_residue_id=np.asarray(structure.author_residue_ids))
                manifests.append(dict(sample_id=r.sample_id,variant=variant,fold=fold,seed=seed,
                    checkpoint_sha256=checkpoint_hash,input_sha256=sha(r.structure_path),label_sha256=sha(r.label_path),
                    forward_seconds=seconds,residue_count=len(mean),cache_path=str(sub/(r.sample_id+".npz"))))
                for m,metric in enumerate(c["model"]["metrics"]):
                    if valid[:,m].any(): per_record.append(dict(sample_id=r.sample_id,variant=variant,fold=fold,seed=seed,
                        role="individual_model",metric=metric,**statistics(y[valid[:,m],m],mean[valid[:,m],m])))
                if seed==7:
                    means[i].append(mean)
                    variances[i].append(variance)
            del model
            torch.cuda.empty_cache()
            print("External forward",variant,fold,seed,len(dataset),flush=True)
        ensemble=out/variant/"ensemble"
        ensemble.mkdir(exist_ok=True)
        for i,r in enumerate(dataset.records):
            ex=dataset[i]
            structure=parse_pdb(r.structure_path)
            members=np.asarray(means[i],dtype=np.float64)
            member_var=np.asarray(variances[i],dtype=np.float64)
            assert members.shape[0]==5
            mean=members.mean(0)
            within=member_var.mean(0)
            between=members.var(0,ddof=0)
            variance=within+between
            y=ex["label_values"].numpy()
            valid=valid_masks[i]
            np.savez_compressed(ensemble/(r.sample_id+".npz"),mean=mean,model_variance=variance,
                within_model_variance=within,between_fold_variance=between,member_means=members,
                target=y,valid=valid,label_variance=ex["label_variances"].numpy(),
                amino_acid=np.asarray(list(structure.sequence)),chain=np.asarray(structure.chain_ids),
                author_residue_id=np.asarray(structure.author_residue_ids))
            row=dict(sample_id=r.sample_id,variant=variant,residue_index=np.arange(len(mean)),
                     amino_acid=list(structure.sequence),chain=structure.chain_ids,author_residue_id=structure.author_residue_ids)
            for m,metric in enumerate(["rmsf","chi_entropy"]):
                if valid[:,m].any(): per_record.append(dict(sample_id=r.sample_id,variant=variant,fold=-1,seed=7,
                    role="five_fold_ensemble",metric=metric,**statistics(y[valid[:,m],m],mean[valid[:,m],m])))
                row.update({f"target_{metric}":y[:,m],f"pred_{metric}":mean[:,m],f"valid_{metric}":valid[:,m],
                    f"sd_{metric}":np.sqrt(variance[:,m]),f"between_fold_sd_{metric}":np.sqrt(between[:,m])})
            per_residue.append(pd.DataFrame(row))
    df=pd.DataFrame(per_record)
    df.to_csv(out/"per_record_metrics.csv",index=False)
    df.groupby(["variant","role","fold","seed","metric"])[["mae","rmse","bias","spearman","pearson"]].mean().to_csv(out/"macro_metrics.csv")
    pd.concat(per_residue,ignore_index=True).to_csv(out/"per_residue.csv.gz",compression="gzip",index=False)
    (out/"prediction_manifest.jsonl").write_text("\n".join(json.dumps(x) for x in manifests)+"\n",encoding="utf-8")
    audit={"status":"complete",**contract,"prespecified_records":protocol["accepted_count"],
        "scored_records":len(dataset),"extraction_failures":extraction["failed"],"checkpoints":14,
        "epoch":40,"external_prediction":"Equal average of five fold models per arm; each external input is absent from all training folds.",
        "seed_stability":"Fixed fold0 membership, seeds7/11/23; reported separately from ensemble predictions and fold SD.",
        "uncertainty":"Within-model variance plus population between-fold mean variance; uncalibrated, not a biological confidence interval.",
        "not_blind":"Earlier structure/pair exploration is disclosed in the frozen external protocol.",
        "homology":"Must join chain-level training overlap audit before claims of unseen-family generalization."}
    (out/"evaluation_audit.json").write_text(json.dumps(audit,indent=2),encoding="utf-8")

if __name__=="__main__":main()
