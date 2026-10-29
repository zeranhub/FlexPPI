"""Separate source/family-purged experiment; preserve the primary 14-fit freeze."""
from pathlib import Path
import json
import shutil
import tarfile
import yaml
from prepare_corrected_labels import sha
HERE=Path(__file__).resolve().parent
REMOTE="/root/flexppi_revision_20261002/corrected_training_v1"

def main():
    audit=json.loads((HERE/"sequence_audit/sequence_audit_summary.json").read_text(encoding="utf-8"))
    if not audit["strict_fit_feasible"]:raise RuntimeError("Prespecified strict fit feasibility rule failed")
    target=HERE/"strict_external_payload"
    if (target/"strict_external_freeze.json").exists():raise RuntimeError("Strict experiment already frozen")
    target.mkdir(exist_ok=True)
    (target/"splits").mkdir(exist_ok=True)
    (target/"configs").mkdir(exist_ok=True)
    rows=[json.loads(x) for x in (HERE/"sequence_audit/purged_training_local.jsonl").read_text(encoding="utf-8").splitlines()]
    assert len(rows)==2431
    for r in rows:
        r["structure_path"]=REMOTE+"/dataset/structures/"+Path(r["structure_path"]).name
        r["label_path"]=REMOTE+"/dataset/labels/"+Path(r["label_path"]).name
    (target/"splits/strict_external_train.jsonl").write_text("\n".join(json.dumps(r) for r in rows)+"\n",encoding="utf-8")
    for variant in ["full","intrachain"]:
        c=yaml.safe_load((HERE/f"training_payload/configs/{variant}_fold0.yaml").read_text(encoding="utf-8"))
        c["output_dir"]=f"runs/strict_{variant}/fold_0"
        c["data"].update(state_manifest="splits/strict_external_train.jsonl",test_state_manifest=None)
        (target/f"configs/strict_{variant}.yaml").write_text(yaml.safe_dump(c,sort_keys=False),encoding="utf-8")
    evaluator=(HERE/"evaluate_external_models.py").read_text(encoding="utf-8")
    evaluator=evaluator.replace('out=ROOT/"external_evaluation"','out=ROOT/"strict_external_evaluation"')
    evaluator=evaluator.replace('[(k,7) for k in range(5)]+[(0,11),(0,23)]','[(0,7)]')
    evaluator=evaluator.replace('checkpoint_variant=variant if seed==7 else f"{variant}_seed{seed}"','checkpoint_variant=f"strict_{variant}"')
    evaluator=evaluator.replace('assert members.shape[0]==5','assert members.shape[0]==1')
    evaluator=evaluator.replace('role="five_fold_ensemble"','role="single_strict_model"')
    evaluator=evaluator.replace('"checkpoints":14','"checkpoints":2')
    evaluator=evaluator.replace('"external_prediction":"Equal average of five fold models per arm; each external input is absent from all training folds."',
        '"external_prediction":"One fixed epoch40 model per arm trained only on the source/family-purged 2431 records; external components excluded."')
    evaluator=evaluator.replace('"seed_stability":"Fixed fold0 membership, seeds7/11/23; reported separately from ensemble predictions and fold SD."',
        '"seed_stability":"This separate strict experiment uses seed7 only; primary fixed-fold stability experiment remains separate."')
    evaluator=evaluator.replace('"training_freeze_sha256":sha(ROOT/"training_freeze.json")','"training_freeze_sha256":sha(ROOT/"strict_external_freeze.json")')
    evaluator=evaluator.replace('"homology":"Must join chain-level training overlap audit before claims of unseen-family generalization."',
        '"homology":"All external-linked sequence/biological components excluded from strict training; one external connected component, no family-independent CI."')
    (target/"evaluate_strict_external_models.py").write_text(evaluator,encoding="utf-8")
    shutil.copy2(HERE/"run_strict_external_queue.py",target/"run_strict_external_queue.py")
    shutil.copy2(HERE/"strict_external_prespecified_plan.json",target/"strict_external_prespecified_plan.json")
    report={"training_records":2431,"training_source_counts":audit["retained_by_source"],
        "training_multichain_records":408,"excluded_training_records":632,"external_records":32,"external_components":1,
        "seed":7,"epoch":40,"fits":2,"no_test_tuning":True,"primary_freeze_sha256":sha(HERE/"training_payload/training_freeze.json"),
        "prespecified_plan_sha256":sha(HERE/"strict_external_prespecified_plan.json"),
        "sequence_search_sha256":audit["sequence_search_sha256"],
        "execution_hashes":{str(p.relative_to(target)):sha(p) for p in target.rglob("*") if p.is_file()},
        "statistical_limitation":"32 source-out complexes are linked by shared protein families; they are one dependence component, not 32 independent families."}
    (target/"strict_external_freeze.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    with tarfile.open(HERE/"strict_external_payload.tar.gz","w:gz") as tar:
        for p in target.iterdir():tar.add(p,arcname=p.name)
    print(json.dumps(report,indent=2))

if __name__=="__main__":main()
