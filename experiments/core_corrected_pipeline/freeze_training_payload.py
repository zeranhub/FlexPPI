# -*- coding: utf-8 -*-
"""Release corrected inputs only after complete QC; portable server paths."""
from pathlib import Path
import json
import shutil
import tarfile
import yaml
import sys
import time
import hashlib
HERE = Path(__file__).resolve().parent
REV = HERE
sys.path.insert(0,str(HERE))
from prepare_corrected_labels import records, sha

REMOTE = "/root/flexppi_revision_20261002/corrected_training_v1"

def portable(r):
    r = json.loads(json.dumps(r))
    for key, directory in [("structure_path","structures"),("label_path","labels")]:
        r[key] = REMOTE + "/dataset/" + directory + "/" + Path(r[key]).name
    r["trajectory_paths"] = []
    return r

def freeze():
    expected = records()
    audits, missing, failed = [],[],[]
    for r in expected:
        file = HERE / "dataset/audit" / (r["sample_id"]+".json")
        if not file.exists():
            missing.append(r["sample_id"])
            continue
        a = json.loads(file.read_text(encoding="utf-8"))
        if a["status"] == "complete":
            for key, value in a["corrected_hashes"].items():
                if sha(a["record"]["structure_path" if key == "structure" else "label_path"]) != value:
                    raise RuntimeError("Corrected input hash changed: "+r["sample_id"])
            audits.append(a)
        elif "exceeds predeclared 1024" in a.get("error", ""):
            failed.append(a)
        else:
            raise RuntimeError("Unresolved extraction failure: " + r["sample_id"] + " " + a.get("error", ""))
    if missing:
        raise RuntimeError(f"Extraction incomplete: {len(missing)} records")
    dest = HERE / "training_payload"
    if (dest/"training_freeze.json").exists():
        raise RuntimeError("Training dataset already frozen; refusing mutation")
    dest.mkdir(exist_ok=True)
    output_records = []
    for a in audits:
        r = a["record"]
        for key,directory in [("structure_path","structures"),("label_path","labels")]:
            path = dest/"dataset"/directory/Path(r[key]).name
            path.parent.mkdir(parents=True,exist_ok=True)
            shutil.copy2(r[key],path)
        output_records.append(portable(r))
    split = dest/"splits"
    split.mkdir(exist_ok=True)
    def write_manifest(name,rs):
        (split/name).write_text("\n".join(json.dumps(r,ensure_ascii=False,sort_keys=True) for r in rs)+"\n",encoding="utf-8")
    write_manifest("all_ready.jsonl",output_records)
    configs = dest/"configs"
    configs.mkdir(exist_ok=True)
    for fold in range(5):
        for role in ["train","test"]:
            selected = [r for r in output_records if (r["metadata"]["fold"]==fold)==(role=="test")]
            write_manifest(f"fold_{fold}_{role}.jsonl",selected)
        for variant in ["full","intrachain"]:
            c = yaml.safe_load((REV/f"reproducibility/configs/primary_fold{fold}.yaml").read_text(encoding="utf-8"))
            c["device"] = "cuda"
            c["output_dir"] = f"runs/{variant}/fold_{fold}"
            c["data"].update(state_manifest=f"splits/fold_{fold}_train.jsonl",test_state_manifest=f"splits/fold_{fold}_test.jsonl",
                               path_mappings={},num_workers=0)
            c["model"]["information_ablation"] = variant
            c["train"].update(resume=True,log_every=100)
            (configs/f"{variant}_fold{fold}.yaml").write_text(yaml.safe_dump(c,sort_keys=False),encoding="utf-8")
    # Seed stability uses the same fixed fold-0 training membership. It is
    # separate from five-fold SD and does not select any epoch from test scores.
    for variant in ["full","intrachain"]:
        for seed in [11,23]:
            c=yaml.safe_load((configs/f"{variant}_fold0.yaml").read_text(encoding="utf-8"))
            c["seed"]=seed
            c["output_dir"]=f"runs/{variant}_seed{seed}/fold_0"
            (configs/f"{variant}_seed{seed}_fold0.yaml").write_text(yaml.safe_dump(c,sort_keys=False),encoding="utf-8")
    shutil.copytree(REV/"reproducibility/src",dest/"src",dirs_exist_ok=True,ignore=shutil.ignore_patterns("__pycache__"))
    for name in ["run_training_queue.py","evaluate_corrected_models.py"]:
        shutil.copy2(HERE/name,dest/name)
    report = {"timestamp_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),
              "records":len(output_records),"excluded_length_qc":[{"sample_id":a["sample_id"],"reason":a["error"]} for a in failed],
              "source_records":dict(__import__('collections').Counter(r["dataset"] for r in output_records)),
              "valid_counts":[sum(a["valid_counts"][m] for a in audits) for m in [0,1]],
              "actual_fold_membership":"Inherited from frozen original biological combination split; length exclusions before scoring. Not all-family isolation.",
              "training":"from scratch, seed7, epoch40, original source-system sampling, no test-driven checkpoint selection",
              "expected_neural_fits":14,"arms":["full","intrachain"],"seed_stability":"fixed fold0 seeds7,11,23, separate from fold variation",
              "external_source":"ThermoPCD separately prespecified, never in training",
              "input_hashes":{str(p.relative_to(dest)):sha(p) for p in sorted((dest/"dataset").rglob("*")) if p.is_file()},
              "execution_hashes":{str(p.relative_to(dest)):sha(p) for p in list((dest/"src").rglob("*.py"))+list(configs.glob("*.yaml"))+list(split.glob("*.jsonl"))},
              "old_originals_modified":False,"inputs_complete_before_training":True}
    (dest/"training_freeze.json").write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding="utf-8")
    with tarfile.open(HERE/"training_payload.tar.gz","w:gz") as tar:
        for path in dest.iterdir():tar.add(path,arcname=path.name)
    print("Corrected training dataset frozen:",len(output_records),report["valid_counts"],flush=True)
    return report

if __name__=="__main__": freeze()
