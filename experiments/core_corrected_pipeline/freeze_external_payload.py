"""Portable external inputs; frozen cohort and every extraction outcome retained."""
from pathlib import Path
import json
import shutil
import tarfile
from prepare_corrected_labels import sha
HERE=Path(__file__).resolve().parent
REMOTE="/root/flexppi_revision_20261002/corrected_training_v1"

def main():
    source=HERE/"external_thermopcd"
    protocol=json.loads((source/"prespecified_protocol.json").read_text(encoding="utf-8"))
    extraction=json.loads((source/"extraction_summary.json").read_text(encoding="utf-8"))
    assert sha(source/"frozen_input_manifest.jsonl")==protocol["manifest_sha256"]
    assert extraction["prespecified"]==protocol["accepted_count"]==32
    rows=[json.loads(x) for x in (source/"ready_test.jsonl").read_text(encoding="utf-8").splitlines()]
    assert len(rows)==extraction["completed"]
    dest=HERE/"external_payload"
    external=dest/"external"
    if (external/"external_payload_freeze.json").exists():
        raise RuntimeError("External scoring payload already frozen; refusing mutation")
    external.mkdir(parents=True,exist_ok=True)
    for name in ["prespecified_protocol.json","frozen_input_manifest.jsonl","extraction_summary.json",
                 "protocol_addendum_20261002.json","protocol_addendum_gro_20261002.json"]:
        if (source/name).exists(): shutil.copy2(source/name,external/name)
    shutil.copytree(source/"audit",external/"audit",dirs_exist_ok=True)
    portable=[]
    for original in rows:
        r=json.loads(json.dumps(original))
        for key,directory,hash_name in [("structure_path","structures","input_sha256"),("label_path","labels","label_sha256")]:
            assert sha(r[key])==r["metadata"][hash_name]
            target=external/directory/Path(r[key]).name
            target.parent.mkdir(exist_ok=True)
            shutil.copy2(r[key],target)
            r[key]=REMOTE+"/external/"+directory+"/"+target.name
        portable.append(r)
    (external/"ready_test.jsonl").write_text("\n".join(json.dumps(r,ensure_ascii=False) for r in portable)+"\n",encoding="utf-8")
    report={"records":len(rows),"prespecified":32,"failures":extraction["failed"],
        "original_cohort_sha256":protocol["manifest_sha256"],
        "external_protocol_sha256":sha(external/"prespecified_protocol.json"),
        "data_sha256":{str(p.relative_to(external)):sha(p) for d in ["structures","labels"] for p in sorted((external/d).iterdir())},
        "purpose":"Portable source-out test inputs, frozen before corrected model scores; no change to cohort membership or initial model structures."}
    (external/"external_payload_freeze.json").write_text(json.dumps(report,indent=2),encoding="utf-8")
    shutil.copy2(HERE/"evaluate_external_models.py",dest/"evaluate_external_models.py")
    with tarfile.open(HERE/"external_payload.tar.gz","w:gz") as tar:
        for path in dest.iterdir(): tar.add(path,arcname=path.name)
    print(json.dumps({"scored_input_records":len(rows),"failures":len(extraction["failed"]),"payload":str(HERE/"external_payload.tar.gz")}))

if __name__=="__main__":main()
