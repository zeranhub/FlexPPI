from pathlib import Path
import json
import shutil
import tarfile
import yaml
import sys
HERE = Path(__file__).resolve().parent
REV = HERE.parents[1]
sys.path.insert(0, str(HERE))
from prepare_corrected_labels import records, resolve

dest = HERE / "gpu_smoke_payload"
dest.mkdir(exist_ok=True)
record = max(records(), key=lambda r:r["metadata"]["residue_count"])
for key, name in [("structure_path", "smoke.pdb"), ("label_path", "smoke.npz")]:
    shutil.copy2(resolve(record[key]), dest / name)
    record[key] = "/root/flexppi_revision_20261002/gpu_smoke/" + name
record["trajectory_paths"] = []
(dest / "smoke.jsonl").write_text(json.dumps(record)+"\n",encoding="utf-8")
shutil.copytree(REV/"reproducibility/src",dest/"src",dirs_exist_ok=True,ignore=shutil.ignore_patterns("__pycache__"))
config = yaml.safe_load((REV/"reproducibility/configs/primary_fold0.yaml").read_text())
config["data"].update(state_manifest="/root/flexppi_revision_20261002/gpu_smoke/smoke.jsonl",test_state_manifest=None)
(dest / "smoke_config.yaml").write_text(yaml.safe_dump(config),encoding="utf-8")
shutil.copy2(HERE/"gpu_training_smoke.py",dest/"gpu_training_smoke.py")
with tarfile.open(HERE/"gpu_smoke_payload.tar.gz","w:gz") as tar:
    for path in dest.iterdir():
        tar.add(path,arcname=path.name)
print("Real longest frozen sample for discarded GPU memory test:",record["sample_id"],record["metadata"]["residue_count"])
