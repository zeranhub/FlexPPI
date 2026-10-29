"""Serial corrected-data fits, resumable, fixed final epoch, with explicit state."""
from pathlib import Path
import subprocess
import sys
import json
import time
import os
ROOT = Path(__file__).resolve().parent
os.chdir(ROOT)
env = dict(os.environ, PYTHONPATH=str(ROOT/"src"), OMP_NUM_THREADS="2", MKL_NUM_THREADS="2")
jobs = [(v,k) for v in ["full","intrachain"] for k in range(5)]
jobs += [(f"{v}_seed{s}",0) for v in ["full","intrachain"] for s in [11,23]]
statefile = ROOT/"queue_state.json"
state = {"started_utc":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),"jobs":[],"status":"running"}
for variant,fold in jobs:
    ckpt = ROOT/f"runs/{variant}/fold_{fold}/epoch_040.pt"
    item = {"variant":variant,"fold":fold,"status":"running","started":time.time()}
    state["jobs"].append(item)
    statefile.write_text(json.dumps(state,indent=2))
    if ckpt.exists():
        item["status"] = "already_complete"
        continue
    log = ROOT/f"logs/{variant}_fold{fold}.log"
    log.parent.mkdir(exist_ok=True)
    print("Starting",variant,fold,flush=True)
    with log.open("a") as output:
        result = subprocess.run([sys.executable,"-u","-m","flexppi","train","--config",f"configs/{variant}_fold{fold}.yaml"],env=env,stdout=output,stderr=subprocess.STDOUT)
    item.update(exit_code=result.returncode,finished=time.time(),status="complete" if result.returncode==0 and ckpt.exists() else "failed")
    statefile.write_text(json.dumps(state,indent=2))
    if item["status"] == "failed":
        state["status"] = "failed"
        statefile.write_text(json.dumps(state,indent=2))
        raise SystemExit("Training failed; retained log and checkpoints: "+str(log))
state["status"] = "training_complete"
statefile.write_text(json.dumps(state,indent=2))
print("All 14 corrected-data fits complete. Starting fixed OOF export.",flush=True)
result = subprocess.run([sys.executable,"-u","evaluate_corrected_models.py"],env=env)
state.update(status="oof_complete" if result.returncode==0 else "evaluation_failed",evaluation_exit_code=result.returncode)
statefile.write_text(json.dumps(state,indent=2))
if result.returncode==0:
    if not (ROOT/"external/ready_test.jsonl").exists():
        state["status"]="awaiting_external_payload"
        statefile.write_text(json.dumps(state,indent=2))
        raise SystemExit("External payload absent; OOF retained, external evaluation pending")
    state["status"]="external_evaluation_running"
    statefile.write_text(json.dumps(state,indent=2))
    result=subprocess.run([sys.executable,"-u","evaluate_external_models.py"],env=env)
    state.update(status="complete" if result.returncode==0 else "external_evaluation_failed",external_exit_code=result.returncode)
    statefile.write_text(json.dumps(state,indent=2))
raise SystemExit(result.returncode)
