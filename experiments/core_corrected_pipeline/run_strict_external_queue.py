"""Two additional serial GPU fits after the immutable primary queue finishes."""
from pathlib import Path
import os
import subprocess
import sys
import json
import time
ROOT=Path(__file__).resolve().parent
os.chdir(ROOT)
statefile=ROOT/"strict_queue_state.json"
while True:
    state=json.loads((ROOT/"queue_state.json").read_text())
    if state["status"]=="complete":break
    if "failed" in state["status"]:
        statefile.write_text(json.dumps({"status":"blocked_by_primary_failure","primary":state["status"]}))
        raise SystemExit("Primary job failed; strict training not started")
    statefile.write_text(json.dumps({"status":"waiting_for_primary_queue","checked_epoch":time.time()}))
    time.sleep(60)
env=dict(os.environ,PYTHONPATH=str(ROOT/"src"),OMP_NUM_THREADS="2",MKL_NUM_THREADS="2")
state={"status":"running","jobs":[]}
for variant in ["full","intrachain"]:
    item={"variant":variant,"started":time.time(),"status":"running"}
    state["jobs"].append(item)
    statefile.write_text(json.dumps(state,indent=2))
    ckpt=ROOT/f"runs/strict_{variant}/fold_0/epoch_040.pt"
    if not ckpt.exists():
        log=ROOT/f"logs/strict_{variant}.log"
        with log.open("a") as output:
            result=subprocess.run([sys.executable,"-u","-m","flexppi","train","--config",f"configs/strict_{variant}.yaml"],env=env,stdout=output,stderr=subprocess.STDOUT)
        if result.returncode or not ckpt.exists():
            state["status"]="training_failed"
            item["status"]="failed"
            statefile.write_text(json.dumps(state,indent=2))
            raise SystemExit("Strict training failed; preserve log")
    item.update(status="complete",finished=time.time())
    statefile.write_text(json.dumps(state,indent=2))
result=subprocess.run([sys.executable,"-u","evaluate_strict_external_models.py"],env=env)
state.update(status="complete" if result.returncode==0 else "evaluation_failed",evaluation_exit_code=result.returncode)
statefile.write_text(json.dumps(state,indent=2))
raise SystemExit(result.returncode)
