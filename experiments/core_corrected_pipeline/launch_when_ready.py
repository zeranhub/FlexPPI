"""Wait for QC-complete labels, freeze hashes, upload, then launch serial GPU fits."""
from pathlib import Path
import sys
import json
import time
import subprocess
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
from prepare_corrected_labels import records
from freeze_training_payload import freeze, REMOTE
statefile=HERE/"launch_state.json"
state={"status":"waiting_for_corrected_labels","created_epoch":time.time(),"training_started":False}
expected=records()
while True:
    complete=0
    excluded=[]
    unresolved=[]
    missing=0
    for r in expected:
        p=HERE/"dataset/audit"/(r["sample_id"]+".json")
        if not p.exists():missing+=1;continue
        try:a=json.loads(p.read_text(encoding="utf-8"))
        except (json.JSONDecodeError,PermissionError):missing+=1;continue
        if a["status"]=="complete":complete+=1
        elif "exceeds predeclared 1024" in a.get("error",""):excluded.append(r["sample_id"])
        else:unresolved.append({"sample_id":r["sample_id"],"error":a.get("error")})
    state.update(checked_epoch=time.time(),complete_records=complete,pending_records=missing,
                 predeclared_length_exclusions=excluded,unresolved_failures=unresolved)
    statefile.write_text(json.dumps(state,indent=2),encoding="utf-8")
    if missing==0:
        if unresolved:
            state["status"]="qc_requires_repair"
            statefile.write_text(json.dumps(state,indent=2),encoding="utf-8")
            raise SystemExit("QC failed; training not started. See launch_state.json")
        break
    time.sleep(60)
print("All corrected input audits resolved; freezing portable training corpus.",flush=True)
freeze()
state["status"]="uploading_frozen_training_payload"
statefile.write_text(json.dumps(state,indent=2),encoding="utf-8")
subprocess.run(["scp","-q",str(HERE/"training_payload.tar.gz"),"hit-gpu:/root/flexppi_revision_20261002/"],check=True)
command="mkdir -p " + REMOTE + " && tar -xzf /root/flexppi_revision_20261002/training_payload.tar.gz -C " + REMOTE
subprocess.run(["ssh","-o","BatchMode=yes","hit-gpu",command],check=True)
launch_code=("import subprocess,pathlib; p=pathlib.Path('"+REMOTE+"'); "
    "f=(p/'queue.log').open('a'); "
    "s=subprocess.Popen(['/root/flexppi_revision_20261002/training_env/bin/python','-u',str(p/'run_training_queue.py')],"
    "cwd=p,stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True); "
    "(p/'queue.pid').write_text(str(s.pid)); print(s.pid)")
# Pass one shell-quoted Python program, with no interpolation of user content.
shell_code="/root/flexppi_revision_20261002/training_env/bin/python -c '"+launch_code.replace("'","'\"'\"'")+"'"
pid=subprocess.check_output(["ssh","-o","BatchMode=yes","hit-gpu",shell_code],text=True).strip()
state.update(status="remote_training_queue_running",training_started=True,remote_pid=pid,
             remote_directory=REMOTE,launched_epoch=time.time(),fits=14,epochs_per_fit=40)
statefile.write_text(json.dumps(state,indent=2),encoding="utf-8")
print(json.dumps(state,indent=2),flush=True)
