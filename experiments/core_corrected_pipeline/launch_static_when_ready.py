"""Launch fixed-budget CPU controls once both corrected input freezes exist."""
from pathlib import Path
import json
import subprocess
import sys
import time
HERE=Path(__file__).resolve().parent
statefile=HERE/"static_launch_state.json"
while not ((HERE/"training_payload/training_freeze.json").exists() and (HERE/"external_payload/external/external_payload_freeze.json").exists()):
    statefile.write_text(json.dumps({"status":"waiting_for_input_freezes","checked_epoch":time.time()}),encoding="utf-8")
    time.sleep(60)
statefile.write_text(json.dumps({"status":"running","started_epoch":time.time()}),encoding="utf-8")
python=sys.executable
with (HERE/"static_corrected_training.log").open("a",encoding="utf-8") as log:
    result=subprocess.run([python,"-u",str(HERE/"run_corrected_static.py"),"--external"],stdout=log,stderr=subprocess.STDOUT)
statefile.write_text(json.dumps({"status":"complete" if result.returncode==0 else "failed","exit_code":result.returncode,"finished_epoch":time.time()}),encoding="utf-8")
raise SystemExit(result.returncode)
