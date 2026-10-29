from pathlib import Path
from importlib.metadata import version
import datetime,json,platform,subprocess
import torch
value={'python':platform.python_version(),'captured_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),
       'packages':{name:version(name) for name in ['torch','numpy','scipy','pandas','PyYAML','scikit-learn']},
       'cuda_runtime':torch.version.cuda,'gpu_read_only_query':subprocess.check_output(
           ['nvidia-smi','--query-gpu=name','--format=csv,noheader'],text=True).strip(),
       'cuda_context_initialized':torch.cuda.is_initialized()}
assert value['cuda_context_initialized'] is False
Path(__file__).with_name('current_runtime_environment.json').write_text(json.dumps(value,indent=2),encoding='utf8')
print(json.dumps(value))
