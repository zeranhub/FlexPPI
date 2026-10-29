"""CPU smoke wrapper; sets two Torch threads without modifying the adapter."""
from pathlib import Path
import runpy
import torch

torch.set_num_threads(2)
print('CPU SMOKE: torch_threads=2; released float32 weights; no dtype conversion.', flush=True)
runpy.run_path(str(Path('/root/flexppi_completion_20261010/baselines/run_planet_md_3063.py')), run_name='__main__')
