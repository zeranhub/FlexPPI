"""Read-only capability inventory; creates no remote jobs or installations."""
import importlib.util
import json
import os
from pathlib import Path
import platform
import shutil
import sys

ROOT = Path(__file__).resolve().parent
handles = []
dll = Path(sys.prefix) / 'Library' / 'bin'
if os.name == 'nt' and dll.is_dir():
    handles.append(os.add_dll_directory(str(dll)))
import numpy as np
import torch

modules = ['transformers', 'esm', 'openmm', 'simtk', 'pdbfixer', 'pyrosetta', 'modeller', 'Bio', 'pymol']
environments = [Path('D:/Anaconda/envs') / name for name in ['MixTCRpred', 'pvb_cpu', 'PVB_env', 'pytorch']]
environments.append(Path('E:/Project/FlexPPI/.venv-md'))
report = {'python_executable': sys.executable, 'python': sys.version, 'platform': platform.platform(),
    'cpu': platform.processor(), 'logical_cores': os.cpu_count(), 'torch': torch.__version__,
    'torch_cuda_available': torch.cuda.is_available(), 'default_torch_threads': torch.get_num_threads(),
    'numpy': np.__version__, 'module_specs_in_current_environment': {name: bool(importlib.util.find_spec(name)) for name in modules},
    'commands_on_PATH': {name: shutil.which(name) for name in ['foldx', 'foldx.exe', 'rosetta_scripts', 'modeller', 'obabel', 'nvidia-smi']},
    'other_environment_package_directories': {},
    'known_esm_cache': {'path': 'E:/Project/protein_flexibility_models/torch_cache/hub/checkpoints/esm2_t33_650M_UR50D.pt',
        'exists': Path('E:/Project/protein_flexibility_models/torch_cache/hub/checkpoints/esm2_t33_650M_UR50D.pt').is_file(),
        'needed_by_target_model': False},
    'target_sequence_backend': 'learned', 'ESM_embeddings_are_not_inputs': True,
    'installs_performed': [], 'remote_computation_performed': False}
for env in environments:
    site = env / 'Lib' / 'site-packages'
    names = [p.name for p in site.iterdir()] if site.exists() else []
    report['other_environment_package_directories'][str(env)] = {
        module: [name for name in names if name.lower() == module.lower() or name.lower().startswith(module.lower()+'-')]
        for module in ['torch'] + modules}
(ROOT / 'environment_inventory.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
print(json.dumps({'status': 'complete', 'report': str(ROOT / 'environment_inventory.json'),
    'torch': torch.__version__, 'cuda': torch.cuda.is_available()}))
