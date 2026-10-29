#!/usr/bin/env bash
set -euo pipefail
export PIP_CONFIG_FILE=/dev/null PIP_EXTRA_INDEX_URL= PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple PIP_DEFAULT_TIMEOUT=180
task_root=/root/flexppi_completion_20261010
if [ ! -x "$task_root/planet_env/bin/python" ]; then "$task_root/python311/bin/python" -m venv "$task_root/planet_env"; fi
# Share the already frozen CUDA/Torch install as a fallback path. The new
# environment's own scientific packages precede it, meeting PLANET's numpy<2.
printf '%s\n' "$task_root/env/lib/python3.11/site-packages" > "$task_root/planet_env/lib/python3.11/site-packages/frozen_torch_fallback.pth"
"$task_root/planet_env/bin/python" -m pip install 'numpy==1.26.4' 'scipy==1.15.3' 'pandas==2.3.3' 'scikit-learn==1.6.1'
"$task_root/planet_env/bin/python" -m pip install --no-deps "$task_root/planet_esm"
"$task_root/planet_env/bin/python" -m pip install -r "$task_root/baselines/planet_inference_requirements.txt" 'lightning==2.4.0' 'torchmetrics==1.4.3' 'h5py==3.14.0' 'numpy-indexed==0.3.7' 'wandb==0.21.4' 'seaborn==0.13.2' 'MDAnalysis==2.9.0' 'statsmodels==0.14.5' -c "$task_root/planet_constraints.txt"
tar -xf "$task_root/baselines/planet_payload.tar" -C "$task_root/baselines"
"$task_root/planet_env/bin/python" -m pip freeze > "$task_root/planet_environment_freeze.txt"
PYTHONPATH="$task_root/baselines/PLANET-MD" "$task_root/planet_env/bin/python" -c 'import torch,numpy; from planet_md import PlanetMD; from planet_md.esm3 import get_model; assert torch.__version__=="2.6.0+cu124"; assert numpy.__version__=="1.26.4"; print(torch.__version__,numpy.__version__)' > "$task_root/planet_environment_validation.txt"
date -u +%FT%TZ > "$task_root/PLANET_ENV_READY"
