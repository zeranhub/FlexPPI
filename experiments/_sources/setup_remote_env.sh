#!/usr/bin/env bash
set -euo pipefail
export PIP_CONFIG_FILE=/dev/null
export PIP_EXTRA_INDEX_URL=
export PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple
export PIP_DEFAULT_TIMEOUT=180
export PIP_DISABLE_PIP_VERSION_CHECK=1
task_root=/root/flexppi_completion_20261010
mkdir -p "$task_root"
cd "$task_root"
if [ ! -x "$task_root/python311/bin/python" ]; then
  curl --fail --location --retry 3 --connect-timeout 20 --max-time 600 \
    https://repo.anaconda.com/miniconda/Miniconda3-py311_24.11.1-0-Linux-x86_64.sh \
    --output Miniconda3-py311_24.11.1-0-Linux-x86_64.sh
  sha256sum Miniconda3-py311_24.11.1-0-Linux-x86_64.sh > installer_sha256.txt
  bash Miniconda3-py311_24.11.1-0-Linux-x86_64.sh -b -p "$task_root/python311"
fi
if [ ! -x "$task_root/env/bin/python" ]; then
  "$task_root/python311/bin/python" -m venv "$task_root/env"
fi
"$task_root/env/bin/python" -m pip install --upgrade pip
"$task_root/env/bin/python" -m pip install torch==2.6.0 --index-url https://download.pytorch.org/whl/cu124
"$task_root/env/bin/python" -m pip install numpy==2.4.6 scipy==1.17.1 pandas==2.3.3 PyYAML==6.0.3 biopython joblib scikit-learn==1.6.1 mrcfile
"$task_root/env/bin/python" -m pip freeze > environment_freeze.txt
"$task_root/env/bin/python" -c 'import torch,numpy,scipy,pandas,yaml,sklearn; assert torch.cuda.is_available(); print(torch.__version__,torch.cuda.is_available(),torch.version.cuda,numpy.__version__,scipy.__version__,pandas.__version__,sklearn.__version__); print(torch.cuda.get_device_name(0))' > environment_validation.txt
date -u +%FT%TZ > ENV_READY
