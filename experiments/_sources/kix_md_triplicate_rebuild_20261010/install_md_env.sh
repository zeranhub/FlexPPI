#!/usr/bin/env bash
# Only installs into the new experiment's dedicated MD environment.
# Supply the absolute Python 3.11 executable as argument 1.
set -euo pipefail
export PIP_CONFIG_FILE=/dev/null
export PIP_EXTRA_INDEX_URL=
export PIP_INDEX_URL=https://pypi.tuna.tsinghua.edu.cn/simple
export PIP_DEFAULT_TIMEOUT=180
export PIP_DISABLE_PIP_VERSION_CHECK=1
TASK_BASE=/root/flexppi_completion_20261010
MD_ENV="$TASK_BASE/md_env"
PYTHON311="${1:?Supply an absolute Python 3.11 executable}"
export PIP_CACHE_DIR="$TASK_BASE/md_pip_cache"
export PYTHONNOUSERSITE=1
mkdir -p "$TASK_BASE/md_install_sources"
if [[ "$($PYTHON311 -c 'import sys; print("%d.%d" % sys.version_info[:2])')" != '3.11' ]]; then
  echo 'The supplied interpreter must be Python 3.11' >&2; exit 2
fi
if [[ ! -x "$MD_ENV/bin/python" ]]; then
  "$PYTHON311" -m venv "$MD_ENV"
fi
MD_PYTHON="$MD_ENV/bin/python"
"$MD_PYTHON" -m pip install --upgrade pip setuptools wheel
"$MD_PYTHON" -m pip install 'numpy==2.2.6' 'scipy==1.15.3' 'OpenMM==8.4.0.post2' 'OpenMM-CUDA-12==8.4.0.post2' 'MDAnalysis==2.9.0'
# Official pdbfixer v1.12 commit; setup.py explicitly declares version 1.12.0.
PDBFIXER_COMMIT=94cfa4c0ca551cdc5f13320f9a658efd59f2b881
PDBFIXER_ARCHIVE="$TASK_BASE/md_install_sources/pdbfixer_${PDBFIXER_COMMIT}.tar.gz"
curl --fail --location --retry 3 "https://github.com/openmm/pdbfixer/archive/${PDBFIXER_COMMIT}.tar.gz" --output "$PDBFIXER_ARCHIVE"
sha256sum "$PDBFIXER_ARCHIVE" > "$PDBFIXER_ARCHIVE.sha256"
"$MD_PYTHON" -m pip install --no-deps "$PDBFIXER_ARCHIVE"
"$MD_PYTHON" -m pip check
"$MD_PYTHON" -m pip freeze > "$TASK_BASE/md_install_sources/environment_freeze.txt"
"$MD_PYTHON" -c 'import importlib.metadata as m,openmm; print({p:m.version(p) for p in ["numpy","scipy","OpenMM","OpenMM-CUDA-12","MDAnalysis","pdbfixer"]}); print("Platforms:",[openmm.Platform.getPlatform(i).getName() for i in range(openmm.Platform.getNumPlatforms())])'
date -u +%FT%TZ > "$TASK_BASE/MD_ENV_READY"
