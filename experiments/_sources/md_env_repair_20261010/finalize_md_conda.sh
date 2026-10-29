#!/usr/bin/env bash
set -euo pipefail
TASK=/root/flexppi_completion_20261010
PY="$TASK/md_conda/bin/python"
test -x "$PY"
export PIP_CONFIG_FILE=/dev/null
export PIP_CACHE_DIR="$TASK/md_pip_cache"
export PIP_DISABLE_PIP_VERSION_CHECK=1
"$PY" -m pip install --no-index --find-links "$TASK/md_install_sources/wheels" numpy==2.2.6 scipy==1.15.3 MDAnalysis==2.9.0
"$PY" -m pip install --no-deps --no-build-isolation "$TASK/md_install_sources/pdbfixer_94cfa4c0ca551cdc5f13320f9a658efd59f2b881.tar.gz"
"$PY" -m pip check
"$PY" "$TASK/md_env_repair_20261010/verify_md_environment.py"
"$PY" -m pip freeze > "$TASK/md_environment_pip_freeze.txt"
CONDA_PKGS_DIRS="$TASK/md_conda_pkgs" CONDA_ENVS_PATH="$TASK/md_conda_envs" "$TASK/python311/bin/conda" list -p "$TASK/md_conda" --explicit > "$TASK/md_environment_conda_explicit.txt"
"$PY" "$TASK/kix_md_triplicate_rebuild_20261010/code/analyze_job.py" --self-test
"$PY" "$TASK/kix_md_triplicate_rebuild_20261010/code/quality_one_repeat.py" --self-test
"$PY" "$TASK/kix_md_triplicate_rebuild_20261010/code/run_md.py" --help
"$PY" "$TASK/kix_md_triplicate_rebuild_20261010/code/queue_triplicate.py" --help
"$PY" "$TASK/kix_md_triplicate_rebuild_20261010/code/audit_rebuilt_prepared.py" --help
"$PY" "$TASK/kix_md_triplicate_rebuild_20261010/code/analyze_triplicate_delta.py" --help
