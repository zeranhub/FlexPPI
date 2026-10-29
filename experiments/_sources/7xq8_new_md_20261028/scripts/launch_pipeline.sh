#!/usr/bin/env bash
set -euo pipefail
cd /ranze/runs/flexppi_7xq8_md_20261029
exec 9>pipeline.lock
flock -n 9 || { echo 'Another pipeline owns pipeline.lock'; exit 2; }
export LD_LIBRARY_PATH=/usr/local/cuda/compat/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}
export OPENMM_CPU_THREADS=8
export OMP_NUM_THREADS=8
PY=/ranze/envs/flexppi-kix-md-20260914/bin/python
printf '%s\n' "$(date -u) Pipeline start; 7XQ8 one 100 ns production trajectory"
if [[ ! -s protein_prepared.pdb ]]; then
  "$PY" -u prepare_membrane.py --run-dir "$PWD" --stage protein --gpu 1
fi
if [[ ! -s membrane_preparation_audit.json ]]; then
  "$PY" -u prepare_membrane.py --run-dir "$PWD" --stage membrane --gpu 1
fi
"$PY" -u run_membrane_md.py --run-dir "$PWD" --gpu 1 --production-ns 100
printf '%s\n' "$(date -u) Pipeline finished"
