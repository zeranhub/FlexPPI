#!/usr/bin/env bash
# Launch one information-ablation arm (sequence|geometry|intrachain) five-fold
# training in the background, one process per fold.  Idempotent and
# crash-restartable: a fold is skipped only when its process is running or its
# final epoch_040.pt exists.  The GPU is chosen by FLEXPPI_GPU.
# The container may be restarted by the scheduler; nohup+setsid keeps processes
# alive across SSH disconnects, and the driver relaunches killed folds.
set -euo pipefail

VARIANT="${1:?usage: launch_information_ablation.sh sequence|geometry|intrachain}"
case "$VARIANT" in
  sequence|geometry|intrachain) ;;
  *) echo "unknown variant: $VARIANT" >&2; exit 1 ;;
esac

REPO=/ranze/repos/FlexPPI
VENV="$REPO/.venv-linux/bin"
GPU="${FLEXPPI_GPU:-2}"
LOGDIR="$REPO/artifacts/information_ablation/logs"
mkdir -p "$LOGDIR"
cd "$REPO"

for fold in 0 1 2 3 4; do
  CONFIG="configs/information_ablation_${VARIANT}_fold${fold}.yaml"
  LOG="$LOGDIR/${VARIANT}_fold${fold}.log"
  FINAL="$REPO/artifacts/information_ablation/${VARIANT}/fold_${fold}/epoch_040.pt"
  if pgrep -f "[f]lexppi train --config $CONFIG" >/dev/null; then
    echo "$VARIANT fold ${fold}: already running (skipping)"
    continue
  fi
  if [ -f "$FINAL" ]; then
    echo "$VARIANT fold ${fold}: already finished ($FINAL exists, skipping)"
    continue
  fi
  nohup setsid env -u _CUDA_COMPAT_PATH \
    LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libcuda.so.1 \
    CUDA_VISIBLE_DEVICES=$GPU PYTHONUNBUFFERED=1 \
    PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    "$VENV/flexppi" train --config "$CONFIG" >> "$LOG" 2>&1 < /dev/null &
  echo "$VARIANT fold ${fold}: launched pid $! on GPU $GPU"
  sleep 2
done
