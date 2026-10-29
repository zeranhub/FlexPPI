#!/usr/bin/env bash
# Evaluate the epoch-40 model of every fold on its own untouched test fold.
# The frozen protocol fixes epoch 40: no model selection from test results.
set -euo pipefail

REPO=/ranze/repos/FlexPPI
VENV="$REPO/.venv-linux/bin"
OUT="$REPO/artifacts/performance_first_v2/evaluation"
mkdir -p "$OUT"

for fold in 0 1 2 3 4; do
  CKPT="$REPO/artifacts/performance_first_v2/fold_${fold}/epoch_040.pt"
  TARGET="$OUT/fold_${fold}_test.json"
  if [ ! -f "$CKPT" ]; then
    echo "fold ${fold}: checkpoint $CKPT missing (training not finished)"
    continue
  fi
  if [ -f "$TARGET" ]; then
    echo "fold ${fold}: evaluation already exists at $TARGET"
    continue
  fi
  (cd "$REPO" && "$VENV/flexppi" evaluate \
    --config "configs/performance_first_v2_fold${fold}.yaml" \
    --checkpoint "$CKPT" \
    --split test \
    --bootstrap 200 \
    --output "$TARGET" 2>&1 | tee "$OUT/fold_${fold}_eval.log")
  echo "fold ${fold}: wrote $TARGET"
done
