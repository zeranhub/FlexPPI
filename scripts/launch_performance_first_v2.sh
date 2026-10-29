#!/usr/bin/env bash
# Launch the FlexPPI performance-first v2 five-fold training inside one tmux
# session (one window per fold).  Idempotent: skips folds that already started.
set -euo pipefail

REPO=/ranze/repos/FlexPPI
VENV="$REPO/.venv-linux/bin"
SESSION=flexppi_v2
LOGDIR="$REPO/artifacts/performance_first_v2/logs"
mkdir -p "$LOGDIR"

tmux has-session -t "$SESSION" 2>/dev/null || tmux new-session -d -s "$SESSION" -n status
tmux rename-window -t "$SESSION:0" status

for fold in 0 1 2 3 4; do
  LOG="$LOGDIR/fold_${fold}.log"
  if [ -f "$LOG" ]; then
    echo "fold ${fold}: log already exists at $LOG (skipping launch)"
    continue
  fi
  tmux new-window -t "$SESSION" -n "fold$fold" \
    "cd $REPO && PYTHONUNBUFFERED=1 $VENV/flexppi train --config configs/performance_first_v2_fold${fold}.yaml 2>&1 | tee $LOG"
  echo "fold ${fold}: launched in tmux window $SESSION:fold$fold"
done

echo "done. attach with: tmux attach -t $SESSION"
