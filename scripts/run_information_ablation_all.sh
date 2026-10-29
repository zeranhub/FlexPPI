#!/usr/bin/env bash
# Run all three information-ablation arms in parallel, one arm per GPU, and
# post-process each arm (epoch-40 evaluation + verified OOF export) as soon as
# its five folds finish.  Every five minutes the launcher is re-run: folds whose
# processes died before epoch_040.pt existed are restarted and resume from
# their last checkpoint.  Stops when all folds are trained and all arms are
# post-processed.
set -uo pipefail

REPO=/ranze/repos/FlexPPI
cd "$REPO"
GPU_SEQUENCE="${GPU_SEQUENCE:-2}"
GPU_GEOMETRY="${GPU_GEOMETRY:-0}"
GPU_INTRACHAIN="${GPU_INTRACHAIN:-1}"

post_arm() {
  local variant="$1"
  local marker="$REPO/artifacts/information_ablation/$variant/.post_complete"
  [ -f "$marker" ] && return 0
  local finished
  finished=$(ls artifacts/information_ablation/$variant/fold_*/epoch_040.pt 2>/dev/null | wc -l)
  [ "$finished" -ne 5 ] && return 0
  echo "$(date -u) post-processing $variant"
  if FLEXPPI_GPU="$GPU_SEQUENCE" bash scripts/post_information_ablation.sh "$variant"; then
    touch "$marker"
    echo "$(date -u) post-processing $variant complete"
  else
    echo "$(date -u) post-processing $variant failed; will retry next cycle"
  fi
}

while true; do
  FLEXPPI_GPU="$GPU_SEQUENCE" bash scripts/launch_information_ablation.sh sequence || true
  FLEXPPI_GPU="$GPU_GEOMETRY" bash scripts/launch_information_ablation.sh geometry || true
  FLEXPPI_GPU="$GPU_INTRACHAIN" bash scripts/launch_information_ablation.sh intrachain || true
  finished=$(ls artifacts/information_ablation/*/fold_*/epoch_040.pt 2>/dev/null | wc -l)
  echo "$(date -u) finished folds: $finished/15"
  for variant in sequence geometry intrachain; do
    post_arm "$variant"
  done
  done_markers=$(ls artifacts/information_ablation/*/.post_complete 2>/dev/null | wc -l)
  if [ "$finished" -eq 15 ] && [ "$done_markers" -eq 3 ]; then
    echo "=== all arms trained and post-processed $(date -u) ==="
    break
  fi
  sleep 300
done
