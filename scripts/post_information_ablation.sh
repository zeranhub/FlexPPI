#!/usr/bin/env bash
# Evaluate and export one information-ablation arm after its training finishes.
# Step 1: epoch-40 checkpoint of every completed fold is evaluated on its own
#         untouched test fold (writes artifacts/information_ablation/<arm>/evaluation/).
# Step 2: the provenance-checked OOF exporter runs with --verify against those
#         server evaluation JSONs and writes output/information_ablation_oof/<arm>/.
set -euo pipefail

VARIANT="${1:?usage: post_information_ablation.sh sequence|geometry|intrachain}"
case "$VARIANT" in
  sequence|geometry|intrachain) ;;
  *) echo "unknown variant: $VARIANT" >&2; exit 1 ;;
esac

REPO=/ranze/repos/FlexPPI
VENV="$REPO/.venv-linux/bin"
GPU="${FLEXPPI_GPU:-2}"
BASE="$REPO/artifacts/information_ablation/$VARIANT"
EVALDIR="$BASE/evaluation"
mkdir -p "$EVALDIR"
cd "$REPO"

export LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libcuda.so.1
unset _CUDA_COMPAT_PATH 2>/dev/null || true

for fold in 0 1 2 3 4; do
  CKPT="$BASE/fold_${fold}/epoch_040.pt"
  TARGET="$EVALDIR/fold_${fold}_test.json"
  if [ ! -f "$CKPT" ]; then
    echo "fold ${fold}: checkpoint missing (training not finished)"
    continue
  fi
  if [ -f "$TARGET" ]; then
    echo "fold ${fold}: evaluation already exists"
    continue
  fi
  CUDA_VISIBLE_DEVICES="$GPU" PYTORCH_CUDA_ALLOC_CONF=expandable_segments:True \
    "$VENV/flexppi" evaluate --config "configs/information_ablation_${VARIANT}_fold${fold}.yaml" \
    --checkpoint "$CKPT" --split test --bootstrap 200 --output "$TARGET" \
    >> "$EVALDIR/fold_${fold}_eval.log" 2>&1
  echo "fold ${fold}: wrote $TARGET"
done

if ! ls "$EVALDIR"/fold_*_test.json >/dev/null 2>&1; then
  echo "no evaluations written; stopping" >&2
  exit 1
fi

CUDA_VISIBLE_DEVICES="$GPU" "$VENV/python" scripts/export_information_ablation_oof.py \
  --variant "$VARIANT" --device cuda --threads 4 --verify \
  >> "$REPO/output/information_ablation_oof_${VARIANT}.log" 2>&1
echo "export complete for $VARIANT"
