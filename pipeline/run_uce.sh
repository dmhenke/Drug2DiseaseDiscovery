#!/usr/bin/env bash
# run_uce.sh -- Stage 1: embed each species file of a study.
#   UCE_MODE=uce (default): UCE embedding (UCE_NLAYERS, default 33). Needs CUDA or MPS.
#   UCE_MODE=pca:           PCA stand-in (scripts/fake_uce.py). NOT UCE; plumbing and baseline only.
# Usage: pipeline/run_uce.sh studies/<name>
# Env:   UCE_MODE, UCE_PYTHON, UCE_REPO, UCE_MODEL, UCE_NLAYERS (33), BATCH (25),
#        SHARD_CELLS (unset = embed each file whole), UCE_SEED (1234)
# Output: <study>/interim/uce/<study>_<species>.h5ad with obsm["X_uce"]
set -euo pipefail

STUDY_DIR="${1:?Usage: run_uce.sh <study_dir>}"
STUDY_DIR="$(cd "$STUDY_DIR" && pwd)"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
STUDY_NAME="$(basename "$STUDY_DIR")"
MODE="${UCE_MODE:-uce}"
PY="${UCE_PYTHON:-python3}"
UCE_REPO="${UCE_REPO:-$ROOT/pipeline/UCE}"
MODEL="${UCE_MODEL:-$ROOT/shared/models/33l_8ep_1024t_1280.torch}"
NLAYERS="${UCE_NLAYERS:-33}"
BATCH="${BATCH:-25}"
OUTDIR="$STUDY_DIR/interim/uce/"   # trailing slash: UCE concatenates paths, it does not join them
SHARDS="$STUDY_DIR/interim/shards/"

mkdir -p "$OUTDIR"

if [ "$MODE" = pca ]; then
  for SP in human mouse; do
    IN="$STUDY_DIR/interim/${STUDY_NAME}_${SP}.h5ad"
    [ -f "$IN" ] || continue
    echo ">> PCA stand-in (NOT UCE) for $SP: $IN"
    "$PY" "$ROOT/scripts/fake_uce.py" "$IN" "$OUTDIR${STUDY_NAME}_${SP}.h5ad"
  done
  echo ">> done (pca mode)."
  exit 0
fi

[ -f "$MODEL" ] || { echo "ERROR: weights not found: $MODEL" >&2; exit 1; }

# UCE reuses *_proc.h5ad and *_counts.npz if present, without checking their content
# (report finding B4). Clear them so every run is fresh.
rm -f "$OUTDIR"*_proc.h5ad "$OUTDIR"*_counts.npz "$OUTDIR"*_shapes_dict.pkl \
      "$OUTDIR"*_chroms.pkl "$OUTDIR"*_starts.pkl "$OUTDIR"*_pe_idx.torch \
      "$OUTDIR"*_uce_adata.h5ad
if [ -n "${SHARD_CELLS:-}" ]; then
  rm -rf "$SHARDS"
  mkdir -p "$SHARDS"
fi

embed() {  # embed <input h5ad> <species>
  (cd "$UCE_REPO" && "$PY" "$ROOT/scripts/uce_seeded.py" \
      --adata_path "$1" --dir "$OUTDIR" --species "$2" \
      --model_loc "$MODEL" --nlayers "$NLAYERS" --batch_size "$BATCH")
}

for SP in human mouse; do
  IN="$STUDY_DIR/interim/${STUDY_NAME}_${SP}.h5ad"
  [ -f "$IN" ] || { echo ">> no $SP input; skipping"; continue; }
  if [ -n "${SHARD_CELLS:-}" ]; then
    "$PY" "$ROOT/scripts/split_shards.py" "$IN" "$SHARDS" "$SHARD_CELLS"
    for S in "$SHARDS${STUDY_NAME}_${SP}"_s[0-9][0-9].h5ad; do
      echo ">> embedding shard $(basename "$S") ($NLAYERS-layer, batch $BATCH)"
      embed "$S" "$SP"
    done
    "$PY" "$ROOT/scripts/merge_uce_shards.py" "$IN" "$OUTDIR" "$OUTDIR${STUDY_NAME}_${SP}.h5ad"
  else
    echo ">> embedding $SP ($NLAYERS-layer, batch $BATCH): $IN"
    embed "$IN" "$SP"
    mv "$OUTDIR${STUDY_NAME}_${SP}_uce_adata.h5ad" "$OUTDIR${STUDY_NAME}_${SP}.h5ad"
  fi
  [ -f "$OUTDIR${STUDY_NAME}_${SP}.h5ad" ] || { echo "ERROR: no output for $SP" >&2; exit 1; }
done
echo ">> done. Next: Stage 2."
