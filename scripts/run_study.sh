#!/usr/bin/env bash
# run_study.sh -- preflight, then Stages 0-3 for one study directory, logged with a manifest.
# Usage: scripts/run_study.sh studies/<name> [diagnostic|full]
# Env:   UCE_MODE (uce|pca), UCE_PYTHON, UCE_MODEL, UCE_NLAYERS, BATCH, SHARD_CELLS,
#        UCE_SEED, CHECK (a script to run at the end, e.g. scripts/check_t0.py)
# Output: <study>/results/run_<kind>_<timestamp>.log. The kind is on the first and last lines.
# Keep test and real runs in separate study folders (e.g. lupus_test, lupus_full).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
STUDY_DIR="$(cd "${1:?usage: run_study.sh studies/<name> [diagnostic|full]}" && pwd)"
KIND="${2:-diagnostic}"
case "$KIND" in diagnostic|full) ;; *) echo "kind must be diagnostic or full" >&2; exit 1 ;; esac
MODE="${UCE_MODE:-uce}"
PY="${UCE_PYTHON:-python3}"
mkdir -p "$STUDY_DIR/results"
LOG="$STUDY_DIR/results/run_${KIND}_$(date +%Y%m%d-%H%M%S).log"

main() {
  echo "===== RUN KIND: $KIND | study: $(basename "$STUDY_DIR") | UCE_MODE: $MODE ====="
  echo "host: $(uname -srm) | started: $(date '+%Y-%m-%dT%H:%M:%S%z')"
  echo "code: $(cd "$ROOT" && cat pipeline/*.R pipeline/*.sh scripts/*.py scripts/*.sh | shasum -a 256 | cut -c1-16)"
  echo "params: BATCH=${BATCH:-25} UCE_NLAYERS=${UCE_NLAYERS:-33} SHARD_CELLS=${SHARD_CELLS:-none} UCE_SEED=${UCE_SEED:-1234}"
  if git -C "$ROOT" rev-parse --git-dir >/dev/null 2>&1; then
    local dirty; dirty="$(git -C "$ROOT" status --porcelain --untracked-files=no -- pipeline scripts studies)"
    echo "git: $(git -C "$ROOT" rev-parse --short HEAD)${dirty:+ DIRTY (uncommitted changes below)}"
    if [ -n "$dirty" ]; then echo "WARNING: dirty tree; numbers from this run do not match any commit"; echo "$dirty"; fi
  else
    echo "git: not a repository, so no dirty-tree check; code hash above is the version"
  fi
  "$ROOT/scripts/preflight.sh" "$STUDY_DIR"
  echo "--- stage 0"
  Rscript "$ROOT/pipeline/01_build_species_anndata.R" "$STUDY_DIR"
  echo "--- stage 1 ($MODE)"
  "$ROOT/pipeline/run_uce.sh" "$STUDY_DIR"
  echo "--- stage 2"
  Rscript "$ROOT/pipeline/02_integrate_and_cluster.R" "$STUDY_DIR"
  echo "--- stage 3"
  Rscript "$ROOT/pipeline/03_pseudobulk_de.R" "$STUDY_DIR"
  if [ -n "${CHECK:-}" ]; then
    echo "--- check: $CHECK"
    "$PY" "$ROOT/$CHECK" "$STUDY_DIR"
  fi
  echo "===== END RUN KIND: $KIND | finished: $(date '+%Y-%m-%dT%H:%M:%S%z') ====="
}

main 2>&1 | tee -a "$LOG"
