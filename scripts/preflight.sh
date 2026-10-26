#!/usr/bin/env bash
# preflight.sh -- read-only checks before a run. Exits 1 on any FAIL.
# Usage: scripts/preflight.sh [study_dir]
# Env:   UCE_MODE=uce|pca (default uce), UCE_PYTHON, UCE_MODEL, UCE_REPO
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODE="${UCE_MODE:-uce}"; PY="${UCE_PYTHON:-python3}"
FAIL=0
ok()   { echo "  [ok]   $*"; }
warn() { echo "  [warn] $*"; }
bad()  { echo "  [FAIL] $*"; FAIL=1; }

echo "== host"
if [ "$(uname -s)" = Darwin ]; then
  RAM_GB=$(( $(sysctl -n hw.memsize) / 1073741824 ))
  echo "  $(sysctl -n machdep.cpu.brand_string), ${RAM_GB} GB RAM"
else
  RAM_GB=$(awk '/MemTotal/ {print int($2/1048576)}' /proc/meminfo)
  echo "  $(uname -m), ${RAM_GB} GB RAM, GPU: $(nvidia-smi --query-gpu=name,memory.total --format=csv,noheader 2>/dev/null || echo none)"
fi
FREE_GB=$(df -Pk "$ROOT" | awk 'NR==2 {print int($4/1048576)}')
if [ "$FREE_GB" -ge 20 ]; then ok "disk: ${FREE_GB} GB free"; else bad "disk: ${FREE_GB} GB free; need >= 20 GB (more for real data)"; fi

echo "== R"
if command -v Rscript >/dev/null 2>&1; then
  MISSING=$(Rscript -e 'need <- c("zellkonverter","SingleCellExperiment","S4Vectors","scuttle","DESeq2","bluster","BiocNeighbors","RcppAnnoy","igraph","Matrix"); cat(setdiff(need, rownames(installed.packages())))' 2>/dev/null)
  if [ -z "$MISSING" ]; then ok "R packages present"; else bad "R packages missing: $MISSING"; fi
else
  bad "Rscript not found"
fi

echo "== python (UCE_MODE=$MODE, $PY)"
if ! command -v "$PY" >/dev/null 2>&1; then
  bad "python not found: $PY"
elif [ "$MODE" = pca ]; then
  if "$PY" -c "import scanpy, anndata, pandas" >/dev/null 2>&1; then ok "scanpy, anndata, pandas importable"
  else bad "scanpy/anndata/pandas missing for $PY"; fi
else
  if "$PY" -c "import torch, accelerate" >/dev/null 2>&1; then ok "torch, accelerate importable"
  else bad "torch/accelerate missing for $PY (see scripts/setup_uce_mac.sh or pipeline/setup_uce_gpu.sh)"; fi
  DEV=$("$PY" -c "import torch; print('cuda' if torch.cuda.is_available() else 'mps' if torch.backends.mps.is_available() else 'cpu')" 2>/dev/null || echo unknown)
  case "$DEV" in
    cuda|mps) ok "accelerator: $DEV" ;;
    cpu)      bad "no CUDA or MPS: the 33-layer model is not practical on CPU" ;;
    *)        bad "could not query torch device" ;;
  esac
  MODEL="${UCE_MODEL:-$ROOT/shared/models/33l_8ep_1024t_1280.torch}"
  if [ -f "$MODEL" ]; then ok "weights: $MODEL"; else bad "weights missing: $MODEL (download needs approval)"; fi
  REPO="${UCE_REPO:-$ROOT/pipeline/UCE}"
  if [ -f "$REPO/eval_single_anndata.py" ]; then ok "UCE repo: $REPO"; else bad "UCE repo missing: $REPO"; fi
  if [ "$RAM_GB" -lt 24 ]; then warn "${RAM_GB} GB RAM: model load may peak near 10-12 GB (inferred); close other apps; use BATCH=4"; fi
fi

if [ -n "${1:-}" ]; then
  echo "== study $1"
  if [ -f "$1/config.R" ]; then ok "config.R"; else bad "config.R missing in $1"; fi
  if ls "$1"/data/*.h5ad >/dev/null 2>&1; then ok "data/*.h5ad present"; else bad "no data/*.h5ad in $1"; fi
fi

echo
if [ "$FAIL" -eq 0 ]; then echo "PREFLIGHT: PASS"; else echo "PREFLIGHT: FAIL"; fi
exit "$FAIL"
