#!/usr/bin/env bash
# benchmark_uce.sh -- measure UCE throughput and peak memory on THIS machine.
# Diagnostic only, not a result. Runs UCE on N and 2N cells; the marginal rate
# (N cells / extra seconds) removes the fixed model-load time.
# Usage: scripts/benchmark_uce.sh <interim_species.h5ad> <N> [BATCH]
# Output: bench/<timestamp>/summary.txt
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IN="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
N="${2:?usage: benchmark_uce.sh <species.h5ad> <N> [BATCH]}"
BATCH="${3:-${BATCH:-4}}"
PY="${UCE_PYTHON:-python3}"
UCE_REPO="${UCE_REPO:-$ROOT/pipeline/UCE}"
MODEL="${UCE_MODEL:-$ROOT/shared/models/33l_8ep_1024t_1280.torch}"
NLAYERS="${UCE_NLAYERS:-33}"
OUT="$ROOT/bench/$(date +%Y%m%d-%H%M%S)"
mkdir -p "$OUT"
SP=human; case "$(basename "$IN")" in *mouse*) SP=mouse ;; esac
if [ "$(uname -s)" = Darwin ]; then TIME=(/usr/bin/time -l); else TIME=(/usr/bin/time -v); fi
[ -x /usr/bin/time ] || TIME=()  # minimal Ubuntu lacks GNU time (apt install time); host RSS then unreported
trap 'kill $(jobs -p) 2>/dev/null || true' EXIT

{
  echo "device: $("$PY" -c 'from accelerate import Accelerator; print(Accelerator().device)')"
  echo "layers=$NLAYERS batch=$BATCH species=$SP input=$IN"
} | tee "$OUT/summary.txt"

for M in "$N" "$((2 * N))"; do
  SUB="$OUT/sub_${M}.h5ad"; RUN="$OUT/run_${M}"; mkdir -p "$RUN"
  "$PY" -c "import sys, anndata as ad; a = ad.read_h5ad(sys.argv[1], backed='r'); a[:int(sys.argv[3])].to_memory().write_h5ad(sys.argv[2])" "$IN" "$SUB" "$M"
  # GPU memory sampled every second; whole card, so other processes on it are included
  if command -v nvidia-smi >/dev/null; then
    nvidia-smi --query-gpu=memory.used,memory.total --format=csv,noheader,nounits -lms 1000 >"$RUN/gpu.txt" &
  fi
  START=$(date +%s)
  (cd "$UCE_REPO" && "${TIME[@]}" "$PY" "$ROOT/scripts/uce_seeded.py" \
      --adata_path "$SUB" --dir "$RUN/" --species "$SP" --model_loc "$MODEL" \
      --nlayers "$NLAYERS" --batch_size "$BATCH") >"$RUN/stdout.txt" 2>"$RUN/time.txt" \
    || { echo "FAILED at $M cells; see $RUN/time.txt and $RUN/stdout.txt" | tee -a "$OUT/summary.txt"; exit 1; }
  SECS=$(( $(date +%s) - START ))
  kill $(jobs -p) 2>/dev/null || true
  RSS=$(grep -i "maximum resident" "$RUN/time.txt" | head -1 | sed 's/^[[:space:]]*//')
  GPU=$( [ -s "$RUN/gpu.txt" ] && sort -t, -k1 -n "$RUN/gpu.txt" | tail -1 | awk -F', *' '{print $1 "/" $2 " MiB"}' || true )
  echo "cells=$M seconds=$SECS peak=[$RSS] gpu_peak=[${GPU:-n/a}]" | tee -a "$OUT/summary.txt"
done

"$PY" - "$OUT/summary.txt" <<'PYEOF' | tee -a "$OUT/summary.txt"
import re, sys
rows = [(int(re.search(r"cells=(\d+)", l).group(1)), int(re.search(r"seconds=(\d+)", l).group(1)))
        for l in open(sys.argv[1]) if l.startswith("cells=")]
(n1, t1), (n2, t2) = rows[0], rows[1]
if t2 <= t1:
    print("marginal rate: not measurable (2N ran no slower than N); raise N")
else:
    rate = (n2 - n1) / (t2 - t1)
    print(f"marginal rate: {rate:.2f} cells/s")
    for target in (10_000, 20_000, 1_263_676):
        print(f"  est. {target:>9,} cells: {target / rate / 3600:.1f} h")
PYEOF
