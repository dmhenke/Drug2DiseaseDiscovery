#!/usr/bin/env bash
# op3_report.sh -- run the OP3 content checks and pack the evidence into one small tarball to carry
# off the node. Builds the `chem` conda env if it is missing. A failed check still produces a tarball,
# because the failure is the thing worth reading. Contains no h5ad data.
# Usage: scripts/op3_report.sh [datadir]     default: shared/perturbation/op3
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA="${1:-$ROOT/shared/perturbation/op3}"
STAMP="$(date +%Y%m%d-%H%M%S)"
OUT="$HOME/op3_report_$STAMP.tar.gz"
CHECK="$DATA/check"
command -v conda >/dev/null || { echo "conda not found on PATH" >&2; exit 1; }
rm -rf "$CHECK"   # only this script's own output; a stale report must never travel with a new one
mkdir -p "$CHECK"

{
  echo "date:    $(date)"
  echo "host:    $(hostname)"
  echo "commit:  $(git -C "$ROOT" rev-parse --short HEAD)"
  echo "dirty:   $(git -C "$ROOT" status --porcelain | wc -l | tr -d ' ') changed files (should be 0)"
  echo "data:    $DATA"
  ls -l "$DATA"/*.h5ad "$DATA"/id_map.csv
  conda env list | grep -q '^chem ' || {
    echo ">> creating conda env chem"
    conda create -y -n chem python=3.11 rdkit anndata scikit-learn pandas -c conda-forge
  }
  conda run --no-capture-output -n chem python -c "
import sys, numpy, pandas, anndata, rdkit
print('python', sys.version.split()[0], 'numpy', numpy.__version__, 'pandas', pandas.__version__,
      'anndata', anndata.__version__, 'rdkit', rdkit.__version__)"
  conda run --no-capture-output -n chem python "$ROOT/scripts/check_op3.py" "$DATA"
  echo "check exit code: $?"
} 2>&1 | tee "$CHECK/run.log"

tar -czf "$OUT" -C "$DATA" check || { echo "FAILED to write $OUT" >&2; exit 1; }
echo
grep -h "OP3 CHECK" "$CHECK/run.log" || echo "OP3 CHECK: no verdict, the check did not finish; see run.log"
echo "Report: $OUT ($(wc -c <"$OUT" | tr -d ' ') bytes)"
tar -tzf "$OUT"
