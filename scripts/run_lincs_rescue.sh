#!/usr/bin/env bash
# Score LINCS compounds against one study's DE on this machine, then bundle the evidence. Downloads nothing: a missing
# input is named as a FAIL with the command that creates it. A tarball is written whatever happens.
#   24 h is primary because OP3 is a 24 h assay; all times is the sensitivity run.
# Usage: scripts/run_lincs_rescue.sh [study]    default lupus_c4_full
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
STUDY="${1:-lupus_c4_full}"
STAMP="$(date +%Y%m%d-%H%M%S)"
NAME="lincs_rescue_$STAMP"; WORK="analysis/$NAME"; mkdir -p "$WORK"
DE="studies/$STUDY/results/DE_all_human.csv"
GCTX=shared/lincs/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx
OP3=shared/perturbation/op3
BAD=0
need() { if [ -e "$1" ]; then echo "  [ok]   $1"; else echo "  [FAIL] $1 is missing; create it with: $2"; BAD=1; fi; }

main() {
  echo "commit $(git rev-parse --short HEAD); $(git status --porcelain | wc -l | tr -d ' ') changed files; host $(hostname)"
  git status --porcelain
  need "$DE" "DISCOVERY_ONLY=1 scripts/run_expansion.sh"
  need "studies/$STUDY/results/.complete" "DISCOVERY_ONLY=1 scripts/run_expansion.sh"
  need shared/lincs/meta/GSE92742_Broad_LINCS_sig_info.txt.gz scripts/fetch_lincs_meta.sh
  need "$GCTX" "scripts/fetch_lincs_level5.sh phase1"
  need "$OP3/de_train.h5ad" scripts/fetch_op3.sh
  conda env list | grep -q '^chem ' && echo "  [ok]   conda env chem" || { echo "  [FAIL] conda env chem is missing; create it with: scripts/op3_report.sh"; BAD=1; }
  if [ "$BAD" != 0 ]; then echo "RESULT: not run, inputs missing"; return; fi
  ls -l "$DE"; cp studies/"$STUDY"/results/run_*.log "$WORK/" 2>/dev/null
  PY=(conda run --no-capture-output -n chem python)
  "${PY[@]}" scripts/test_lincs_rescue.py || { echo "RESULT: not run, self-test failed"; return; }
  for spec in "t24 --time 24" "all_times"; do
    set -- $spec; out="$WORK/$1"; shift
    echo; echo "===== $out | $(date '+%F %T')"
    "${PY[@]}" scripts/lincs_rescue.py "$GCTX" shared/lincs/meta "$out" --de "$DE" --op3 "$OP3" "$@" \
      || echo "FAILED: $out"
  done
  echo "RESULT: finished"
}

main 2>&1 | tee "$WORK/run_log.txt"
OUT="$HOME/$NAME.tar.gz"
tar -czf "$OUT" --exclude=lincs_rescue.tsv -C analysis "$NAME" || { echo "could not write $OUT"; exit 1; }
echo "Report: $OUT ($(wc -c <"$OUT" | tr -d ' ') bytes). Full scores stay on the node in $WORK."
