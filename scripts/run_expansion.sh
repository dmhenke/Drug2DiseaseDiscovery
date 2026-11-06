#!/usr/bin/env bash
# run_expansion.sh -- add cohorts to the SLE signature, then pack the evidence into one tarball.
#   1 Perez 2022 cohort 4, every cell, author labels (discovery)
#   2 GSE135779 children and adults, marker-based labels (replication, each analysed alone)
#   3 replication tables and a consensus signature; compounds rescored on it if OP3 is present
# PCA stand-in throughout: with labels set, the DE does not depend on the embedding, and no GPU is needed.
# A failed step is logged and the rest still run. Finished studies are skipped; set FORCE=1 to redo them.
# Usage: scripts/run_expansion.sh     Env: SRC (Perez h5ad), FORCE, DISCOVERY_ONLY=1 (step 1 only)
set -uo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"; cd "$ROOT"
SRC="${SRC:-$ROOT/data/cellxgene/perez_pbmc.h5ad}"
STAMP="$(date +%Y%m%d-%H%M%S)"
WORK="$ROOT/analysis/expansion_$STAMP"; mkdir -p "$WORK"
LOG="$WORK/run_log.txt"
ENV=scrna

step() { echo; echo "===== $1 | $(date '+%F %T')"; }
fail() { echo "FAILED: $1"; FAILED="${FAILED:-} [$1]"; }

run_de() {  # run_de <study>
  local d="studies/$1"
  if [ -z "${FORCE:-}" ] && [ -s "$d/results/DE_all_human.csv" ] && [ -f "$d/results/.complete" ]; then
    echo "$1: already complete, skipping"; return 0
  fi
  rm -f "$d/results/.complete" "$d/results/DE_"*.csv
  UCE_MODE=pca UCE_PYTHON="$PY" scripts/run_study.sh "$d" full && touch "$d/results/.complete"
}

main() {
  echo "commit $(git rev-parse --short HEAD); $(git status --porcelain | wc -l | tr -d ' ') changed files; host $(hostname)"
  git status --porcelain
  df -h "$ROOT" | tail -1; free -g 2>/dev/null | sed -n 2p
  command -v conda >/dev/null || { echo "conda not on PATH"; exit 1; }
  conda env list | grep -q "^$ENV " || conda create -y -n "$ENV" python=3.11 scanpy python-igraph leidenalg -c conda-forge \
    || { echo "could not create conda env $ENV"; exit 1; }
  PY="$(conda run -n "$ENV" which python)"
  "$PY" -c "import scanpy, leidenalg, scipy; print('scanpy', scanpy.__version__)" || { echo "env $ENV is broken"; exit 1; }

  step "1 Perez cohort 4, full"
  if [ ! -f "$SRC" ]; then
    fail "Perez source missing: $SRC (set SRC=...)"
  else
    if [ ! -f studies/lupus_c4_full/data/lupus_c4_full.h5ad ]; then
      "$PY" scripts/subsample_h5ad.py "$SRC" studies/lupus_c4_full/data/lupus_c4_full.h5ad \
        --donor-col donor_id --condition-col disease --healthy normal --diseased "systemic lupus erythematosus" \
        --all-donors --max-cells-per-donor 0 --counts raw --where sex=female \
        --where "self_reported_ethnicity=European American" --where Processing_Cohort=4.0 || fail "Perez subset"
    fi
    run_de lupus_c4_full || fail "lupus_c4_full DE"
  fi

  if [ -n "${DISCOVERY_ONLY:-}" ]; then echo "DISCOVERY_ONLY: stopping after step 1"; echo "FAILED STEPS:${FAILED:- none}"; return; fi

  step "2a marker annotation against author labels (20 v 20 test subset)"
  if [ -f studies/lupus_test/data/lupus_c4.h5ad ]; then
    "$PY" scripts/annotate_markers.py studies/lupus_test/data/lupus_c4.h5ad "$WORK/annotation_check" \
      --truth cell_type --symbols feature_name || fail "annotation check"
  else
    fail "annotation check: studies/lupus_test/data/lupus_c4.h5ad missing"
  fi

  step "2b GSE135779 download and build"
  if [ ! -f studies/gse135779_child/data/gse135779_child.h5ad ] || [ ! -f studies/gse135779_adult/data/gse135779_adult.h5ad ]; then
    scripts/fetch_gse135779.sh || fail "GSE135779 download"
    "$PY" scripts/build_gse135779.py shared/gse135779 studies || fail "GSE135779 build"
  fi
  run_de gse135779_child || fail "gse135779_child DE"
  run_de gse135779_adult || fail "gse135779_adult DE"

  step "3 replication"
  local disc=lupus_c4_full
  if "$PY" scripts/compare_cohorts.py "$WORK/replication" "$disc" gse135779_child gse135779_adult; then
    if [ -f shared/perturbation/op3/de_train.h5ad ] && conda env list | grep -q '^chem '; then
      conda run --no-capture-output -n chem python scripts/op3_rescue.py shared/perturbation/op3 \
        "$WORK/replication/consensus_signatures.tsv" "$WORK/rescue_consensus" || fail "rescue on consensus"
    else
      echo "rescue skipped: OP3 files or the chem env are absent"
    fi
  else
    fail "replication"
  fi

  step "packing"
  for s in lupus_c4_full gse135779_child gse135779_adult; do
    mkdir -p "$WORK/$s"
    cp studies/"$s"/results/*.csv studies/"$s"/results/*.tsv studies/"$s"/results/warnings.txt studies/"$s"/results/run_*.log "$WORK/$s/" 2>/dev/null
  done
  echo "FAILED STEPS:${FAILED:- none}"
}

main 2>&1 | tee "$LOG"
OUT="$HOME/expansion_report_$STAMP.tar.gz"
tar -czf "$OUT" -C "$ROOT/analysis" "expansion_$STAMP" || { echo "could not write $OUT"; exit 1; }
echo; grep -E "^(FAILED STEPS|ANNOTATION CHECK|  \[(ok|FAIL)\])" "$LOG"
echo "Report: $OUT ($(wc -c <"$OUT" | tr -d ' ') bytes)"
