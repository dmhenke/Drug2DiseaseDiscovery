#!/usr/bin/env bash
# Build studies/lupus_test/data/lupus_c<cohort>.h5ad from the Perez 2022 PBMC object (CELLxGENE c55dc602).
# Female, European American; one file per Processing_Cohort (cohort 1 is healthy-only, excluded).
# Default is cohort 4 alone: one batch, so no batch term is needed and no donor can be counted twice.
# With several cohorts, donors sequenced in more than one are dropped -- which removes 20 of the
# 21 healthy donors in cohort 2 (dry run 2026-10-27), so multi-cohort needs a different dedup.
# Usage: scripts/make_lupus_test.sh [--dry-run]   Env: SRC, DONORS (5), CELLS (600), COHORTS ("4")
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${SRC:-$ROOT/data/cellxgene/perez_pbmc.h5ad}"
COHORTS="${COHORTS:-4}"
SPAN=(); [ "$(wc -w <<<"$COHORTS")" -gt 1 ] && SPAN=(--drop-donors-spanning Processing_Cohort)
for k in $COHORTS; do
  echo "== cohort $k"
  python3 "$ROOT/scripts/subsample_h5ad.py" "$SRC" "$ROOT/studies/lupus_test/data/lupus_c$k.h5ad" \
    --donor-col donor_id --condition-col disease \
    --healthy normal --diseased "systemic lupus erythematosus" \
    --donors-per-condition "${DONORS:-5}" --max-cells-per-donor "${CELLS:-600}" --counts raw \
    ${SPAN[@]+"${SPAN[@]}"} \
    --where sex=female --where "self_reported_ethnicity=European American" \
    --where "Processing_Cohort=$k.0" "$@"
done
