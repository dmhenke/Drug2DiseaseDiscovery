#!/usr/bin/env bash
# Build studies/lupus_test/data/lupus_c{2,3,4}.h5ad from the Perez 2022 PBMC object (CELLxGENE c55dc602).
# Female, European American; one file per Processing_Cohort (cohort 1 is healthy-only, excluded).
# Donors sequenced in more than one cohort are dropped so no donor is counted twice.
# Usage: scripts/make_lupus_test.sh [--dry-run]   Env: SRC, DONORS (3), CELLS (600), COHORTS ("2 3 4")
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SRC="${SRC:-$ROOT/data/cellxgene/perez_pbmc.h5ad}"
for k in ${COHORTS:-2 3 4}; do
  echo "== cohort $k"
  python3 "$ROOT/scripts/subsample_h5ad.py" "$SRC" "$ROOT/studies/lupus_test/data/lupus_c$k.h5ad" \
    --donor-col donor_id --condition-col disease \
    --healthy normal --diseased "systemic lupus erythematosus" \
    --donors-per-condition "${DONORS:-3}" --max-cells-per-donor "${CELLS:-600}" --counts raw \
    --drop-donors-spanning Processing_Cohort \
    --where sex=female --where "self_reported_ethnicity=European American" \
    --where "Processing_Cohort=$k.0" "$@"
done
