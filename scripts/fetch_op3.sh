#!/usr/bin/env bash
# fetch_op3.sh -- OP3 PBMC drug-perturbation data (146 compounds x 3 donors x 4 cell types).
# Pulls the Open Problems benchmark build (~1.1 GB) rather than the 14.2 GB GEO h5ad: same
# experiment, already QC'd, pseudobulked and carrying SMILES per compound.
# GEO GSE279945 (raw) / Szalata et al. NeurIPS 2024 D&B. CC BY 4.0.
# Usage: scripts/fetch_op3.sh [outdir]       default: shared/perturbation/op3
# Sizes are checked; S3 ETags are multipart so they are not md5s. Content is checked by
# scripts/check_op3.py, which is what decides whether a file is reusable.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$ROOT/shared/perturbation/op3}"
BASE="https://openproblems-data.s3.amazonaws.com/resources/task_perturbation_prediction/datasets/neurips-2023-data"
mkdir -p "$OUT"

# file                        bytes        what it holds
while read -r name bytes _; do
  [ -z "$name" ] && continue
  dest="$OUT/$name"
  if [ -f "$dest" ] && [ "$(wc -c <"$dest")" -eq "$bytes" ]; then echo "ok $name (cached)"; continue; fi
  echo ">> $name ($bytes bytes)"
  curl -fL --retry 3 -C - -o "$dest" "$BASE/$name"
  have="$(wc -c <"$dest")"
  [ "$have" -eq "$bytes" ] || { echo "FAIL $name: $have bytes, expected $bytes" >&2; exit 1; }
  echo "ok $name"
done <<'EOF'
id_map.csv	3860	test-split compound x cell type index
de_train.h5ad	183168750	limma DE per (compound, cell type), train split; has SMILES + sm_lincs_id
de_test.h5ad	109139040	same, held-out B and Myeloid cells
sc_test.h5ad	182823707	single cells, test split
sc_train.h5ad	893316406	single cells, train split (largest; needed only for the single-cell model)
EOF
echo ">> done. Next: python scripts/check_op3.py $OUT"
