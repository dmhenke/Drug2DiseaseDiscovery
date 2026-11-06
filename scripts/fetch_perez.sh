#!/usr/bin/env bash
# Perez 2022 PBMC object (CELLxGENE dataset 218acb0f), 12,218,105,530 bytes per the curation API, checked 2026-11-06.
# Resumable. Skips the download when the file is already complete.
set -euo pipefail
cd "$(dirname "$0")/.."
out=data/cellxgene/perez_pbmc.h5ad
want=12218105530
mkdir -p "$(dirname "$out")"
have() { stat -c %s "$out" 2>/dev/null || stat -f %z "$out" 2>/dev/null || echo 0; }
if [ "$(have)" != "$want" ]; then
  df -h . | tail -1
  curl -fL -C - --retry 5 --retry-delay 10 -o "$out" \
    https://datasets.cellxgene.cziscience.com/c55dc602-d168-4d15-acc1-5de4f2f5d551.h5ad
fi
[ "$(have)" = "$want" ] || { echo "FAILED: $out is $(have) bytes, expected $want" >&2; exit 1; }
echo "ok: $out"
