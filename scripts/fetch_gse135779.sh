#!/usr/bin/env bash
# fetch_gse135779.sh -- Nehar-Belaid 2020 SLE PBMC count matrices from GEO (1.3 GB tar, 56 samples).
# Usage: scripts/fetch_gse135779.sh [outdir]       default: shared/gse135779
# Resumable, size-checked. Content is checked by scripts/build_gse135779.py.
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$ROOT/shared/gse135779}"
BASE="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE135nnn/GSE135779/suppl"
BYTES=1299783680
mkdir -p "$OUT/raw"
TAR="$OUT/GSE135779_RAW.tar"
if [ ! -f "$TAR" ] || [ "$(wc -c <"$TAR")" -ne "$BYTES" ]; then
  curl -fL --retry 3 -C - -o "$TAR" "$BASE/GSE135779_RAW.tar"
fi
[ "$(wc -c <"$TAR")" -eq "$BYTES" ] || { echo "FAIL: $TAR is not $BYTES bytes" >&2; exit 1; }
curl -fL --retry 3 -o "$OUT/GSE135779_genes.tsv.gz" "$BASE/GSE135779_genes.tsv.gz"
tar -xf "$TAR" -C "$OUT/raw"
N="$(find "$OUT/raw" -name '*.gz' | wc -l)"
[ "$N" -eq 112 ] || { echo "FAIL: expected 112 files (56 samples x 2), found $N" >&2; exit 1; }
echo "ok: $N files in $OUT/raw"
