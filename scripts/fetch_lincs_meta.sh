#!/usr/bin/env bash
# LINCS L1000 metadata only (~15 MB): enough to decide which Level 5 matrix is worth downloading.
# Sizes are the GEO FTP listing of 2026-11-06 (gzipped). No expression data is fetched here.
set -euo pipefail
cd "$(dirname "$0")/.."
out=shared/lincs/meta
mkdir -p "$out"

base1=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl
base2=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl
files=(
  "$base1/GSE92742_Broad_LINCS_sig_info.txt.gz"
  "$base1/GSE92742_Broad_LINCS_cell_info.txt.gz"
  "$base1/GSE92742_Broad_LINCS_pert_info.txt.gz"
  "$base1/GSE92742_Broad_LINCS_gene_info.txt.gz"
  "$base2/GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz"
  "$base2/GSE70138_Broad_LINCS_cell_info_2017-04-28.txt.gz"
  "$base2/GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz"
  "$base2/GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz"
)
for u in "${files[@]}"; do
  f="$out/$(basename "$u")"
  if gzip -t "$f" 2>/dev/null; then echo "complete: $f"; continue; fi
  curl -fL -C - --retry 3 -o "$f" "$u"
  gzip -t "$f"   # content test: a truncated or HTML error page fails here, not three steps later
done
ls -lh "$out"
echo "ok: metadata in $out"
