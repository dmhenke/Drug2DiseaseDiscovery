#!/usr/bin/env bash
# LINCS Level 5 (replicate-consensus z-scores). Run only after scripts/lincs_triage.py says a phase is usable.
#   phase2 (default)  5.0 GB gzipped, ~5.8 GB unzipped, 118,050 signatures   <- start here
#   phase1            20 GB gzipped, ~23 GB unzipped, 473,647 signatures
# Both are resumable. The unzipped .gctx is kept (HDF5 cannot be sliced while gzipped) and the .gz removed.
# Sizes are the GEO FTP listing of 2026-11-06.
set -euo pipefail
cd "$(dirname "$0")/.."
out=shared/lincs
mkdir -p "$out"

case "${1:-phase2}" in
  phase2) url=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE70nnn/GSE70138/suppl/GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz ;;
  phase1) url=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92742/suppl/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz ;;
  *) echo "usage: $0 [phase2|phase1]" >&2; exit 2 ;;
esac

gz="$out/$(basename "$url")"
gctx="${gz%.gz}"
if [ -s "$gctx" ]; then echo "already unzipped: $gctx"; else
  df -h "$out" | tail -1
  curl -fL -C - --retry 5 --retry-delay 10 -o "$gz" "$url"
  gzip -t "$gz"                      # content test before spending disk on the unzip
  gunzip -k "$gz" && rm -f "$gz"
fi
python3 - "$gctx" <<'EOF'
import sys, h5py
with h5py.File(sys.argv[1], "r") as f:   # content test: a gctx has these three datasets or it is not usable
    m = f["0/DATA/0/matrix"]
    print("matrix", m.shape, "rows", f["0/META/ROW/id"].shape[0], "cols", f["0/META/COL/id"].shape[0])
EOF
echo "ok: $gctx"
