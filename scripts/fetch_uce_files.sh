#!/usr/bin/env bash
# Download the UCE 33-layer weights and support files, verify md5, extract protein embeddings.
# Uses ndownloader.figshare.com: figshare.com/ndownloader answers scripts with an AWS WAF
# challenge (HTTP 202, empty body), which wget/curl save as a 0-byte file and UCE's own
# downloader (requests, no status check) saves as an empty file it then never re-fetches.
# ~11.4 GB total. Usage: scripts/fetch_uce_files.sh   Env: UCE_MODEL, UCE_REPO
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MODEL="${UCE_MODEL:-$ROOT/shared/models/33l_8ep_1024t_1280.torch}"
REPO="${UCE_REPO:-$ROOT/pipeline/UCE}"
md5() { if command -v md5sum >/dev/null; then md5sum "$1" | cut -d' ' -f1; else md5 -q "$1"; fi; }
grep -v '^#' "$ROOT/scripts/uce_files.tsv" | while IFS=$'\t' read -r id where name bytes sum; do
  if [ "$where" = model ]; then out="$MODEL"; else out="$REPO/model_files/$name"; fi
  mkdir -p "$(dirname "$out")"
  if [ -f "$out" ] && [ "$(md5 "$out")" = "$sum" ]; then echo "ok (cached) $name"; continue; fi
  # a wrong-size file cannot be resumed (e.g. a 0-byte WAF response, or a larger corrupt file)
  if [ -f "$out" ] && [ "$(wc -c <"$out" | tr -d ' ')" -ge "$bytes" ]; then rm -f "$out"; fi
  echo ">> $name ($bytes bytes)"
  curl -fL --retry 3 -C - -o "$out" "https://ndownloader.figshare.com/files/$id"
  got="$(md5 "$out")"
  if [ "$got" != "$sum" ]; then rm -f "$out"; echo "FAIL: md5 mismatch for $name ($got != $sum); removed, rerun" >&2; exit 1; fi
  echo "ok $name"
  if [[ "$name" == *.tar.gz ]]; then tar -xzf "$out" -C "$(dirname "$out")"; echo "   extracted into $(dirname "$out")"; fi
done
