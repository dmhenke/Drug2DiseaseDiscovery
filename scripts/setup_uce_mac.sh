#!/usr/bin/env bash
# setup_uce_mac.sh -- Python env for UCE on Apple Silicon (no CUDA). Idempotent.
# Creates $UCE_VENV (default ~/uce-venv) with UCE's pinned requirements.
# Needs python3.10 or python3.11 (UCE pins torch 2.1.1).
# UNTESTED on this machine: the pins predate current macOS builds. preflight.sh checks MPS.
# Note: this clones UCE from GitHub and installs packages (external downloads; get approval first).
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
VENV="${UCE_VENV:-$HOME/uce-venv}"
PYBIN="$(command -v python3.10 || command -v python3.11 || true)"
[ -n "$PYBIN" ] || { echo "Need python3.10 or 3.11 (e.g. brew install python@3.10)" >&2; exit 1; }
if [ ! -d "$ROOT/pipeline/UCE" ]; then
  git clone https://github.com/snap-stanford/UCE.git "$ROOT/pipeline/UCE"
fi
[ -d "$VENV" ] || "$PYBIN" -m venv "$VENV"
"$VENV/bin/pip" install --upgrade pip
"$VENV/bin/pip" install torch==2.1.1
"$VENV/bin/pip" install -r "$ROOT/pipeline/UCE/requirements.txt"
echo ">> done. Then: export UCE_PYTHON=$VENV/bin/python"
