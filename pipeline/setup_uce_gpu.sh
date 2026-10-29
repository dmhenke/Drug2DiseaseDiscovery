#!/usr/bin/env bash
# ===========================================================================
# setup_uce_gpu.sh  --  one-time setup for the UCE embedding step (Stage 1).
# Run once on the GPU machine; shared across ALL studies (lung, kidney, ...).
# Clones UCE into pipeline/UCE and stores model weights in shared/models/.
# ===========================================================================
set -euo pipefail
PROJECT_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# 0. GPU + driver visible? (preinstalled on cloud A100/H100 images)
nvidia-smi

# 1. Isolated Python env (Python 3.10)
conda env list | grep -q '^uce ' || conda create -y -n uce python=3.10
# conda's activate scripts read unset vars (PS1), so relax -u around them
set +u
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate uce
set -u
[ "$(command -v python)" = "$(conda info --base)/envs/uce/bin/python" ] || { echo "conda activate uce failed" >&2; exit 1; }

# 2. CUDA-matched PyTorch (match the CUDA version nvidia-smi prints; cu121 shown)
pip install torch==2.1.1 --index-url https://download.pytorch.org/whl/cu121
python -c "import torch; print('CUDA available:', torch.cuda.is_available())"   # must be True

# 3. Clone UCE into pipeline/UCE and install deps
cd "$PROJECT_ROOT/pipeline"
[ -d UCE ] || git clone https://github.com/snap-stanford/UCE.git
cd UCE && pip install -r requirements.txt

# 4. 33-layer weights + UCE support files (~11.4 GB), md5-verified. Not UCE's own downloader:
#    figshare.com answers scripts with an empty 202 and UCE would cache the empty files.
"$PROJECT_ROOT/scripts/fetch_uce_files.sh"

# 5. Accelerate config (single-GPU defaults are fine)
accelerate config default || accelerate config

# 6. No demo smoke test: with no args UCE fetches the 4-layer model and demo data through the
#    blocked figshare host. scripts/preflight.sh + scripts/benchmark_uce.sh are the smoke test.

echo ">> Setup complete. Embed a study with:  pipeline/run_uce.sh studies/<name>"
echo ">> If CUDA out-of-memory: lower BATCH (e.g. BATCH=8 pipeline/run_uce.sh studies/lung)"
