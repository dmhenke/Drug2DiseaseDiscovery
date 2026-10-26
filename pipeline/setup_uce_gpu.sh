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
conda create -y -n uce python=3.10
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"; conda activate uce

# 2. CUDA-matched PyTorch (match the CUDA version nvidia-smi prints; cu121 shown)
pip install torch==2.1.1 --index-url https://download.pytorch.org/whl/cu121
python -c "import torch; print('CUDA available:', torch.cuda.is_available())"   # must be True

# 3. Clone UCE into pipeline/UCE and install deps
cd "$PROJECT_ROOT/pipeline"
[ -d UCE ] || git clone https://github.com/snap-stanford/UCE.git
cd UCE && pip install -r requirements.txt

# 4. Put the 33-layer paper weights in the SHARED model dir (used by every study)
mkdir -p "$PROJECT_ROOT/shared/models"
echo ">> Download 33l_8ep_1024t_1280.torch into $PROJECT_ROOT/shared/models/"
echo "   figshare: https://figshare.com/articles/dataset/Universal_Cell_Embedding_Model_Files/24320806"

# 5. Accelerate config (single-GPU defaults are fine)
accelerate config default || accelerate config

# 6. Smoke test on the bundled demo (downloads ESM2 protein + gene tables)
python eval_single_anndata.py

echo ">> Setup complete. Embed a study with:  pipeline/run_uce.sh studies/<name>"
echo ">> If CUDA out-of-memory: lower BATCH (e.g. BATCH=8 pipeline/run_uce.sh studies/lung)"
