#!/usr/bin/env bash
# Stage 1: embed each species file with UCE (33-layer paper model).
# This is the ONLY Python step -- UCE has no R port. Everything before (01_*.R)
# and after (02_*.R, 03_*.R) is R.
#
# Prereqs:
#   git clone https://github.com/snap-stanford/UCE.git && cd UCE
#   pip install -r requirements.txt
#   # download 33-layer weights (README link) into model_files/:
#   #   https://figshare.com/articles/dataset/Universal_Cell_Embedding_Model_Files/24320806
#
# UCE embeds ONE species per run, so mouse and human are embedded separately and
# re-joined in R (Stage 2). Auxiliary files (ESM2 protein embeddings, gene tables)
# download automatically on first run. Needs a GPU.
#
# Outputs: uce_out/human_lung.h5ad and uce_out/mouse_lung.h5ad,
#          each with cell embeddings in .obsm["X_uce"].

set -euo pipefail

OUTDIR="uce_out"
MODEL="model_files/33l_8ep_1024t_1280.torch"
BATCH=25                              # 33-layer model on an 80GB GPU

# --- human ---
python eval_single_anndata.py \
  --adata_path "../human_lung.h5ad" \
  --dir "${OUTDIR}" \
  --species human \
  --model_loc "${MODEL}" \
  --nlayers 33 \
  --batch_size "${BATCH}"

# --- mouse ---  (only if you have mouse data)
python eval_single_anndata.py \
  --adata_path "../mouse_lung.h5ad" \
  --dir "${OUTDIR}" \
  --species mouse \
  --model_loc "${MODEL}" \
  --nlayers 33 \
  --batch_size "${BATCH}"

# No local GPU? Run just these two commands on a cloud GPU (one A100/H100/A10),
# then copy uce_out/*_lung.h5ad back. Stages 02 and 03 (R) are CPU-only.
