#!/usr/bin/env python
"""PCA stand-in for UCE. NOT UCE: baseline and plumbing test only.

Writes obsm["X_uce"] = PCA of log-normalised counts; X stays raw counts (UCE's contract).
Lets Stages 0, 2 and 3 run with no model, GPU, or download.

Usage: python scripts/fake_uce.py IN.h5ad OUT.h5ad [n_comps=50]
"""
import sys

import anndata as ad
import scanpy as sc

inp, out = sys.argv[1], sys.argv[2]
n = int(sys.argv[3]) if len(sys.argv) > 3 else 50
adata = ad.read_h5ad(inp)
adata.layers["counts"] = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
n = min(n, adata.n_obs - 1, adata.n_vars - 1)
sc.tl.pca(adata, n_comps=n, svd_solver="arpack", random_state=1234)
adata.obsm["X_uce"] = adata.obsm.pop("X_pca")
adata.X = adata.layers.pop("counts")
adata.write_h5ad(out)
print(f"fake_uce (PCA, NOT UCE): {adata.n_obs} cells -> X_uce {adata.obsm['X_uce'].shape}")
