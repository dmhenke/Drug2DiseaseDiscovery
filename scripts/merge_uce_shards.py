#!/usr/bin/env python
"""Concatenate per-shard UCE outputs into one species file, in the input's cell order.

Usage: python scripts/merge_uce_shards.py IN.h5ad UCE_DIR OUT.h5ad
Expects UCE_DIR/<stem>_sNN_uce_adata.h5ad. Fails if a shard is missing or the cell order
or count differs from IN.h5ad (a recorded failure, not a silent skip).
"""
import glob
import os
import sys

import anndata as ad
import numpy as np

inp, uce_dir, out = sys.argv[1:4]
stem = os.path.basename(inp)[:-len(".h5ad")]
shards = sorted(glob.glob(os.path.join(uce_dir, f"{stem}_s[0-9][0-9]_uce_adata.h5ad")))
if not shards:
    sys.exit(f"ERROR: no shard outputs for {stem} in {uce_dir}")
ref = ad.read_h5ad(inp, backed="r").obs_names.to_numpy()
parts, names = [], []
for f in shards:
    s = ad.read_h5ad(f, backed="r")
    parts.append(np.asarray(s.obsm["X_uce"]))
    names.append(s.obs_names.to_numpy())
    s.file.close()
names = np.concatenate(names)
if len(names) != len(ref) or not np.array_equal(names, ref):
    sys.exit(f"ERROR: shard cells ({len(names)}) do not match input ({len(ref)}) in order")
full = ad.read_h5ad(inp)
full.obsm["X_uce"] = np.vstack(parts)
full.write_h5ad(out)
print(f"merged {len(shards)} shards -> {out} ({full.n_obs} cells x {full.obsm['X_uce'].shape[1]} dims)")
