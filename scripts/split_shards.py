#!/usr/bin/env python
"""Split an interim species h5ad into contiguous cell shards (deterministic order).

UCE takes one species per call and holds a dense copy of each input in RAM, so large
species files are embedded shard by shard. Shards keep all genes and obs.

Usage: python scripts/split_shards.py IN.h5ad OUTDIR CELLS_PER_SHARD
Writes OUTDIR/<stem>_s00.h5ad, _s01.h5ad, ...
"""
import os
import sys

import anndata as ad

inp, outdir, per = sys.argv[1], sys.argv[2], int(sys.argv[3])
stem = os.path.basename(inp)[:-len(".h5ad")]
src = ad.read_h5ad(inp, backed="r")
n = src.n_obs
starts = list(range(0, n, per))
for i, s in enumerate(starts):
    src[s:min(s + per, n)].to_memory().write_h5ad(os.path.join(outdir, f"{stem}_s{i:02d}.h5ad"))
print(f"split {n} cells into {len(starts)} shards of <= {per}")
