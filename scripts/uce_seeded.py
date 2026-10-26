#!/usr/bin/env python
"""Run UCE's eval_single_anndata.py with fixed seeds. Run it from the UCE repo directory.

UCE samples genes with np.random and sets no seed upstream, so repeated runs differ
(report finding C2). Arguments are passed through unchanged.
Env: UCE_SEED (default 1234)

Usage: python /path/to/scripts/uce_seeded.py --adata_path ... --dir ... --species human ...
"""
import os
import random
import runpy
import sys

import numpy as np
import torch

seed = int(os.environ.get("UCE_SEED", "1234"))
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
print(f"uce_seeded: seed={seed}", flush=True)
sys.path.insert(0, os.getcwd())
sys.argv = ["eval_single_anndata.py"] + sys.argv[1:]
runpy.run_path("eval_single_anndata.py", run_name="__main__")
