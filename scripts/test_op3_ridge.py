#!/usr/bin/env python
"""Planted signal must beat the baselines; pure noise must not. Run: python scripts/test_op3_ridge.py"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from op3_ridge import Experiment, folds, summarise, FOLDS  # noqa: E402
import op3_ridge  # noqa: E402

op3_ridge.PERMUTATIONS, op3_ridge.BOOTSTRAPS = 5, 200
seed = pd.read_csv(Path(__file__).parent.parent / "shared/compounds/seed_compounds.tsv", sep="\t").dropna(subset=["smiles"]).drop_duplicates("smiles")
cells = ["T cells", "B cells"]
obs = pd.DataFrame([(n, c, s, "train") for n, s in zip(seed.name, seed.smiles) for c in cells],
                   columns=["sm_name", "cell_type", "SMILES", "file"]).drop_duplicates(["sm_name", "cell_type"])


def run(planted):
    rng = np.random.default_rng(1)
    e = Experiment(obs, np.zeros((len(obs), 300)), 0)
    beta = rng.normal(size=(e.fp.shape[1], 300)) * (0.1 if planted else 0.0)
    e.Y = e.fp[e.rc] @ beta + rng.normal(size=e.Y.shape) + np.array([[1.0], [-1.0]])[e.ct] * rng.normal(size=300)
    rows = np.arange(len(obs))
    per = pd.concat([e.evaluate("whole_compound", rows[~m], rows[m])
                     for m in (np.isin(e.groups, h) for h in folds(e.groups, FOLDS, e.rng))])
    s = summarise(per, e.rng).query("metric == 'r'").set_index("comparison").ci_lo
    return s["ridge_minus_cell_type_mean"], s["ridge_minus_permuted"]


assert all(x > 0 for x in run(True)), "planted signal not recovered"
assert not all(x > 0 for x in run(False)), "ridge beat the baselines on pure noise"
print("op3_ridge self-test PASS")
