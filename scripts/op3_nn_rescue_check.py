#!/usr/bin/env python
"""If a drug had not been measured, would a predicted profile have ranked it correctly as a lupus rescuer?

Holds each OP3 compound out, predicts its profile on the SLE signature genes from the structurally nearest
other compound (same cell type), scores it like a measured one, and compares with its measured score.
  rescue = -(t - m).d / |d|^2   m = median compound; 1 = the specific part of t cancels the disease exactly
Also: same compound measured in another cell type, which approximates the CMap route (measured elsewhere).

Expectation, written before the run: Stage P r_resid 0.095 is weak, so rank agreement should be low overall,
better above Tanimoto 0.4, and JAK inhibitors should not be recovered from structure (ruxolitinib's nearest is
crizotinib at 0.18). Cross-cell-type agreement should be clearly higher; if not, CMap will not rescue this.
Run on signature genes only, so this needs no OP3 matrix. Intervals resample compounds; pairs are not independent.
Usage: python scripts/op3_nn_rescue_check.py outdir [--seed 0]
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from mol_fingerprints import fingerprints  # noqa: E402

BOOTS, PERMS, TOP = 2000, 2000, 10
BINS = [(0, .25), (.25, .4), (.4, 1.01)]


def rescue(t, m, d):
    return -((t - m) @ d) / (d @ d)


def rho(a, b):
    return spearmanr(a, b).statistic


def interval(x, y, rng):
    n = len(x)
    v = [rho(x[i], y[i]) for i in (rng.integers(0, n, n) for _ in range(BOOTS))]
    return np.nanpercentile(v, [2.5, 97.5])


def self_check():
    d = np.array([2.0, -1.0])
    assert np.isclose(rescue(-d, 0 * d, d), 1) and np.isclose(rescue(d, d, d), 0)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("out", type=Path); p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    self_check()
    rng = np.random.default_rng(a.seed)
    root = Path(__file__).parent.parent / "analysis/op3_rescue"
    genes = pd.read_csv(root / "rescue_genes.tsv", sep="\t")
    smiles = pd.read_csv(root / "rescue.tsv", sep="\t").drop_duplicates("sm_name").set_index("sm_name").SMILES
    fp, _, bad = fingerprints(smiles.to_list())
    assert not bad
    inter = fp @ fp.T
    n = fp.sum(1)
    tan = pd.DataFrame(inter / (n[:, None] + n[None, :] - inter), smiles.index, smiles.index)
    np.fill_diagonal(tan.values, -1)

    score, rows, per = {}, [], []
    for pop, g in genes.groupby("population", sort=False):
        T = g.pivot(index="sm_name", columns="gene", values="drug")
        assert not T.isna().any().any()
        d = g.drop_duplicates("gene").set_index("gene").disease[T.columns].to_numpy()
        m = np.median(T.to_numpy(), 0)
        meas = pd.Series([rescue(t, m, d) for t in T.to_numpy()], T.index)
        score[pop] = meas
        S = tan.loc[T.index, T.index]
        near = S.idxmax(1)
        pred = pd.Series([rescue(T.loc[near[c]].to_numpy(), m, d) for c in T.index], T.index)
        x, y, sim = pred.to_numpy(), meas.to_numpy(), S.max(1).to_numpy()
        top = set(meas.nlargest(TOP).index)
        hits = len(top & set(pred.nlargest(TOP).index))
        null = [rho(rng.permutation(x), y) for _ in range(PERMS)]
        r = rho(x, y)
        rows.append(dict(population=pop, route="structure_nearest", subset="all", n=len(x), spearman=r,
                         ci_lo=(ci := interval(x, y, rng))[0], ci_hi=ci[1], p_perm=(1 + np.sum(np.abs(null) >= abs(r))) / (PERMS + 1),
                         top10_recovered=hits, top10_expected_by_chance=TOP ** 2 / len(x)))
        for lo, hi in BINS:
            k = (sim >= lo) & (sim < hi)
            if k.sum() >= 5:
                ci = interval(x[k], y[k], rng)
                rows.append(dict(population=pop, route="structure_nearest", subset=f"tanimoto {lo}-{hi}", n=int(k.sum()),
                                 spearman=rho(x[k], y[k]), ci_lo=ci[0], ci_hi=ci[1]))
        rk = meas.rank(ascending=False)
        rp = pred.rank(ascending=False)
        per.append(pd.DataFrame(dict(population=pop, measured=meas, predicted=pred, nearest=near, tanimoto=sim,
                                     measured_rank=rk, predicted_rank=rp)))
    # Same compound, other cell type. Non-classical shares Myeloid cells with classical, so it is not a transfer.
    other = {"Classical mono": ["NK", "CD4 T"], "NK": ["Classical mono", "CD4 T"], "CD4 T": ["Classical mono", "NK"]}
    for pop, srcs in other.items():
        for s in srcs:
            both = score[pop].index.intersection(score[s].index)
            x, y = score[s][both].to_numpy(), score[pop][both].to_numpy()
            ci = interval(x, y, rng)
            rows.append(dict(population=pop, route=f"measured_in_{s}", subset="all", n=len(x), spearman=rho(x, y), ci_lo=ci[0], ci_hi=ci[1]))
    res, per = pd.DataFrame(rows), pd.concat(per).reset_index()
    a.out.mkdir(parents=True, exist_ok=True)
    res.to_csv(a.out / "nn_rescue_summary.tsv", sep="\t", index=False, float_format="%.3f")
    per.to_csv(a.out / "nn_rescue_per_compound.tsv", sep="\t", index=False, float_format="%.3f")
    print(res.round(3).to_string(index=False))
    print(f"\nTests: {res.p_perm.notna().sum()} permutation; bin and cross-type intervals are descriptive, unadjusted.")
    ctrl = per[per.sm_name.isin(["Ruxolitinib", "GLPG0634", "Prednisolone", "Mometasone Furoate", "Methotrexate"]) & (per.population == "Classical mono")]
    print(ctrl.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
