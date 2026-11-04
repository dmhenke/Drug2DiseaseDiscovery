#!/usr/bin/env python
"""Stage P of planning/op3_pilot.md: can a compound's structure predict its effect on PBMCs?

Cell-type mean + ridge on Morgan fingerprint -> per-gene DE, against four baselines. The cell-type term is
not penalised, so a large alpha falls back to the cell-type mean rather than below it. Additive, so it
cannot learn a compound effect that differs by cell type; that is deliberate for a baseline.

Two evaluations:
  whole_compound  PRIMARY. 5-fold, every row of a held-out compound removed. Unseen chemistry.
  official        OP3's split. All 76 test compounds have T and NK rows in train, so this measures
                  transfer between cell types, not unseen chemistry. Adds a same_compound baseline.

Expectation, written before any run: compound-specific signal is weak and chemistry is sparse (median
nearest-neighbour Tanimoto 0.26), so a narrow win or a loss to nearest_compound would not be surprising.
A large win on unseen compounds would be, and should be checked for leakage before it is believed.
Usage: python scripts/op3_ridge.py shared/perturbation/op3 outdir [--layer logFC] [--seed 0]
"""
import argparse
import socket
import subprocess
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.linear_model import Ridge

sys.path.insert(0, str(Path(__file__).parent))
from mol_fingerprints import fingerprints  # noqa: E402

ALPHAS = 10.0 ** np.arange(0, 7)
FOLDS, PERMUTATIONS, BOOTSTRAPS = 5, 20, 2000
TARGETS = ("TSC22D3", "FKBP5")


def load(d, layer):
    tr, te = ad.read_h5ad(d / "de_train.h5ad"), ad.read_h5ad(d / "de_test.h5ad")
    assert list(tr.var_names) == list(te.var_names), "train and test genes differ"
    dense = lambda m: m.toarray() if hasattr(m, "toarray") else np.asarray(m)  # noqa: E731
    Y = np.vstack([dense(tr.layers[layer]), dense(te.layers[layer])]).astype(np.float64)
    n_nan = int(np.isnan(Y).sum())
    Y = np.nan_to_num(Y)
    obs = pd.concat([tr.obs.assign(file="train"), te.obs.assign(file="test")], ignore_index=True)
    obs = obs[["sm_name", "cell_type", "SMILES", "file"]].astype(str)
    return obs, Y, list(tr.var_names), n_nan


def design(obs):
    """Compound fingerprints, Tanimoto matrix, row->compound index and leakage groups."""
    comp = obs.drop_duplicates("sm_name").reset_index(drop=True)
    fp, keys, bad = fingerprints(comp.SMILES)
    assert not bad, f"unparseable SMILES: {comp.sm_name[bad].tolist()}"
    assert comp.SMILES.nunique() == len(comp), "one structure under two names"
    inter = fp @ fp.T
    n = fp.sum(1)
    tan = inter / (n[:, None] + n[None, :] - inter)
    row_comp = obs.sm_name.map({s: i for i, s in enumerate(comp.sm_name)}).to_numpy()
    # InChIKey connectivity block, so stereoisomers of one molecule can never straddle a split
    block = pd.factorize([k[:14] for k in keys])[0]
    return comp, fp, tan, row_comp, block[row_comp]


def rowcorr(a, b):
    a, b = a - a.mean(1, keepdims=True), b - b.mean(1, keepdims=True)
    den = np.sqrt((a * a).sum(1) * (b * b).sum(1))
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(den > 0, (a * b).sum(1) / den, np.nan)


def folds(groups, k, rng):
    return np.array_split(rng.permutation(np.unique(groups)), k)


class Experiment:
    def __init__(self, obs, Y, seed):
        self.obs, self.Y, self.rng = obs, Y, np.random.default_rng(seed)
        self.comp, self.fp, self.tan, self.rc, self.groups = design(obs)
        self.cells = sorted(obs.cell_type.unique())
        self.ct = obs.cell_type.map({c: i for i, c in enumerate(self.cells)}).to_numpy()

    def ridge(self, trn, tst, alpha, fp=None):
        fp = self.fp if fp is None else fp
        _, ctm = self.ct_mean(trn)
        fit = Ridge(alpha=alpha).fit(fp[self.rc[trn]], self.Y[trn] - ctm[self.ct[trn]])
        return ctm[self.ct[tst]] + fit.predict(fp[self.rc[tst]])

    def select_alpha(self, trn):
        """Inner compound-grouped CV on training rows only, so no held-out compound picks alpha."""
        scores = np.zeros(len(ALPHAS))
        for held in folds(self.groups[trn], 4, self.rng):
            mask = np.isin(self.groups[trn], held)
            a, b = trn[~mask], trn[mask]
            for i, alpha in enumerate(ALPHAS):
                scores[i] += np.nanmean(rowcorr(self.ridge(a, b, alpha), self.Y[b]))
        return float(ALPHAS[scores.argmax()])

    def ct_mean(self, trn):
        glob = self.Y[trn].mean(0)
        return glob, np.vstack([self.Y[trn][self.ct[trn] == i].mean(0) if (self.ct[trn] == i).any() else glob
                                for i in range(len(self.cells))])

    def nearest(self, trn, tst):
        """The most similar training compound that was measured in the same cell type."""
        pred, who, sim = np.zeros((len(tst), self.Y.shape[1])), [], []
        for k, r in enumerate(tst):
            cand = trn[self.ct[trn] == self.ct[r]]
            sims = self.tan[self.rc[r], self.rc[cand]]
            best = cand[sims.argmax()]
            pred[k], sim_k = self.Y[best], sims.max()
            who.append(self.comp.sm_name[self.rc[best]])
            sim.append(float(sim_k))
        return pred, who, sim

    def evaluate(self, name, trn, tst):
        assert not (name == "whole_compound" and set(self.groups[trn]) & set(self.groups[tst])), "compound leakage"
        alpha = self.select_alpha(trn)
        glob, ctm = self.ct_mean(trn)
        truth, base = self.Y[tst], ctm[self.ct[tst]]
        preds = {"global_mean": np.tile(glob, (len(tst), 1)), "cell_type_mean": base,
                 "ridge": self.ridge(trn, tst, alpha)}
        preds["nearest_compound"], who, sim = self.nearest(trn, tst)
        if name == "official":
            preds["same_compound"] = np.vstack([self.Y[trn][self.rc[trn] == self.rc[r]].mean(0) for r in tst])
        out = pd.DataFrame({"evaluation": name, "sm_name": self.obs.sm_name.to_numpy()[tst],
                            "cell_type": self.obs.cell_type.to_numpy()[tst],
                            "nearest": who, "tanimoto": sim, "alpha": alpha})
        for m, p in preds.items():
            out[f"r.{m}"] = rowcorr(p, truth)
            # Correlation after removing the cell-type mean: is any compound-specific signal captured?
            out[f"r_resid.{m}"] = rowcorr(p - base, truth - base)
            out[f"rmse.{m}"] = np.sqrt(((p - truth) ** 2).mean(1))
        # Null: break the structure -> response link among training compounds, keep everything else.
        perm = {k: [] for k in ("r", "r_resid", "rmse")}
        ids = np.unique(self.rc[trn])
        for _ in range(PERMUTATIONS):
            fp = self.fp.copy()
            fp[ids] = self.fp[self.rng.permutation(ids)]
            p = self.ridge(trn, tst, alpha, fp)
            perm["r"].append(rowcorr(p, truth))
            perm["r_resid"].append(rowcorr(p - base, truth - base))
            perm["rmse"].append(np.sqrt(((p - truth) ** 2).mean(1)))
        for k, v in perm.items():
            out[f"{k}.permuted"] = np.mean(v, axis=0)
        return out

    def biology(self, genes):
        """Prednisolone x myeloid with prednisolone wholly unseen; then with no glucocorticoid at all."""
        names = self.obs.sm_name.str.lower().to_numpy()
        pred_rows = names == "prednisolone"
        tst = np.flatnonzero(pred_rows & (self.obs.cell_type.str.contains("Myeloid").to_numpy()))
        rows = []
        for label, drop in (("mometasone_in_training", pred_rows), ("no_glucocorticoid", pred_rows | (names == "mometasone furoate"))):
            trn = np.flatnonzero(~drop)
            ridge = self.ridge(trn, tst, self.select_alpha(trn))[0]
            _, ctm = self.ct_mean(trn)
            near, who, sim = self.nearest(trn, tst)
            for g in TARGETS:
                j = genes.index(g)
                rows.append({"training": label, "gene": g, "measured": self.Y[tst[0], j], "ridge": ridge[j],
                             "cell_type_mean": ctm[self.ct[tst[0]], j], "nearest_compound": near[0, j],
                             "nearest_name": who[0], "tanimoto": sim[0]})
        return pd.DataFrame(rows)


def boot(x, rng):
    x = x[~np.isnan(x)]
    means = rng.choice(x, (BOOTSTRAPS, len(x))).mean(1)
    return x.mean(), *np.percentile(means, [2.5, 97.5]), len(x)


def summarise(per_row, rng):
    """Every interval resamples compounds, not rows: rows of one compound are not independent."""
    rows = []
    for ev, df in per_row.groupby("evaluation"):
        byc = df.groupby("sm_name").mean(numeric_only=True)
        for metric in ("r", "r_resid", "rmse"):
            cols = [c for c in byc if c.startswith(metric + ".")]
            for c in cols:
                m = c.split(".", 1)[1]
                est, lo, hi, n = (boot(byc[c].to_numpy(), rng) if byc[c].notna().any() else (np.nan,) * 4)
                rows.append(dict(evaluation=ev, metric=metric, comparison=m, estimate=est, ci_lo=lo, ci_hi=hi, n_compounds=n))
                if m != "ridge" and byc[c].notna().any():
                    est, lo, hi, n = boot((byc[f"{metric}.ridge"] - byc[c]).to_numpy(), rng)
                    rows.append(dict(evaluation=ev, metric=metric, comparison=f"ridge_minus_{m}", estimate=est, ci_lo=lo, ci_hi=hi, n_compounds=n))
    return pd.DataFrame(rows)


def verdict(summary, bio):
    """Criteria from planning/op3_pilot.md, fixed before the first run. Beats = paired CI above zero."""
    s = summary[(summary.evaluation == "whole_compound") & (summary.metric == "r")].set_index("comparison")
    won = lambda m: bool(s.loc[f"ridge_minus_{m}", "ci_lo"] > 0)  # noqa: E731
    b = bio[bio.training == "mometasone_in_training"]
    return {"1 beats global mean, cell-type mean and permuted": won("global_mean") and won("cell_type_mean") and won("permuted"),
            "2 beats nearest training compound": won("nearest_compound"),
            "3 unseen prednisolone: TSC22D3 and FKBP5 predicted up, above the cell-type mean":
                bool(((b.ridge > 0) & (b.ridge > b.cell_type_mean)).all())}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("data", type=Path); p.add_argument("out", type=Path)
    p.add_argument("--layer", default="logFC"); p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    obs, Y, genes, n_nan = load(a.data, a.layer)
    e = Experiment(obs, Y, a.seed)
    root = Path(__file__).parent.parent
    git = lambda *c: subprocess.run(["git", "-C", str(root), *c], capture_output=True, text=True).stdout.strip()  # noqa: E731
    head = [f"STAGE P run on {a.data}: {len(obs)} rows, {len(e.comp)} compounds, {len(genes)} genes, layer {a.layer}",
            f"commit {git('rev-parse', '--short', 'HEAD')}; {len(git('status', '--porcelain').splitlines())} changed files; "
            f"host {socket.gethostname()}; seed {a.seed}; NaN values set to 0: {n_nan}",
            f"folds {FOLDS}, permutations {PERMUTATIONS}, bootstraps {BOOTSTRAPS}, alphas {ALPHAS.tolist()}"]
    per = []
    rows = np.arange(len(obs))
    for held in folds(e.groups, FOLDS, e.rng):
        mask = np.isin(e.groups, held)
        per.append(e.evaluate("whole_compound", rows[~mask], rows[mask]))
    per.append(e.evaluate("official", rows[obs.file == "train"], rows[obs.file == "test"]))
    per = pd.concat(per, ignore_index=True)
    summary = summarise(per, e.rng)
    bio = e.biology(genes)
    per.to_csv(a.out / "per_row.tsv", sep="\t", index=False, float_format="%.4f")
    summary.to_csv(a.out / "summary.tsv", sep="\t", index=False, float_format="%.4f")
    bio.to_csv(a.out / "prednisolone.tsv", sep="\t", index=False, float_format="%.3f")
    crit = verdict(summary, bio)
    show = summary[summary.metric.isin(["r", "r_resid"])]
    lines = head + ["", show.round(3).to_string(index=False), "", bio.round(2).to_string(index=False), "",
                    f"alphas chosen: {sorted(per.alpha.unique().tolist())}", ""]
    lines += [f"  [{'ok' if v else 'FAIL'}] {k}" for k, v in crit.items()]
    lines += ["Criteria use whole_compound, metric r; CIs over compounds, unadjusted, all 5 comparisons must clear zero.",
              f"STAGE P ({a.layer}): " + ("PASS" if all(crit.values()) else "FAIL -> go to the pilot's fail branch")]
    (a.out / "verdict.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
