#!/usr/bin/env python
"""How far does each measured OP3 compound move an SLE cell population toward its healthy state?

Per population, on the signature genes OP3 measured: d = log2FC SLE vs healthy, t = log2FC compound vs DMSO.
  match_gain = 1 - ||d + t|| / ||d||    1 = exactly healthy, 0 = unchanged, < 0 = further from healthy
This ASSUMES disease and drug effects add in log space. OP3 measured healthy cells, so nothing here observes
a treated SLE cell. Empirical p: same d on random non-signature genes. BH within each population.

Expectation, written before the run: interferon-dominated monocyte signatures reward compounds that lower ISGs.
Cultured healthy PBMCs have low baseline interferon, so JAK inhibitors may show little (Inferred). Broad
transcriptional poisons could score well by lowering everything: read match_gain with off_target_rms.
Usage: python scripts/op3_rescue.py shared/perturbation/op3 analysis/sle_signatures/signatures.tsv outdir [--seed 0]
"""
import argparse
import socket
import subprocess
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

NULLS, BOOTS = 20000, 1000  # NULLS must allow q < 0.05 for one hit among 140: p floor x 140 < 0.05
CONTROLS = ["Ruxolitinib", "GLPG0634", "Prednisolone", "Mometasone Furoate", "Methotrexate"]


def gain(d, t):
    """d is 1-D (genes); t is (..., genes)."""
    return 1 - np.linalg.norm(d + t, axis=-1) / np.linalg.norm(d)


def bh(p):
    p = np.asarray(p)
    order = np.argsort(p)
    q = p[order] * len(p) / np.arange(1, len(p) + 1)
    out = np.empty(len(p))
    out[order] = np.minimum.accumulate(q[::-1])[::-1].clip(max=1)
    return out


def self_check():
    d = np.array([2.0, -1.0, 0.5])
    assert np.isclose(gain(d, -d), 1) and np.isclose(gain(d, 0 * d), 0) and np.isclose(gain(d, d), -1)
    assert np.allclose(bh([0.01, 0.04, 0.03]), [0.03, 0.04, 0.04])


def load(data):
    tr, te = ad.read_h5ad(data / "de_train.h5ad"), ad.read_h5ad(data / "de_test.h5ad")
    assert list(tr.var_names) == list(te.var_names), "train and test genes differ"
    dense = lambda m: m.toarray() if hasattr(m, "toarray") else np.asarray(m)  # noqa: E731
    Y = np.vstack([dense(tr.layers["logFC"]), dense(te.layers["logFC"])]).astype(float)
    obs = pd.concat([tr.obs, te.obs], ignore_index=True)[["sm_name", "cell_type", "SMILES"]].astype(str)
    assert not np.isnan(Y).any(), "NaN in logFC"
    return obs, Y, list(tr.var_names)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("data", type=Path); p.add_argument("signatures", type=Path); p.add_argument("out", type=Path)
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    self_check()
    a.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    obs, Y, genes = load(a.data)
    sig = pd.read_csv(a.signatures, sep="\t")
    col = {g: i for i, g in enumerate(genes)}
    # ISGs are the answer key of the lupus run: if these are not up, the signature's sign is wrong.
    mono = sig[sig.population == "Classical mono"].set_index("gene").log2fc_sle_vs_healthy
    assert mono.get("IFI44L", 0) > 0, "IFI44L not up in SLE classical monocytes: signature sign is wrong"

    rows, long, cover = [], [], []
    for pop, s in sig.groupby("population", sort=False):
        s = s[s.gene.isin(col)]
        n_all = int((sig.population == pop).sum())
        idx = s.gene.map(col).to_numpy()
        d = s.log2fc_sle_vs_healthy.to_numpy()
        cell = s.op3_cell_type.iat[0]
        keep = np.flatnonzero(obs.cell_type == cell)
        T = Y[keep]
        others = np.setdiff1d(np.arange(Y.shape[1]), idx)
        score = gain(d, T[:, idx])
        null = np.stack([gain(d, T[:, rng.choice(others, len(idx), replace=False)]) for _ in range(NULLS)], 1)
        pval = (1 + (null >= score[:, None]).sum(1)) / (NULLS + 1)
        boot = rng.integers(0, len(idx), (BOOTS, len(idx)))
        ci = np.stack([gain(d[b], T[:, idx[b]]) for b in boot], 1)
        cover.append(dict(population=pop, op3_cell_type=cell, signature_genes=n_all, measured_in_op3=len(idx),
                          compounds=len(keep), tests_for_bh=len(keep)))
        moved = np.abs(d + T[:, idx]) < np.abs(d)
        q = bh(pval)
        for i, r in enumerate(keep):
            t = T[i, idx]
            rows.append(dict(population=pop, op3_cell_type=cell, sm_name=obs.sm_name[r], SMILES=obs.SMILES[r],
                             evidence="measured", match_gain=score[i], ci_lo=np.percentile(ci[i], 2.5),
                             ci_hi=np.percentile(ci[i], 97.5), p_perm=pval[i], q_bh=q[i],
                             genes_toward_healthy=int(moved[i].sum()), genes_scored=len(idx),
                             reversal_cosine=-float(d @ t / (np.linalg.norm(d) * np.linalg.norm(t))),
                             off_target_rms=float(np.sqrt((T[i, others] ** 2).mean()))))
            long += [dict(population=pop, sm_name=obs.sm_name[r], gene=g, disease=x, drug=y) for g, x, y in zip(s.gene, d, t)]
    res = pd.DataFrame(rows).sort_values(["population", "match_gain"], ascending=[True, False])
    res.to_csv(a.out / "rescue.tsv", sep="\t", index=False, float_format="%.4f")
    pd.DataFrame(long).to_csv(a.out / "rescue_genes.tsv", sep="\t", index=False, float_format="%.3f")
    cov = pd.DataFrame(cover)
    cov.to_csv(a.out / "coverage.tsv", sep="\t", index=False)

    root = Path(__file__).parent.parent
    git = lambda *c: subprocess.run(["git", "-C", str(root), *c], capture_output=True, text=True).stdout.strip()  # noqa: E731
    lines = [f"Run on {a.data}; commit {git('rev-parse', '--short', 'HEAD')}; "
             f"{len(git('status', '--porcelain').splitlines())} changed files; host {socket.gethostname()}; seed {a.seed}; "
             f"nulls {NULLS}; gene bootstraps {BOOTS}", "", cov.to_string(index=False), ""]
    for c, n in zip(cov.population, cov.measured_in_op3 / cov.signature_genes):
        lines.append(f"  [{'ok' if n >= 0.5 else 'FAIL'}] {c}: {n:.0%} of signature genes measured (want at least 50%)")
    lines += ["", "Population, compounds with q < 0.05 and match_gain > 0:"]
    for pop, g in res.groupby("population", sort=False):
        lines.append(f"  {pop}: {int(((g.q_bh < .05) & (g.match_gain > 0)).sum())} of {len(g)}; median match_gain {g.match_gain.median():.3f}")
    lines += ["", "Positive controls, rank by match_gain (recorded, not a gate):"]
    for pop, g in res.groupby("population", sort=False):
        g = g.reset_index(drop=True)
        for c in CONTROLS:
            for i in g.index[g.sm_name == c]:
                lines.append(f"  {pop}: {c} rank {i + 1}/{len(g)}, match_gain {g.match_gain[i]:.3f}, q {g.q_bh[i]:.3f}")
    (a.out / "run_summary.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()
