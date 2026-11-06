#!/usr/bin/env python
"""Does the SLE signature replicate in independent cohorts? Writes replication tables and a consensus signature.

Discovery: padj < 0.05 in the discovery study, ribosomal and mitochondrial genes removed. A gene REPLICATES when a
validation cohort has the same sign and p < 0.05. Background: sign agreement of genes with discovery padj > 0.5,
which shows how much agreement a shared artefact alone produces. Consensus log2FC is the inverse-variance mean of
the cohorts that tested the gene; it is written in the format scripts/op3_rescue.py reads.

Expectation, written before the first real run: monocyte interferon genes (IFI27, IFI44L, ISG15) replicate in both
GSE135779 cohorts, and concordance beats background in monocytes. CD4 T ribosomal hits are excluded by the filter,
so they cannot replicate here by construction.
Usage: python scripts/compare_cohorts.py outdir discovery_study validation_study [validation_study ...]
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).parent.parent
POPULATIONS = {"classical monocyte": ("Classical mono", "Myeloid cells"), "non-classical monocyte": ("Non-classical mono", "Myeloid cells"),
               "natural killer cell": ("NK", "NK cells"), "CD4-positive, alpha-beta T cell": ("CD4 T", "T cells"),
               "CD8-positive, alpha-beta T cell": ("CD8 T", "T cells"), "B cell": ("B", "B cells")}
ANSWER_KEY = ("IFI27", "IFI44L", "ISG15")
RIBO = re.compile(r"^(RPL|RPS|MT-)")


def read(study):
    path = ROOT / "studies" / study / "results" / "DE_all_human.csv"
    if not path.exists():
        sys.exit(f"ERROR: {path} is missing; did the {study} run finish?")
    d = pd.read_csv(path)
    return d[d.cell_type.isin(POPULATIONS)].dropna(subset=["log2FoldChange", "lfcSE"])


def wilson(k, n, z=1.96):
    if n == 0:
        return np.nan, np.nan
    p, den = k / n, 1 + z * z / n
    mid, half = (p + z * z / (2 * n)) / den, z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return mid - half, mid + half


def compare(disc, val, cell, name):
    d, v = disc[disc.cell_type == cell].set_index("gene"), val[val.cell_type == cell].set_index("gene")
    v = v[~v.index.duplicated()]
    sig = d[(d.padj < 0.05) & ~d.index.str.match(RIBO)]
    both = sig.index.intersection(v.index)
    same = np.sign(sig.log2FoldChange[both]) == np.sign(v.log2FoldChange[both])
    bg = d[(d.padj > 0.5) & ~d.index.str.match(RIBO)].index.intersection(v.index)
    bg_same = (np.sign(d.log2FoldChange[bg]) == np.sign(v.log2FoldChange[bg])).mean() if len(bg) else np.nan
    lo, hi = wilson(int(same.sum()), len(both))
    return dict(population=POPULATIONS[cell][0], validation=name, discovery_genes=len(sig), tested=len(both),
                same_sign=int(same.sum()), concordance=same.mean() if len(both) else np.nan, ci_lo=lo, ci_hi=hi,
                background_concordance=bg_same,
                replicated_p05=int((same & (v.pvalue[both] < 0.05)).sum()),
                spearman_log2fc=stats.spearmanr(sig.log2FoldChange[both], v.log2FoldChange[both])[0] if len(both) > 2 else np.nan)


def consensus(disc, vals, cell):
    d = disc[disc.cell_type == cell].set_index("gene")
    d = d[~d.index.duplicated()]
    sig = d[(d.padj < 0.05) & ~d.index.str.match(RIBO)]
    rows = []
    for g, r in sig.iterrows():
        lfc, se, hit = [r.log2FoldChange], [r.lfcSE], []
        for name, v in vals.items():
            w = v[(v.cell_type == cell) & (v.gene == g)]
            if len(w):
                w = w.iloc[0]
                lfc.append(w.log2FoldChange); se.append(w.lfcSE)
                if np.sign(w.log2FoldChange) == np.sign(r.log2FoldChange) and w.pvalue < 0.05:
                    hit.append(name)
        if hit:
            weight = 1 / np.square(se)
            rows.append(dict(population=POPULATIONS[cell][0], op3_cell_type=POPULATIONS[cell][1], gene=g,
                             log2fc_sle_vs_healthy=round(float((np.array(lfc) * weight).sum() / weight.sum()), 4),
                             padj=f"{r.padj:.3g}", base_mean=round(r.baseMean, 1), replicated_in=",".join(hit)))
    return rows


def main():
    out, disc_name, val_names = Path(sys.argv[1]), sys.argv[2], sys.argv[3:]
    if not val_names:
        sys.exit(__doc__)
    out.mkdir(parents=True, exist_ok=True)
    disc = read(disc_name)
    vals = {n: read(n) for n in val_names}
    table = pd.DataFrame([compare(disc, v, c, n) for n, v in vals.items() for c in POPULATIONS if (disc.cell_type == c).any()])
    table.to_csv(out / "replication.tsv", sep="\t", index=False, float_format="%.3f")
    cons = pd.DataFrame([r for c in POPULATIONS for r in consensus(disc, vals, c)])
    cons.to_csv(out / "consensus_signatures.tsv", sep="\t", index=False)
    print(f"Discovery {disc_name}; {len(table)} comparisons, intervals unadjusted\n")
    print(table.round(3).to_string(index=False))
    print("\nConsensus genes per population:", cons.groupby("population").size().to_dict() if len(cons) else "NONE")
    for name, v in vals.items():
        m = v[(v.cell_type == "classical monocyte") & v.gene.isin(ANSWER_KEY)].set_index("gene")
        got = {g: round(float(m.log2FoldChange[g]), 2) for g in ANSWER_KEY if g in m.index}
        ok = len(got) == len(ANSWER_KEY) and all(x > 0 for x in got.values())
        print(f"  [{'ok' if ok else 'FAIL'}] {name}: interferon genes up in classical monocytes: {got}")


if __name__ == "__main__":
    main()
