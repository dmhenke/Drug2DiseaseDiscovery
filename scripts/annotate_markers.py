#!/usr/bin/env python
"""Label PBMC clusters from canonical markers, for datasets that ship without cell-type labels.

Leiden clusters on log-normalised counts, each assigned the type whose markers score highest. Cluster-level, so an
interferon-high monocyte is still a monocyte. Method is marker-based: record it, and do not read it as author labels.

Validate against author labels before trusting it (expectation, written before the first real run: at least 90 %
agreement for classical monocytes, NK and B cells; CD4 vs CD8 T cells will be weaker):
  python scripts/annotate_markers.py IN.h5ad OUTDIR --truth cell_type --symbols feature_name
"""
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# Names match the Perez 2022 CELLxGENE labels, so cohorts compare by name.
MARKERS = {
    "classical monocyte": ["CD14", "LYZ", "S100A8", "S100A9", "VCAN"],
    "non-classical monocyte": ["FCGR3A", "MS4A7", "CDKN1C", "LST1"],
    "natural killer cell": ["GNLY", "NKG7", "KLRD1", "KLRF1"],
    "CD4-positive, alpha-beta T cell": ["CD3D", "CD3E", "IL7R", "CD4"],
    "CD8-positive, alpha-beta T cell": ["CD3D", "CD3E", "CD8A", "CD8B"],
    "B cell": ["MS4A1", "CD79A", "CD79B"],
    "plasmablast": ["MZB1", "JCHAIN", "XBP1"],
    "conventional dendritic cell": ["FCER1A", "CD1C", "CLEC10A"],
    "plasmacytoid dendritic cell": ["LILRA4", "CLEC4C", "IRF7"],
    "platelet": ["PPBP", "PF4"],
    "erythrocyte": ["HBB", "HBA1"],
}


def annotate(a, resolution=2.0, seed=1234):
    """Return (labels, cluster table). `a` holds raw counts in X and symbols in var_names; it is not modified."""
    b = ad.AnnData(a.X.copy(), var=pd.DataFrame(index=a.var_names))
    b.var_names_make_unique()
    sc.pp.normalize_total(b, target_sum=1e4)
    sc.pp.log1p(b)
    log = b.copy()
    sc.pp.highly_variable_genes(b, n_top_genes=2000, subset=True)
    sc.pp.scale(b, max_value=10)
    sc.tl.pca(b, n_comps=30, random_state=seed)
    sc.pp.neighbors(b, n_neighbors=15, random_state=seed)
    sc.tl.leiden(b, resolution=resolution, random_state=seed, flavor="igraph", n_iterations=2, directed=False)
    clusters = b.obs.leiden.astype(str).to_numpy()
    ids = sorted(set(clusters), key=int)
    genes = sorted({g for m in MARKERS.values() for g in m if g in log.var_names})
    mean = pd.DataFrame({c: np.asarray(log[clusters == c][:, genes].X.mean(0)).ravel() for c in ids}, index=genes).T
    z = (mean - mean.mean()) / mean.std().replace(0, 1)
    score = pd.DataFrame({t: z[[g for g in m if g in genes]].mean(1) for t, m in MARKERS.items()})
    score = score.fillna(-99)  # a type with no measured marker can never win
    top = score.idxmax(1)
    ordered = np.sort(score.to_numpy(), axis=1)
    table = pd.DataFrame({"cluster": ids, "cells": [int((clusters == c).sum()) for c in ids], "label": top.to_numpy(),
                          "score": ordered[:, -1], "margin_to_second": ordered[:, -1] - ordered[:, -2]})
    return pd.Series(clusters).map(top).to_numpy(), table


def main():
    p = argparse.ArgumentParser()
    p.add_argument("h5ad"); p.add_argument("out", type=Path)
    p.add_argument("--truth", required=True); p.add_argument("--symbols")
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=True)
    x = ad.read_h5ad(a.h5ad)
    if a.symbols:
        x.var_names = x.var[a.symbols].astype(str).to_numpy()
    label, table = annotate(x)
    truth = x.obs[a.truth].astype(str).to_numpy()
    tab = pd.crosstab(pd.Series(truth, name="author"), pd.Series(label, name="marker"))
    tab.to_csv(a.out / "annotation_vs_author.tsv", sep="\t")
    table.to_csv(a.out / "annotation_clusters.tsv", sep="\t", index=False, float_format="%.3f")
    agree = pd.Series({t: (label[truth == t] == t).mean() for t in MARKERS if (truth == t).any()}).round(3)
    agree.rename("agreement").to_csv(a.out / "annotation_agreement.tsv", sep="\t", index_label="author_label")
    print(f"{len(truth)} cells; agreement with author labels (want >= 0.90 for classical monocyte, NK, B):")
    print(agree.to_string())
    print("ANNOTATION CHECK:", "PASS" if min(agree.get(k, 0) for k in
          ("classical monocyte", "natural killer cell", "B cell")) >= 0.9 else "FAIL")


if __name__ == "__main__":
    main()
