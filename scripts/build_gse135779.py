#!/usr/bin/env python
"""Turn the GSE135779 per-sample matrices into one labelled h5ad per age cohort, with QC tables.

Children and adults stay separate: age is fully tied to cohort, and each is a replication test of its own.
Sex and treatment are not in the GEO metadata (they are in the paper's supplement), so neither is adjusted for.
Expectations from the GEO summary: ~276k cells in children, ~132k in adults; integer counts; both groups in each cohort.
Usage: python scripts/build_gse135779.py shared/gse135779 studies [--samples scripts/gse135779_samples.tsv]
"""
import argparse
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
from annotate_markers import annotate  # noqa: E402

COHORTS = {"Children": "gse135779_child", "Adult": "gse135779_adult"}


def read_sample(raw, s, symbols):
    mtx, bc = raw / f"{s.gsm}_{s.library}_matrix.mtx.gz", raw / f"{s.gsm}_{s.library}_barcodes.tsv.gz"
    m = sp.csr_matrix(scipy.io.mmread(mtx)).T.tocsr()  # 10x stores genes x cells
    barcodes = pd.read_csv(bc, header=None, sep="\t")[0].astype(str).to_numpy()
    if m.shape != (len(barcodes), len(symbols)):
        sys.exit(f"ERROR {s.donor}: matrix {m.shape}, expected {len(barcodes)} cells x {len(symbols)} genes")
    if not np.allclose(m.data[:100000], np.round(m.data[:100000])):
        sys.exit(f"ERROR {s.donor}: matrix is not integer counts")
    obs = pd.DataFrame({"sample_id": s.donor, "library": s.library, "age": s.age,
                        "condition": "diseased" if s.group == "SLE" else "healthy"},
                       index=[f"{s.donor}::{b}" for b in barcodes])
    return ad.AnnData(m.astype(np.float32), obs=obs)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("src", type=Path); p.add_argument("studies", type=Path)
    p.add_argument("--samples", type=Path, default=Path(__file__).with_name("gse135779_samples.tsv"))
    a = p.parse_args()
    samples = pd.read_csv(a.samples, sep="\t")
    genes = pd.read_csv(a.src / "GSE135779_genes.tsv.gz", sep="\t", header=None)
    symbols = pd.Index(genes[1].astype(str)).to_numpy()
    for group, name in COHORTS.items():
        sub = samples[samples.age_group == group]
        parts = [read_sample(a.src / "raw", s, symbols) for s in sub.itertuples()]
        x = ad.concat(parts, join="outer", merge="same")
        x.var_names = pd.Index(symbols)
        x.var_names_make_unique()
        mito = np.asarray(x[:, x.var_names.str.startswith("MT-")].X.sum(1)).ravel()
        total = np.asarray(x.X.sum(1)).ravel()
        x.obs["pct_mt"] = 100 * mito / np.maximum(total, 1)
        x.obs["total_counts"] = total
        x.obs["cell_type"], clusters = annotate(x)
        qc = a.studies / name / "results"
        qc.mkdir(parents=True, exist_ok=True)
        clusters.to_csv(qc / "annotation_clusters.tsv", sep="\t", index=False, float_format="%.3f")
        by = x.obs.groupby(["sample_id", "condition"], observed=True)
        by.agg(cells=("age", "size"), median_counts=("total_counts", "median"), median_pct_mt=("pct_mt", "median"),
               age=("age", "first")).reset_index().to_csv(qc / "qc_per_sample.tsv", sep="\t", index=False, float_format="%.2f")
        comp = pd.crosstab(x.obs.sample_id, x.obs.cell_type, normalize="index")
        comp.insert(0, "condition", x.obs.groupby("sample_id", observed=True).condition.first().reindex(comp.index))
        comp.to_csv(qc / "composition_per_sample.tsv", sep="\t", float_format="%.4f")
        data = a.studies / name / "data"
        data.mkdir(parents=True, exist_ok=True)
        x.write_h5ad(data / f"{name}.h5ad")
        print(f"{name}: {x.n_obs} cells x {x.n_vars} genes; donors "
              f"{x.obs.groupby('condition', observed=True).sample_id.nunique().to_dict()}")
        print(x.obs.cell_type.value_counts().to_string())


if __name__ == "__main__":
    main()
