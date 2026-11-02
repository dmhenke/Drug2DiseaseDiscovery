#!/usr/bin/env python
"""Does the single-cell pseudobulk SLE signature replicate in the bulk whole-blood meta-analysis?
Expectation: interferon genes replicate in every cell type; most significant single-cell genes that are
on the bulk arrays are significant in bulk with the same sign. Whole blood includes neutrophils and
mixes cell types, so per-gene fold changes should correlate only moderately.
Usage: python analysis/sle_bulk_meta/compare_with_scrna.py <scRNA DE results dir> [out.tsv]
"""
import re, sys
from pathlib import Path
import numpy as np
import pandas as pd

HERE = Path(__file__).parent
bulk = pd.read_csv(HERE / "outputs/merged_sle_diff_expr.csv")
bulk = bulk[~bulk.gene.str.match(r"^NA(\.|$)")]  # probes without a symbol, merged under "NA" upstream
res = Path(sys.argv[1])
rows = []
for f in sorted(res.glob("DE_human_*.csv")):
    if f.name == "DE_all_human.csv":
        continue
    sc = pd.read_csv(f).dropna(subset=["padj"])
    m = bulk.merge(sc, on="gene")
    sig = sc[sc.padj < 0.05].merge(bulk, on="gene", how="left")
    on_array = sig["adj.P.Val"].notna()
    rep = on_array & (sig["adj.P.Val"] < 0.05) & (np.sign(sig.logFC) == np.sign(sig.log2FoldChange))
    rows.append(dict(cell_type=f.stem.replace("DE_human_", ""), genes_shared=len(m),
                     spearman_lfc=round(m.logFC.corr(m.log2FoldChange, method="spearman"), 2) if len(m) > 2 else np.nan,
                     sc_sig=len(sig), sc_sig_on_bulk_array=int(on_array.sum()), replicated=int(rep.sum()),
                     not_on_array=";".join(sig.gene[~on_array].head(12))))
out = pd.DataFrame(rows)
print(out.to_string(index=False))
r = bulk[bulk.gene.str.match(r"^(RPL|RPS)")]
print(f"\nbulk ribosomal RPL/RPS genes: {len(r)}, down at adj.P<0.05: {((r['adj.P.Val'] < 0.05) & (r.logFC < 0)).sum()}, "
      f"up: {((r['adj.P.Val'] < 0.05) & (r.logFC > 0)).sum()}")
if len(sys.argv) > 2:
    out.to_csv(sys.argv[2], sep="\t", index=False)
