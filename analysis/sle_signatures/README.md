# SLE signatures per cell population

`signatures.tsv` holds SLE vs healthy genes from the lupus UCE test run: Perez 2022 cohort 4, 20 vs 20 donors,
pseudobulk DESeq2 on author cell-type labels, run 2026-10-29. Filter: padj < 0.05, with RPL/RPS/MT- genes
removed because of the unexplained ribosomal artefact. This is the same signature `scripts/l1000_reversal.py` uses.

Not included: cDC (0 significant genes), CD8 T (2 genes), B cells (no DE file from that run).
CD4 T has only 11 genes, so its scores are unstable.
