#!/usr/bin/env python
"""Recorded check for the synthetic T0 fixture. Expectations were fixed before the run:
  - recall of planted genes at padj < 0.05 >= 0.75 in every cell type
  - non-planted calls at padj < 0.05 <= 25 per cell type (about 5% of 460 nulls)
Prints a row per cell type. Exit 1 on any failure, or if results are missing (never a silent pass).

Usage: python scripts/check_t0.py studies/synthetic_t0
"""
import os
import sys

import pandas as pd

d = sys.argv[1]
res_path = os.path.join(d, "results", "DE_all_human.csv")
if not os.path.exists(res_path):
    sys.exit(f"FAIL: {res_path} missing; Stage 3 produced no results")
res = pd.read_csv(res_path)
planted = set(pd.read_csv(os.path.join(d, "data", "planted_truth.tsv"), sep="\t").gene)

ok = True
print(f"{'cell_type':<16}{'tested':>8}{'recall':>9}{'false':>7}  verdict")
for ct, g in res.groupby("cell_type"):
    called = set(g.loc[g.padj < 0.05, "gene"])
    recall = len(called & planted) / len(planted)
    false = len(called - planted)
    passed = recall >= 0.75 and false <= 25
    ok = ok and passed
    print(f"{ct:<16}{len(g):>8}{recall:>9.2f}{false:>7}  {'PASS' if passed else 'FAIL'}")
print("RESULT:", "PASS" if ok else "FAIL")
sys.exit(0 if ok else 1)
