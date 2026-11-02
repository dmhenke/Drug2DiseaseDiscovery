#!/usr/bin/env python
"""Content checks on the OP3 download, with the biological expectation beside each one.

A file is reusable only if it passes these; the fetch script checks size only.
Expectations (written before running, from Szalata et al. NeurIPS 2024 D&B):
  - 4 cell types (T, B, NK, Myeloid), 3 donors, ~140 compounds, 5,317 genes after their filtering.
  - Every compound has a SMILES that RDKit parses.
  - Positive controls belinostat (HDAC) and dabrafenib (BRAF) are present. DMSO is the reference
    these DE values are computed against, so it must NOT appear as a compound row.
  - Glucocorticoid response: prednisolone induces TSC22D3 and FKBP5 in myeloid cells (it lands in
    the test split). If this fails, the labels or the DE sign convention are wrong, not the biology.
Usage: python scripts/check_op3.py shared/perturbation/op3
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError:  # record it, never skip silently; the other checks are independent
    Chem = None
d = Path(sys.argv[1] if len(sys.argv) > 1 else "shared/perturbation/op3")
fails = []


def check(name, ok, detail=""):
    print(f"  [{'ok' if ok else 'FAIL'}] {name}{': ' + detail if detail else ''}")
    if not ok:
        fails.append(name)


for f in ["de_train.h5ad", "de_test.h5ad"]:
    print(f"== {f}")
    a = ad.read_h5ad(d / f)
    obs = a.obs
    check("expected obs columns", {"cell_type", "sm_name", "SMILES"} <= set(obs.columns), ", ".join(obs.columns[:8]))
    check("genes", a.n_vars > 4000, f"{a.n_vars} genes x {a.n_obs} (compound, cell type) rows")
    print("       cell types:", ", ".join(map(str, obs.cell_type.unique())))
    sm = obs.drop_duplicates("sm_name")[["sm_name", "SMILES"]]
    if Chem is None:
        check("every SMILES parses", False, f"{len(sm)} compounds NOT CHECKED: no rdkit (pip install rdkit)")
    else:
        bad = [n for n, s in zip(sm.sm_name, sm.SMILES) if Chem.MolFromSmiles(str(s)) is None]
        check("every SMILES parses", not bad, f"{len(sm)} compounds" + (f"; bad: {bad[:5]}" if bad else ""))

# Compounds live in whichever split they were assigned to, so look across both.
tr, te = ad.read_h5ad(d / "de_train.h5ad"), ad.read_h5ad(d / "de_test.h5ad")
names = set(tr.obs.sm_name.astype(str)) | set(te.obs.sm_name.astype(str))
for c in ["Belinostat", "Dabrafenib"]:
    check(f"positive control present: {c}", any(c.lower() in n.lower() for n in names))
# DMSO is the reference these DE values are computed against, so it must NOT be a row of its own.
check("DMSO is the reference, not a row", not any("dimethyl sulfoxide" in n.lower() for n in names),
      f"{len(names)} compounds across both splits")
hits = sorted(n for n in names if any(k in n.lower() for k in
              ("ruxolitinib", "prednisolone", "dexameth", "mometasone", "tofacitinib", "baricitinib", "methotrexate")))
print("  interferon/steroid axis compounds:", ", ".join(hits) or "NONE")
check("a JAK inhibitor is present", any("ruxolitinib" in h.lower() or "tinib" in h.lower() for h in hits))
check("a glucocorticoid is present", any(k in " ".join(hits).lower() for k in ("prednisolone", "dexameth", "mometasone")))

# Glucocorticoid sanity check: canonical direct targets of the receptor. Positive DE = induced.
for split, a in (("train", tr), ("test", te)):
    row = (a.obs.sm_name.astype(str).str.lower().str.contains("prednisolone")
           & a.obs.cell_type.astype(str).str.contains("Myeloid"))
    if not row.any():
        continue
    v = np.asarray(a[row].X).ravel()
    g = {n: i for i, n in enumerate(a.var_names)}
    vals = {k: float(v[g[k]]) for k in ("TSC22D3", "FKBP5") if k in g}
    check(f"prednisolone induces TSC22D3/FKBP5 in myeloid ({split})",
          bool(vals) and all(x > 0 for x in vals.values()),
          str({k: round(x, 2) for k, x in vals.items()}) + " (signed DE statistic)")
    break
else:
    check("prednisolone x myeloid row found in either split", False, "compound or cell-type label changed")

print("\nOP3 CHECK:", "FAIL -> " + ", ".join(fails) if fails else "PASS")
sys.exit(1 if fails else 0)
