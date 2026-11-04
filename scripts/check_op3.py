#!/usr/bin/env python
"""Content checks on the OP3 download, with the biological expectation beside each one.

A file is reusable only if it passes these; the fetch script checks size only.
Expectations (written before running, from Szalata et al. NeurIPS 2024 D&B):
  - 4 cell types (T, B, NK, Myeloid), ~140 compounds, 5,317 genes after their filtering.
    Donors are not checked: the DE tables have no donor column.
  - Every compound has a SMILES that RDKit parses.
  - Positive controls belinostat (HDAC) and dabrafenib (BRAF) are present. DMSO is the reference
    these DE values are computed against, so it must NOT appear as a compound row.
  - Glucocorticoid response: prednisolone induces TSC22D3 and FKBP5 in myeloid cells (it lands in
    the test split). If this fails, the labels or the DE sign convention are wrong, not the biology.
Usage: python scripts/check_op3.py shared/perturbation/op3
Writes <dir>/check/{check_report.txt,split_by_cell_type.tsv,compounds.tsv,glucocorticoid_targets.tsv}.
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

try:
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
except ImportError:  # record it, never skip silently; the other checks are independent
    Chem = None
d = Path(sys.argv[1] if len(sys.argv) > 1 else "shared/perturbation/op3")
out = d / "check"
out.mkdir(exist_ok=True)
report = open(out / "check_report.txt", "w")
fails = []


def say(line):
    print(line)
    report.write(line + "\n")


def check(name, ok, detail=""):
    say(f"  [{'ok' if ok else 'FAIL'}] {name}{': ' + detail if detail else ''}")
    if not ok:
        fails.append(name)


for f in ["de_train.h5ad", "de_test.h5ad"]:
    say(f"== {f}")
    a = ad.read_h5ad(d / f)
    obs = a.obs
    check("expected obs columns", {"cell_type", "sm_name", "SMILES"} <= set(obs.columns), ", ".join(obs.columns[:8]))
    check("genes", a.n_vars > 4000, f"{a.n_vars} genes x {a.n_obs} (compound, cell type) rows")
    say("       cell types: " + ", ".join(map(str, obs.cell_type.unique())))
    say("       X: " + ("absent" if a.X is None else "present") + "; layers: " + (", ".join(a.layers.keys()) or "none"))
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
say("  interferon/steroid axis compounds: " + (", ".join(hits) or "NONE"))
check("a JAK inhibitor is present", any("ruxolitinib" in h.lower() or "tinib" in h.lower() for h in hits))
check("a glucocorticoid is present", any(k in " ".join(hits).lower() for k in ("prednisolone", "dexameth", "mometasone")))

# What the split really holds out decides how Stage P must split: test compounds that also have
# train rows in another cell type are not unseen chemistry.
obs = pd.concat([tr.obs.assign(file="train"), te.obs.assign(file="test")]).astype({"sm_name": str, "cell_type": str})
obs.groupby(["cell_type", "file"]).sm_name.nunique().unstack(fill_value=0).to_csv(out / "split_by_cell_type.tsv", sep="\t")
cols = {f: obs[obs.file == f].groupby("sm_name").cell_type.agg(lambda x: ",".join(sorted(x))) for f in ("train", "test")}
pd.DataFrame({"SMILES": obs.drop_duplicates("sm_name").set_index("sm_name").SMILES,
              "train_cell_types": cols["train"], "test_cell_types": cols["test"]}).to_csv(out / "compounds.tsv", sep="\t")
shared = set(cols["train"].index) & set(cols["test"].index)
say(f"  {len(shared)} of {len(cols['test'])} test compounds also have train rows (not unseen chemistry)")

# Glucocorticoid sanity check: canonical direct targets of the receptor. Positive = induced.
# The DE values may sit in layers rather than X; sign_log10_pval and logFC are both signed.
TARGETS = ("TSC22D3", "FKBP5")
rows = []
for split, a in (("train", tr), ("test", te)):
    layers = {k: a.layers[k] for k in ("logFC", "sign_log10_pval") if k in a.layers}
    if a.X is not None:
        layers["X"] = a.X
    names = a.obs.sm_name.astype(str).str.lower()
    for i in np.flatnonzero(names.str.contains("prednisolone|mometasone")):
        for gene in TARGETS:
            if gene in a.var_names:
                j = a.var_names.get_loc(gene)
                rows.append({"split": split, "sm_name": a.obs.sm_name.iloc[i], "cell_type": a.obs.cell_type.iloc[i],
                             "gene": gene, **{k: float(np.asarray(v[i, j].todense() if hasattr(v[i, j], "todense") else v[i, j]))
                                              for k, v in layers.items()}})
gc = pd.DataFrame(rows)
gc.to_csv(out / "glucocorticoid_targets.tsv", sep="\t", index=False)
sel = gc[gc.sm_name.str.lower().str.contains("prednisolone") & gc.cell_type.str.contains("Myeloid")] if len(gc) else gc
signed = [c for c in ("logFC", "sign_log10_pval", "X") if c in sel.columns]
if sel.empty or not signed:
    check("prednisolone x myeloid values found", False, "no row, no target gene, or no signed DE layer")
else:
    check("prednisolone induces TSC22D3/FKBP5 in myeloid",
          set(sel.gene) == set(TARGETS) and bool((sel[signed[0]] > 0).all()),
          "; ".join(f"{g} {signed[0]}={v:.2f} ({f})" for g, v, f in zip(sel.gene, sel[signed[0]], sel.split)))

say("\nOP3 CHECK: " + ("FAIL -> " + ", ".join(fails) if fails else "PASS"))
report.close()
sys.exit(1 if fails else 0)
