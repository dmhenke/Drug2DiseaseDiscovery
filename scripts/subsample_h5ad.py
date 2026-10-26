#!/usr/bin/env python
"""Subsample a CELLxGENE-style h5ad by donor. Reads in backed mode; the full file is never loaded.

Picks donors per condition (seeded), caps cells per donor, writes raw integer counts to X,
and writes OUT.manifest.json. Run with --dry-run first to see donor and cell counts.

Usage:
  python scripts/subsample_h5ad.py SRC.h5ad OUT.h5ad --donor-col donor_id --condition-col disease \\
      --healthy normal --diseased "systemic lupus erythematosus" \\
      --donors-per-condition 5 --max-cells-per-donor 2000 [--counts raw|X] [--seed 1234] [--dry-run]
  --all-donors   use every donor (full cohort)
  --max-cells-per-donor 0   no cap
"""
import argparse
import datetime
import json
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp


def main():
    p = argparse.ArgumentParser()
    p.add_argument("src")
    p.add_argument("out")
    p.add_argument("--donor-col", default="donor_id")
    p.add_argument("--condition-col", default="disease")
    p.add_argument("--healthy", required=True)
    p.add_argument("--diseased", required=True)
    p.add_argument("--donors-per-condition", type=int, default=5)
    p.add_argument("--all-donors", action="store_true")
    p.add_argument("--max-cells-per-donor", type=int, default=2000, help="0 = no cap")
    p.add_argument("--counts", choices=["raw", "X"], default="raw")
    p.add_argument("--seed", type=int, default=1234)
    p.add_argument("--dry-run", action="store_true")
    a = p.parse_args()

    src = ad.read_h5ad(a.src, backed="r")
    donor = src.obs[a.donor_col].astype(str).to_numpy()
    raw_cond = src.obs[a.condition_col].astype(str).to_numpy()
    cond = np.where(raw_cond == a.healthy, "healthy",
                    np.where(raw_cond == a.diseased, "diseased", ""))
    df = pd.DataFrame({"donor": donor, "cond": cond})
    n_labels = df[df.cond != ""].groupby("donor")["cond"].nunique()
    mixed = list(n_labels[n_labels > 1].index)
    if mixed:
        print(f"WARNING: dropping {len(mixed)} donors that carry both conditions", file=sys.stderr)
    ok = (cond != "") & ~np.isin(donor, mixed)
    rng = np.random.default_rng(a.seed)

    chosen = {}
    for c in ["healthy", "diseased"]:
        ds = sorted(set(donor[ok & (cond == c)]))
        if a.all_donors:
            chosen[c] = ds
        else:
            if len(ds) < a.donors_per_condition:
                sys.exit(f"ERROR: {len(ds)} {c} donors available, need {a.donors_per_condition}")
            chosen[c] = sorted(rng.choice(ds, a.donors_per_condition, replace=False).tolist())

    keep = []
    for c in chosen:
        for d in chosen[c]:
            idx = np.flatnonzero(ok & (donor == d))
            if a.max_cells_per_donor and len(idx) > a.max_cells_per_donor:
                idx = np.sort(rng.choice(idx, a.max_cells_per_donor, replace=False))
            keep.append(idx)
    keep = np.sort(np.concatenate(keep))

    cond_of = {d: c for c in chosen for d in chosen[c]}
    print(f"donors: healthy={len(chosen['healthy'])} diseased={len(chosen['diseased'])} cells={len(keep)}")
    for d, n in pd.Series(donor[keep]).value_counts().sort_index().items():
        print(f"  {cond_of[d]:<9} {d:<40} {n:>7}")
    if a.dry_run:
        return

    if a.counts == "raw":
        if src.raw is None:
            sys.exit("ERROR: source has no .raw; use --counts X only if X holds integer counts")
        X, var = src.raw.X, src.raw.var
    else:
        X, var = src.X, src.var
    sub = X[keep]
    sub = sub.tocsr() if sp.issparse(sub) else sp.csr_matrix(sub)
    head = sub[:200]
    if (head.data < 0).any() or not np.allclose(head.data, np.round(head.data)):
        sys.exit("ERROR: selected matrix is not non-negative integer counts")

    out = ad.AnnData(X=sub.astype(np.float32), obs=src.obs.iloc[keep].copy(), var=var.copy())
    out.write_h5ad(a.out)
    man = {
        "source": os.path.abspath(a.src), "source_bytes": os.path.getsize(a.src),
        "output": os.path.abspath(a.out), "cells": int(len(keep)),
        "counts_from": a.counts, "seed": a.seed, "max_cells_per_donor": a.max_cells_per_donor,
        "donors": chosen,
        "created": datetime.datetime.now().isoformat(timespec="seconds"),
    }
    with open(a.out + ".manifest.json", "w") as f:
        json.dump(man, f, indent=2)
    print(f"wrote {a.out}: {len(keep)} cells; manifest {a.out}.manifest.json")


if __name__ == "__main__":
    main()
