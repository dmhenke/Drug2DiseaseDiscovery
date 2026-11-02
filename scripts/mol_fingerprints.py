#!/usr/bin/env python
"""SMILES -> Morgan fingerprints, keyed by InChIKey. Stage C1 of planning/chemical_perturbation_model.md.

Morgan r=2, 2048 bits: no extra model, deterministic, and the baseline a pretrained encoder has to beat.
Reads a TSV/CSV with `smiles` and either `inchikey` or `name`; writes <out>.npy (float32, rows in the
order of <out>.index.tsv). Run with --self-check to run the analogue-pair test only.

Expectation: known analogue pairs score higher than unrelated pairs. If that fails, the structures or
the fingerprint settings are wrong.
Usage: python scripts/mol_fingerprints.py in.tsv out_prefix [--radius 2] [--bits 2048]
       python scripts/mol_fingerprints.py --self-check
"""
import argparse
import sys

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator as fpgen
from rdkit.Chem.inchi import MolToInchiKey

RDLogger.DisableLog("rdApp.*")

# (a, b, c): similarity(a, b) must exceed similarity(a, c). Steroid backbone, then JAK inhibitors.
PAIRS = [
    ("C[C@]12C[C@@H]([C@H]3[C@H]([C@@H]1CC[C@@]2(C(=O)CO)O)CCC4=CC(=O)C=C[C@]34C)O",   # prednisolone
     "C[C@H]1C[C@H]2[C@@H]3CC[C@@]([C@]3(C[C@@H]([C@@H]2[C@@]4(C1=CC(=O)C=C4)C)O)C)(C(=O)CO)O",  # methylprednisolone
     "CN(CC1=CN=C2C(=N1)C(=NC(=N2)N)N)C3=CC=C(C=C3)C(=O)N[C@@H](CCC(=O)O)C(=O)O"),      # methotrexate
    ("CCS(=O)(=O)N1CC(C1)(CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3",                           # baricitinib
     "C1CCC(C1)[C@@H](CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3",                               # ruxolitinib
     "CCN(CCCC(C)NC1=C2C=CC(=CC2=NC=C1)Cl)CCO"),                                        # hydroxychloroquine
]


def fingerprints(smiles, radius=2, bits=2048):
    """Return (array, keys, bad_indices). Unparseable SMILES give an all-zero row, never a dropped one."""
    gen = fpgen.GetMorganGenerator(radius=radius, fpSize=bits)
    rows, keys, bad = np.zeros((len(smiles), bits), dtype=np.float32), [], []
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(str(s))
        if m is None:
            bad.append(i); keys.append(None); continue
        DataStructs.ConvertToNumpyArray(gen.GetFingerprint(m), rows[i])
        keys.append(MolToInchiKey(m))
    return rows, keys, bad


def self_check(radius=2, bits=2048):
    ok = True
    for a, b, c in PAIRS:
        f, _, _ = fingerprints([a, b, c], radius, bits)
        tan = lambda x, y: float(x @ y) / float((x + y > 0).sum())  # noqa: E731  (binary Tanimoto)
        near, far = tan(f[0], f[1]), tan(f[0], f[2])
        print(f"  analogue {near:.2f} vs unrelated {far:.2f}  {'ok' if near > far else 'FAIL'}")
        ok &= near > far
    print("FINGERPRINT CHECK:", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("table", nargs="?"); p.add_argument("out", nargs="?")
    p.add_argument("--radius", type=int, default=2); p.add_argument("--bits", type=int, default=2048)
    p.add_argument("--self-check", action="store_true")
    a = p.parse_args()
    if a.self_check or not a.table:
        sys.exit(0 if self_check(a.radius, a.bits) else 1)
    if not self_check(a.radius, a.bits):
        sys.exit("refusing to write fingerprints: analogue check failed")
    df = pd.read_csv(a.table, sep="\t" if a.table.endswith((".tsv", ".txt")) else ",")
    X, keys, bad = fingerprints(df.smiles, a.radius, a.bits)
    df = df.assign(inchikey_rdkit=keys, fingerprint_ok=[i not in set(bad) for i in range(len(df))])
    np.save(a.out + ".npy", X)
    df.to_csv(a.out + ".index.tsv", sep="\t", index=False)
    print(f"wrote {a.out}.npy {X.shape} and {a.out}.index.tsv; {len(bad)} unparseable (zero rows, flagged)")
