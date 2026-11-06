#!/usr/bin/env python
"""Lupus treatment status per compound, from shared/compounds/lupus_use.tsv, the one definition every page reads.

Compounds match on the InChIKey connectivity block of their largest fragment, so salts, stereo labels and
vendor names (FK506, Vincristine sulfate) resolve to one drug. No match means not found in the 2026-11-05
check, not "never used".
Run directly to verify every entry resolves to a compound we hold.
"""
from pathlib import Path

import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem.inchi import MolToInchiKey

RDLogger.DisableLog("rdApp.*")
ROOT = Path(__file__).parent.parent
TABLE = ROOT / "shared/compounds/lupus_use.tsv"
LABELS = {"current": "Current lupus use", "historical": "Historical lupus use", "trial_only": "Lupus trials only"}


def block(smiles):
    m = Chem.MolFromSmiles(str(smiles))
    if m is None:
        return None
    parts = Chem.GetMolFrags(m, asMols=True)
    return MolToInchiKey(max(parts, key=lambda x: x.GetNumHeavyAtoms()))[:14]


def load():
    t = pd.read_csv(TABLE, sep="\t")
    assert t.block.is_unique and set(t.status) <= set(LABELS), "bad lupus_use.tsv"
    return t.set_index("block")


def lookup(smiles, table=None):
    """Return a dict for the compound, or None."""
    t = load() if table is None else table
    b = block(smiles)
    if b not in t.index:
        return None
    r = t.loc[b]
    return dict(status=r.status, label=LABELS[r.status], drug=r["name"], basis=r.basis, source=r.source)


if __name__ == "__main__":
    t = load()
    held = {block(s) for s in pd.read_csv(ROOT / "shared/compounds/seed_compounds.tsv", sep="\t").smiles.dropna()}
    held |= {block(s) for s in pd.read_csv(ROOT / "analysis/op3_rescue/rescue.tsv", sep="\t").SMILES}
    missing = t[~t.index.isin(held)]
    assert missing.empty, f"entries match no compound: {missing['name'].tolist()}"
    assert lookup("CN1CCN(CC1)C")is None
    print(f"lupus_use OK: {len(t)} entries, {t.status.value_counts().to_dict()}")
