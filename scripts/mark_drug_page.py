#!/usr/bin/env python
"""Write lupus treatment status into docs/drug-reversal.html and docs/data/l1000_reversal_hits.csv.

L1000 hits carry a PubChem CID; shared/compounds/seed_compounds.tsv maps CID to structure. Safe to rerun.
"""
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from lupus_use import ROOT, load, lookup  # noqa: E402

seed = pd.read_csv(ROOT / "shared/compounds/seed_compounds.tsv", sep="\t").dropna(subset=["smiles", "pubchem_cid"])
smiles = dict(zip(seed.pubchem_cid.astype(float).astype(int).astype(str), seed.smiles))
hits_path = ROOT / "docs/data/l1000_reversal_hits.csv"
hits = pd.read_csv(hits_path, dtype={"pubchem": str})
table = load()
use = {c: lookup(smiles[c], table) for c in hits.pubchem.dropna().unique() if c in smiles}
hits["lupus_use"] = hits.pubchem.map(lambda c: use[c]["status"] if use.get(c) else "")
hits.to_csv(hits_path, index=False)
unchecked = hits[hits.pubchem.isna() | ~hits.pubchem.isin(smiles)]
assert not hits.compound.str.contains("FK506|sirolimus", na=False).any() or (hits.lupus_use == "current").any()

page = ROOT / "docs/drug-reversal.html"
data = json.dumps({c: u for c, u in use.items() if u}, separators=(",", ":")).replace("</", "<\\/")
html, n = re.subn(r"/\*LUPUS_USE\*/.*?/\*END\*/", lambda _: f"/*LUPUS_USE*/const LUPUS_USE = {data};/*END*/", page.read_text(), flags=re.S)
assert n == 1, "marker block not found"
page.write_text(html)
marked = hits[hits.lupus_use != ""]
print(f"{len(marked)} of {len(hits)} hits marked; {len(unchecked)} rows have no structure to check")
print(marked.groupby(["compound", "lupus_use"], dropna=False).size().to_string())
