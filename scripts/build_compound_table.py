#!/usr/bin/env python
"""Build the seed compound table: SLE positive controls + compounds from the L1000 reversal lookup.

Structures from PubChem (by name or CID); ChEMBL ID and max clinical phase joined on InChIKey,
the key every chemical source shares (ChEMBL, PubChem, Repurposing Hub, L1000, OP3, Tahoe).
Every SMILES is parsed with RDKit; unparseable rows are kept and flagged, not dropped.
Usage: python scripts/build_compound_table.py docs/data/l1000_reversal_hits.csv shared/compounds/seed_compounds.tsv
Sends compound names, CIDs and InChIKeys to PubChem and ChEMBL; nothing else leaves the machine.
"""
import json, subprocess, sys, time, urllib.parse
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")
# Small-molecule SLE drugs and interferon-pathway inhibitors (planning/SLE_drug_discovery_plan.md).
# Antibodies (belimumab, anifrolumab, obinutuzumab, rituximab) have no SMILES and are excluded.
CONTROLS = {
    "hydroxychloroquine": "SLE standard of care", "prednisolone": "glucocorticoid",
    "methylprednisolone": "glucocorticoid", "dexamethasone": "glucocorticoid",
    "mycophenolic acid": "IMPDH inhibitor (active form of mycophenolate)", "azathioprine": "antimetabolite",
    "methotrexate": "antimetabolite", "cyclophosphamide": "alkylating agent",
    "voclosporin": "calcineurin inhibitor, approved for lupus nephritis", "tacrolimus": "calcineurin inhibitor",
    "deucravacitinib": "TYK2 inhibitor, phase 3 in SLE", "baricitinib": "JAK1/2 inhibitor",
    "tofacitinib": "JAK1/3 inhibitor", "ruxolitinib": "JAK1/2 inhibitor", "fedratinib": "JAK2 inhibitor (TG101348)",
}
PUG = "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound"
CHEMBL = "https://www.ebi.ac.uk/chembl/api/data/molecule.json"


def get(url):  # curl, not urllib: python.org builds on macOS ship without a CA bundle
    r = subprocess.run(["curl", "-fsS", "--max-time", "60", url], capture_output=True, text=True)
    time.sleep(0.25)  # PubChem asks for <= 5 requests/s
    return json.loads(r.stdout) if r.returncode == 0 else None


def pubchem(kind, ids):
    j = get(f"{PUG}/{kind}/{urllib.parse.quote(ids)}/property/SMILES,InChIKey,Title/JSON")
    return j["PropertyTable"]["Properties"] if j else []


rows = []
for name, role in CONTROLS.items():
    p = pubchem("name", name)
    rows.append(dict(name=name, source="positive_control", note=role, **({k: p[0].get(k) for k in ("CID", "SMILES", "InChIKey")} if p else {})))
hits = pd.read_csv(sys.argv[1]).dropna(subset=["pubchem"])
cids = sorted({str(int(c)) for c in hits.pubchem})
props = {str(p["CID"]): p for i in range(0, len(cids), 50) for p in pubchem("cid", ",".join(cids[i:i + 50]))}
for cid in cids:
    p = props.get(cid, {})
    h = hits[hits.pubchem.astype(int).astype(str) == cid].iloc[0]
    rows.append(dict(name=h.compound if isinstance(h.compound, str) else (p.get("Title") or h.pert_id), source="l1000_hit",
                     note=f"{h.pert_id}; best {hits[hits.pubchem.astype(int).astype(str) == cid].score.max():.3f}",
                     CID=p.get("CID", int(cid)), SMILES=p.get("SMILES"), InChIKey=p.get("InChIKey")))

df = pd.DataFrame(rows)
keys = df.InChIKey.dropna().unique().tolist()
chembl = {}
for i in range(0, len(keys), 40):
    j = get(f"{CHEMBL}?molecule_structures__standard_inchi_key__in={','.join(keys[i:i + 40])}&limit=100")
    for m in (j or {}).get("molecules", []):
        chembl[m["molecule_structures"]["standard_inchi_key"]] = (m["molecule_chembl_id"], m.get("max_phase"))
df["chembl_id"] = df.InChIKey.map(lambda k: chembl.get(k, (None, None))[0])
df["max_phase"] = df.InChIKey.map(lambda k: chembl.get(k, (None, None))[1])
df["rdkit_ok"] = df.SMILES.map(lambda s: isinstance(s, str) and Chem.MolFromSmiles(s) is not None)
df = df.rename(columns={"CID": "pubchem_cid", "SMILES": "smiles", "InChIKey": "inchikey"})
df.to_csv(sys.argv[2], sep="\t", index=False)
print(f"{len(df)} compounds ({(df.source == 'positive_control').sum()} controls): smiles {df.smiles.notna().sum()}, "
      f"rdkit_ok {df.rdkit_ok.sum()}, in ChEMBL {df.chembl_id.notna().sum()}, approved (phase 4) {(df.max_phase.astype(float) == 4).sum()}")
missing = df[~df.rdkit_ok].name.tolist()
if missing:
    print("WARNING: no parseable structure for:", ", ".join(map(str, missing)))
