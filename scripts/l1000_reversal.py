"""Query L1000CDS2 (reverse mode) with each cell type's SLE-vs-healthy signature from the lupus UCE run.
Signature: padj<0.05, ribosomal/mito genes removed (unexplained artefact). Up in SLE -> up list.
Calibration: 5 random signatures per cell type, same up/down sizes, drawn from that cell type's
non-significant tested genes (seed 1234); a real signature should out-score them.
Expectation, written before the run: classical monocytes (interferon-dominated) return a top-50
whose scores exceed every random set, and at least one JAK inhibitor appears among named hits.
Usage: python scripts/l1000_reversal.py <DE results dir> <out dir>   (writes l1000_hits.csv, l1000_calibration.csv)
Sends gene symbols to the public L1000CDS2 server (maayanlab.cloud); nothing else leaves the machine."""
import json, re, subprocess, sys, time
from pathlib import Path
import numpy as np
import pandas as pd

DL = Path(sys.argv[1])
HERE = Path(sys.argv[2])
URL = "https://maayanlab.cloud/L1000CDS2/query"
TYPES = {
    "Classical mono": "DE_human_classical_monocyte.csv",
    "Non-classical mono": "DE_human_non-classical_monocyte.csv",
    "NK": "DE_human_natural_killer_cell.csv",
    "CD4 T": "DE_human_CD4-positive,_alpha-beta_T_cell.csv",
}
RIBO = re.compile(r"^(RPL|RPS|MT-)")
rng = np.random.default_rng(1234)


def query(up, dn):
    body = {"data": {"upGenes": up, "dnGenes": dn},
            "config": {"aggravate": False, "searchMethod": "geneSet", "share": False,
                       "combination": False, "db-version": "latest"}}
    # curl, not urllib: python.org builds on macOS ship without a CA bundle
    r = subprocess.run(["curl", "-fsS", "--max-time", "120", "-X", "POST", URL, "-H", "Content-Type: application/json",
                        "-d", json.dumps(body)], capture_output=True, text=True, check=True)
    out = json.loads(r.stdout)
    time.sleep(1)
    return out["topMeta"]


hits, calib = [], []
for ct, f in TYPES.items():
    d = pd.read_csv(DL / f).dropna(subset=["padj"])
    d = d[~d.gene.str.match(RIBO)]
    sig = d[d.padj < 0.05]
    up = sig[sig.log2FoldChange > 0].gene.tolist(); dn = sig[sig.log2FoldChange < 0].gene.tolist()
    top = query(up, dn)
    for rank, m in enumerate(top, 1):
        ov = m.get("overlap", {})
        hits.append(dict(cell_type=ct, rank=rank, n_up=len(up), n_dn=len(dn), score=m["score"],
                         compound=m["pert_desc"] if m["pert_desc"] != "-666" else "",
                         pert_id=m["pert_id"], pubchem=m.get("pubchem_id", ""), cell_line=m["cell_id"],
                         dose=f'{float(m["pert_dose"]):g} {m["pert_dose_unit"]}',
                         time=f'{float(m["pert_time"]):g} {m["pert_time_unit"]}',
                         reversed_up=";".join(ov.get("up/dn", [])), reversed_dn=";".join(ov.get("dn/up", []))))
    real = [m["score"] for m in top]
    bg = d[d.padj >= 0.05].gene.to_numpy()
    for k in range(5):
        pick = rng.choice(bg, len(up) + len(dn), replace=False)
        s = [m["score"] for m in query(pick[:len(up)].tolist(), pick[len(up):].tolist())]
        calib.append(dict(cell_type=ct, rep=k, top1=max(s), median50=float(np.median(s))))
    print(f"{ct}: up={len(up)} dn={len(dn)} top1={max(real):.3f} median50={np.median(real):.3f} | random top1 "
          f"{[c['top1'] for c in calib if c['cell_type'] == ct]}", flush=True)

pd.DataFrame(hits).to_csv(HERE / "l1000_hits.csv", index=False)
pd.DataFrame(calib).to_csv(HERE / "l1000_calibration.csv", index=False)
