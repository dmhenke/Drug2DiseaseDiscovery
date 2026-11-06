#!/usr/bin/env python
"""Planted-failure test of lincs_rescue.py, needs no download. Expected outcome per planted case is asserted below."""
import gzip
import subprocess
import sys
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from lincs_rescue import spec  # noqa: E402

G = 80


def tsv(path, df):
    with gzip.open(path, "wt") as fh:
        df.to_csv(fh, sep="\t", index=False)


def build(root):
    rng = np.random.default_rng(1)
    meta, op3 = root / "meta", root / "op3"
    meta.mkdir(); op3.mkdir()
    smiles = pd.read_csv(HERE.parent / "analysis/op3_rescue/rescue.tsv", sep="\t").drop_duplicates("sm_name").SMILES.to_list()[:100]
    names = [f"cpd{i}" for i in range(len(smiles))]
    genes = [f"G{i}" for i in range(G)]
    # OP3: Myeloid and NK carry signal, T cells are what CD4 T is scored on
    cells = ["Myeloid cells", "NK cells", "T cells", "B cells"]
    obs = pd.DataFrame([(n, c, s) for n, s in zip(names, smiles) for c in cells], columns=["sm_name", "cell_type", "SMILES"])
    Y = rng.normal(0, 1, (len(obs), G))
    a = ad.AnnData(np.zeros(Y.shape), obs=obs.set_index(obs.index.astype(str)), var=pd.DataFrame(index=genes), layers={"logFC": Y})
    a.write_h5ad(op3 / "de_train.h5ad")
    a[:0].copy().write_h5ad(op3 / "de_test.h5ad")
    myeloid = pd.DataFrame(Y[obs.cell_type == "Myeloid cells"], index=names, columns=genes)
    # DE: monocytes on G0-59; NK and the rest below the 20-gene floor must be refused, not scored
    d = rng.normal(0, 1, G)
    rows = [dict(gene=g, log2FoldChange=d[i], padj=0.001 if i < 60 else 0.9, lfcSE=.1, pvalue=.001, cell_type=c)
            for c in ("classical monocyte", "non-classical monocyte") for i, g in enumerate(genes)]
    rows += [dict(gene=g, log2FoldChange=d[i], padj=0.001 if i < 5 else 0.9, lfcSE=.1, pvalue=.001, cell_type="natural killer cell")
             for i, g in enumerate(genes)]
    pd.DataFrame(rows).to_csv(root / "DE_all_human.csv", index=False)
    # LINCS: 24 h faithful and 6 h random in THP1, random MCF7, 30 shared compounds in U937
    all_names = names + [f"extra{i}" for i in range(50)]
    sig, prof = [], {}
    for line, keep in (("THP1", all_names), ("MCF7", all_names), ("U937", names[:30] + all_names[100:])):
        for c in keep:
            for t in (6, 24):
                s = f"{line}:{c}:{t}"
                sig.append(dict(sig_id=s, pert_type="trt_cp", cell_id=line, pert_iname=c, pert_id=f"BRD-{c}", pert_dose=10.0, pert_time=t))
                real = line in ("THP1", "U937") and t == 24 and c in myeloid.index
                prof[s] = myeloid.loc[c].to_numpy() * 3 + rng.normal(0, .3, G) if real else rng.normal(0, 3, G)
    sig = pd.DataFrame(sig)
    tsv(meta / "GSE92742_Broad_LINCS_sig_info.txt.gz", sig)
    tsv(meta / "GSE92742_Broad_LINCS_gene_info.txt.gz", pd.DataFrame(dict(pr_gene_id=[str(1000 + i) for i in range(G)], pr_gene_symbol=genes, pr_is_lm=1)))
    tsv(meta / "GSE92742_Broad_LINCS_pert_info.txt.gz", pd.DataFrame(dict(pert_id=[f"BRD-{c}" for c in names], canonical_smiles=smiles)))
    with h5py.File(root / "x.gctx", "w") as f:
        f["0/DATA/0/matrix"] = np.array([prof[s] for s in sig.sig_id], dtype=np.float32)
        f["0/META/ROW/id"] = np.array([str(1000 + i) for i in range(G)], dtype="S")
        f["0/META/COL/id"] = np.array(sig.sig_id.to_list(), dtype="S")


def run(root, out, *extra, de="DE_all_human.csv"):
    return subprocess.run([sys.executable, str(HERE / "lincs_rescue.py"), str(root / "x.gctx"), str(root / "meta"), str(root / out),
                           "--de", str(root / de), "--op3", str(root / "op3"), "--lines", "THP1", "MCF7", "U937", "NOMO1", *extra],
                          capture_output=True, text=True)


def main():
    d = np.array([2.0, -1.0, 0.5])
    assert np.isclose(spec(-d, 0 * d, d), 1) and np.isclose(spec(d, 0 * d, d), -1)
    with tempfile.TemporaryDirectory() as td:
        root = Path(td)
        build(root)
        missing = run(root, "none", de="absent.csv")
        assert missing.returncode != 0 and "missing" in missing.stderr, "missing DE must stop the run, not fall back"
        r = run(root, "t24", "--time", "24")
        assert r.returncode == 0, r.stderr
        c = pd.read_csv(root / "t24/calibration.tsv", sep="\t").set_index(["cell_line", "population"])
        assert c.loc[("THP1", "Classical mono"), "status"] == "transfers" and c.loc[("THP1", "Classical mono"), "spearman_vs_op3"] > .8, c
        assert c.loc[("MCF7", "Classical mono"), "status"] != "transfers", "random line passed"
        assert c.loc[("U937", "Classical mono"), "status"] == "report_only", "a thin line must not pass"
        assert c.loc[("NOMO1", "Classical mono"), "status"] == "no_signatures"
        assert c.loc[("MCF7", "NK"), "status"] == "too_few_de_genes", "5 DE genes must be refused"
        assert (c.status == "transfers").sum() == 2, c.status
        assert (pd.read_csv(root / "t24/passing_pairs.tsv", sep="\t").cell_line == "THP1").all()
        for name in ("disease_basis", "disease_vectors", "context_agreement", "controls"):
            assert (root / f"t24/{name}.tsv").exists(), name
        # A time filter has to change the answer: pooling with random 6 h signatures weakens the faithful line
        assert run(root, "all").returncode == 0
        pooled = pd.read_csv(root / "all/calibration.tsv", sep="\t").set_index(["cell_line", "population"])
        assert pooled.loc[("THP1", "Classical mono"), "spearman_vs_op3"] < c.loc[("THP1", "Classical mono"), "spearman_vs_op3"]
    print("test_lincs_rescue ok")


if __name__ == "__main__":
    main()
