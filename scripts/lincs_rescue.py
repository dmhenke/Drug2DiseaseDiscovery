#!/usr/bin/env python
"""Score every LINCS compound in a PBMC-proxy cell line for how specifically it reverses the SLE signature.

This widens the compound axis with measurements, not predictions: `scripts/op3_nn_rescue_check.py` showed that
predicting an unmeasured compound's profile from structure recovers measured rescue rank at Spearman 0.05-0.33,
so imputed rows would be noise that looks like data. A compound measured in a monocytic leukaemia line is a real
measurement in the wrong cell; the cost of that substitution is measurable here, which is the point.

Per cell line, on the SLE signature genes L1000 covers:
  spec = -((t - m) . d) / (|t - m| |d|)   t = compound z-scores, m = median compound in that line, d = SLE log2FC
Scale-free on purpose: L1000 z-scores and OP3 log2FC are different units, so ||d + t|| is not comparable across
platforms and the match-to-healthy form used for OP3 is not computed here.

CALIBRATION GATE. The OP3 compounds present in the line are scored both ways, and the Spearman between the
LINCS-line score and the OP3 PBMC score is reported with a compound-resampled CI. If that interval includes
zero, this line does not transfer to that PBMC population and its compounds must not be ranked beside OP3's;
the gate is recorded in the output either way.

Expectation, written before the run: calibration near the 0.30 that one OP3 cell type achieves for another is
the best case; a monocytic leukaemia line is further from a classical monocyte than a T cell is, so lower is
expected. Interferon-driven ISGs may respond differently in a line with an altered interferon pathway.

Usage: python scripts/lincs_rescue.py shared/lincs/<file>.gctx shared/lincs/meta analysis/lincs_rescue \
           --lines THP1 HL60 U937 --phase GSE70138
       python scripts/lincs_rescue.py --self-check
"""
import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from lincs_triage import PROXY, norm, read  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
BOOTS, CHUNK = 2000, 2000
SIG_FILES = {"GSE92742": ("GSE92742_Broad_LINCS_sig_info.txt.gz", "GSE92742_Broad_LINCS_gene_info.txt.gz"),
             "GSE70138": ("GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
                          "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")}


def spec(T, m, d):
    """Cosine between the compound-specific part of each profile and the reversed disease vector."""
    R = np.atleast_2d(T) - m
    n = np.linalg.norm(R, axis=1) * np.linalg.norm(d)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, -(R @ d) / n, np.nan)


def gctx_slice(path, col_ids, row_ids):
    """Columns for the wanted signatures, rows for the wanted genes. Orientation is detected, not assumed."""
    with h5py.File(path, "r") as f:
        decode = lambda k: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f[k][:]])  # noqa: E731
        rid, cid, M = decode("0/META/ROW/id"), decode("0/META/COL/id"), f["0/DATA/0/matrix"]
        assert set(M.shape) == {len(rid), len(cid)} or len(rid) == len(cid), \
            f"matrix {M.shape} matches neither {len(rid)} rows nor {len(cid)} cols"
        transposed = M.shape[0] == len(cid) and M.shape[1] == len(rid)  # gctx stores (cols, rows)
        rpos = pd.Series(np.arange(len(rid)), rid).reindex(row_ids).dropna().astype(int)
        cpos = pd.Series(np.arange(len(cid)), cid).reindex(col_ids).dropna().astype(int)
        assert len(rpos) and len(cpos), "no requested gene or signature id is in this matrix"
        gi = np.sort(rpos.to_numpy())
        out = np.empty((len(cpos), len(gi)), dtype=np.float32)
        ci = cpos.to_numpy()
        for s in range(0, len(ci), CHUNK):
            blk = ci[s:s + CHUNK]
            order = np.argsort(blk)
            take = blk[order]
            part = M[take, :][:, gi] if transposed else M[gi, :][:, take].T
            out[np.arange(s, s + len(blk))[order]] = part
        genes_out = pd.Series(rid[gi]).to_list()
        return pd.DataFrame(out, index=cpos.index, columns=genes_out)


def op3_scores():
    """OP3 measured scores in the same scale-free form, per population."""
    g = pd.read_csv(ROOT / "analysis/op3_rescue/rescue_genes.tsv", sep="\t")
    out = {}
    for pop, s in g.groupby("population", sort=False):
        T = s.pivot(index="sm_name", columns="gene", values="drug")
        d = s.drop_duplicates("gene").set_index("gene").disease[T.columns].to_numpy()
        X = T.to_numpy()
        out[pop] = pd.Series(spec(X, np.median(X, 0), d), T.index)
    return out


def interval(x, y, rng):
    n = len(x)
    v = [spearmanr(x[i], y[i]).statistic for i in (rng.integers(0, n, n) for _ in range(BOOTS))]
    return float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))


def self_check(tmp=None):
    d = np.array([2.0, -1.0, 0.5])
    m = np.zeros(3)
    assert np.isclose(spec(-d, m, d), 1) and np.isclose(spec(d, m, d), -1)
    assert np.isnan(spec(m, m, d)).all()
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = Path(tmp or td) / "t.gctx"
        rows, cols = ["7", "8", "9"], ["s1", "s2"]
        M = np.array([[1.0, 2, 3], [4, 5, 6]], dtype=np.float32)  # (cols, rows), as gctx stores it
        with h5py.File(p, "w") as f:
            f.create_dataset("0/DATA/0/matrix", data=M)
            f.create_dataset("0/META/ROW/id", data=np.array(rows, dtype="S"))
            f.create_dataset("0/META/COL/id", data=np.array(cols, dtype="S"))
        got = gctx_slice(p, ["s2", "s1"], ["9", "7"])
        assert list(got.index) == ["s2", "s1"] and got.loc["s1", "9"] == 3 and got.loc["s2", "7"] == 4, got


def main():
    p = argparse.ArgumentParser()
    p.add_argument("gctx", type=Path, nargs="?"); p.add_argument("meta", type=Path, nargs="?")
    p.add_argument("out", type=Path, nargs="?")
    p.add_argument("--phase", default="GSE70138", choices=list(SIG_FILES))
    p.add_argument("--lines", nargs="+", default=["THP1", "HL60", "U937", "JURKAT"])
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--self-check", action="store_true")
    a = p.parse_args()
    self_check()
    if a.self_check:
        print("self-check ok")
        return
    assert a.gctx and a.meta and a.out, "need gctx, meta and out"
    a.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    sig_f, gene_f = SIG_FILES[a.phase]
    sig = read(a.meta / sig_f)
    sig = sig[sig.pert_type == "trt_cp"] if "pert_type" in sig else sig
    cell_col = "cell_id" if "cell_id" in sig else "cell_iname"
    sig["cell_key"] = sig[cell_col].map(norm)
    sid = "sig_id" if "sig_id" in sig else "distil_id"

    gi = read(a.meta / gene_f)
    sym = "pr_gene_symbol" if "pr_gene_symbol" in gi else "gene_symbol"
    gid = "pr_gene_id" if "pr_gene_id" in gi else "gene_id"
    lm = "pr_is_lm" if "pr_is_lm" in gi else "feature_space"
    gi["_lm"] = gi[lm].astype(str).isin(["1", "landmark", "True"])

    sle = pd.read_csv(ROOT / "analysis/sle_signatures/signatures.tsv", sep="\t")
    op3 = op3_scores()
    pop_of_line = {norm(l): pops.split(", ") for pops, lines in PROXY.values() for l in lines}

    rows, calib = [], []
    for line in a.lines:
        s = sig[sig.cell_key == norm(line)]
        if s.empty:
            calib.append(dict(cell_line=line, population=None, note="no trt_cp signatures in this phase"))
            continue
        for pop in pop_of_line.get(norm(line), []):
            sp = sle[sle.population == pop]
            if sp.empty:
                continue
            g = gi[gi[sym].isin(sp.gene)][[sym, gid, "_lm"]].drop_duplicates(sym)
            if len(g) < 0.5 * len(sp):
                calib.append(dict(cell_line=line, population=pop, signature_genes=len(sp), in_l1000=len(g),
                                  note="FAIL under half the signature genes in L1000; not scored"))
                continue
            Z = gctx_slice(a.gctx, s[sid].to_list(), g[gid].astype(str).to_list())
            sym_of = dict(zip(g[gid].astype(str), g[sym]))
            Z.columns = [sym_of[c] for c in Z.columns]
            d = sp.set_index("gene").log2fc_sle_vs_healthy[Z.columns].to_numpy()
            prof = Z.groupby(s.set_index(sid).pert_iname.reindex(Z.index).to_numpy()).median()  # over doses and times
            n_sig = Z.groupby(s.set_index(sid).pert_iname.reindex(Z.index).to_numpy()).size()
            X = prof.to_numpy()
            sc = spec(X, np.median(X, 0), d)
            res = pd.DataFrame(dict(population=pop, cell_line=line, phase=a.phase, sm_name=prof.index,
                                    specific_reversal=sc, signatures=n_sig.to_numpy(),
                                    genes_scored=len(d), genes_measured=int(g._lm.sum()),
                                    genes_inferred=int((~g._lm).sum()), evidence=f"measured_{line}"))
            res["specific_rank"] = res.specific_reversal.rank(ascending=False, method="min")
            rows.append(res)
            # Calibration against the same compounds measured in real PBMCs
            o = op3[pop]
            key = pd.Series(res.specific_reversal.to_numpy(), [norm(x) for x in res.sm_name])
            ok = pd.Series(o.to_numpy(), [norm(x) for x in o.index])
            both = key.index.intersection(ok.index)
            if len(both) >= 10:
                x, y = key[both].to_numpy(), ok[both].to_numpy()
                lo, hi = interval(x, y, rng)
                r = spearmanr(x, y).statistic
                calib.append(dict(cell_line=line, population=pop, op3_overlap=len(both), spearman_vs_op3=round(float(r), 3),
                                  ci_lo=round(lo, 3), ci_hi=round(hi, 3), compounds_added=len(res),
                                  note="transfers" if lo > 0 else "FAIL CI includes zero; do not rank beside OP3"))
            else:
                calib.append(dict(cell_line=line, population=pop, op3_overlap=len(both), compounds_added=len(res),
                                  note="FAIL fewer than 10 OP3 compounds; no calibration"))
    res = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    cal = pd.DataFrame(calib)
    res.to_csv(a.out / "lincs_rescue.tsv", sep="\t", index=False, float_format="%.4f")
    cal.to_csv(a.out / "calibration.tsv", sep="\t", index=False)
    passed = cal[cal.note.astype(str).eq("transfers")] if len(cal) else cal
    lines_out = [f"LINCS {a.phase}: scored {0 if res.empty else res.sm_name.nunique()} distinct compounds "
                 f"across {0 if res.empty else res.cell_line.nunique()} cell line(s); OP3 contributes 140.", "",
                 cal.to_string(index=False), ""]
    for _, r in cal.iterrows():
        lines_out.append(f"  [{'ok' if r.get('note') == 'transfers' else 'FAIL'}] {r.cell_line} -> {r.population}: "
                         f"{r.get('note')}")
    lines_out += ["", f"Usable expansion: {len(passed)} line-population pair(s). Compounds from a FAILed pair stay in "
                      "lincs_rescue.tsv with their evidence tag and must not be pooled into a ranking with OP3.",
                  "L1000 z-scores are not log2FC, so only the rank of specific_reversal is comparable to OP3; the "
                  "match-to-healthy value is not computed for these rows."]
    (a.out / "run_summary.txt").write_text("\n".join(lines_out) + "\n")
    print("\n".join(lines_out))


if __name__ == "__main__":
    main()
