#!/usr/bin/env python
"""Score every LINCS compound in a cell line for how specifically it reverses the SLE signature, and measure
what the cell line costs by comparing the OP3 compounds it shares with real PBMCs.

This widens the compound axis with measurements, not predictions: scripts/op3_nn_rescue_check.py showed that
predicting an unmeasured compound's profile from structure recovers measured rescue rank at Spearman 0.05-0.33,
so imputed rows would be noise that renders like data.

Two things the triage forced (2026-11-06):
  * The disease vector is built from the **L1000 landmark genes of the full DE**, not from the 134-gene SLE
    signature: only 4 of 66 classical-monocyte signature genes are landmarks, so the signature version would
    score genes L1000 infers from the 978 landmarks with a linear model, i.e. score the model. Landmark genes
    with padj<0.05 are primary where at least MIN_SIG of them qualify; `genes_basis` records what was used.
  * Immune lines hold only a few hundred compounds each, so the core lines (far context, ~20k compounds) are
    scored too, and the calibration against OP3 is what decides whether their ranking means anything.

  spec = -((t - m) . d) / (|t - m| |d|)   t = compound z-scores, m = median compound in that line, d = SLE log2FC
Scale-free on purpose: L1000 z-scores and OP3 log2FC are different units, so OP3's match-to-healthy value is not
comparable and is not computed for these rows. Compounds are matched to OP3 on InChIKey connectivity block.

CALIBRATION GATE: Spearman between the line's score and the OP3 PBMC score over shared compounds, with a
compound-resampled CI. CI through zero means the line does not transfer to that population; its compounds keep
their evidence tag and are not ranked beside OP3's. Recorded either way.

Expectation, written before the run: the ~0.30 that one OP3 cell type reaches for another is the ceiling; a
monocytic leukaemia line should be below it, and a melanoma line lower still or at zero. An interferon signature
may transfer especially badly to a line with a disrupted interferon pathway.

Usage: python scripts/lincs_rescue.py <level5.gctx> shared/lincs/meta analysis/lincs_rescue \
           --phase GSE92742 --de-dir studies/lupus_test/results --lines THP1 U937 NOMO1 SKM1 PL21 JURKAT MCF7 PC3
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
from lincs_triage import ALL_POPS, CORE, PROXY, col, norm, read  # noqa: E402
from lupus_use import block  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
BOOTS, CHUNK, MIN_SIG, MIN_OVERLAP = 2000, 2000, 20, 10
SIG_FILES = {"GSE92742": ("GSE92742_Broad_LINCS_sig_info.txt.gz", "GSE92742_Broad_LINCS_gene_info.txt.gz"),
             "GSE70138": ("GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
                          "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")}
# The full DE file per population, as scripts/run_study.sh writes it (see scripts/l1000_reversal.py).
DE_FILES = {"Classical mono": "DE_human_classical_monocyte.csv",
            "Non-classical mono": "DE_human_non-classical_monocyte.csv",
            "NK": "DE_human_natural_killer_cell.csv",
            "CD4 T": "DE_human_CD4-positive,_alpha-beta_T_cell.csv"}


def spec(T, m, d):
    """Cosine between each profile's compound-specific part and the reversed disease vector."""
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
        gi, ci = np.sort(rpos.to_numpy()), cpos.to_numpy()
        out = np.empty((len(cpos), len(gi)), dtype=np.float32)
        for s in range(0, len(ci), CHUNK):
            blk = ci[s:s + CHUNK]
            order = np.argsort(blk)
            part = M[blk[order], :][:, gi] if transposed else M[gi, :][:, blk[order]].T
            out[np.arange(s, s + len(blk))[order]] = part
        return pd.DataFrame(out, index=cpos.index, columns=pd.Series(rid[gi]).to_list())


def disease_vectors(de_dir, landmarks):
    """Per population, the disease direction over L1000 landmark genes, from the full DE table.

    Returns {population: {"basis": label, "d": Series by gene symbol}}. Falls back to the 134-gene signature
    only when the DE table is absent, and says so in the label.
    """
    out = {}
    sig = pd.read_csv(ROOT / "analysis/sle_signatures/signatures.tsv", sep="\t")
    for pop, f in DE_FILES.items():
        p = Path(de_dir) / f if de_dir else None
        if p and p.exists():
            d = pd.read_csv(p).dropna(subset=["log2FoldChange"])
            gcol = col(d, "gene", "symbol", "gene_name")
            d = d[d[gcol].isin(landmarks)].drop_duplicates(gcol).set_index(gcol)
            s = d[d.padj < 0.05].log2FoldChange if "padj" in d else d.log2FoldChange
            if len(s) >= MIN_SIG:
                out[pop] = dict(basis=f"landmark_padj05 ({len(s)} of {len(d)} DE landmarks)", d=s)
            elif len(d):
                out[pop] = dict(basis=f"landmark_all ({len(d)} genes; only {len(s)} reached padj<0.05)",
                                d=d.log2FoldChange)
        if pop not in out:
            s = sig[sig.population == pop].set_index("gene").log2fc_sle_vs_healthy
            n_lm = int(s.index.isin(landmarks).sum())
            out[pop] = dict(d=s, basis=(f"signature_FALLBACK ({n_lm} of {len(s)} signature genes are landmarks; "
                                        "DE table absent, so most genes scored are L1000-inferred)"
                                        if n_lm < 0.5 * len(s) else f"signature ({n_lm} landmark genes)"))
    return out


def op3_scores(vectors):
    """OP3 measured scores on the same gene basis where it overlaps OP3's genes, keyed by InChIKey block."""
    g = pd.read_csv(ROOT / "analysis/op3_rescue/rescue_genes.tsv", sep="\t")
    smi = pd.read_csv(ROOT / "analysis/op3_rescue/rescue.tsv", sep="\t").drop_duplicates("sm_name")
    blk = {n: block(s) for n, s in zip(smi.sm_name, smi.SMILES)}
    out, basis = {}, {}
    for pop, s in g.groupby("population", sort=False):
        T = s.pivot(index="sm_name", columns="gene", values="drug")
        d = vectors[pop]["d"].reindex(T.columns).dropna()
        if len(d) < MIN_OVERLAP:  # OP3 kept 5,317 genes; the landmark basis may barely overlap this signature
            d = s.drop_duplicates("gene").set_index("gene").disease[T.columns]
            basis[pop] = f"op3 scored on its own signature genes ({len(d)}): landmark basis overlapped too few"
        else:
            basis[pop] = f"op3 scored on {len(d)} shared basis genes"
        X = T[d.index].to_numpy()
        sc = pd.Series(spec(X, np.median(X, 0), d.to_numpy()), T.index)
        o = pd.Series(sc.to_numpy(), [blk.get(n) for n in sc.index]).dropna()
        out[pop] = o[~o.index.duplicated()]
    return out, basis


def interval(x, y, rng):
    n = len(x)
    v = [spearmanr(x[i], y[i]).statistic for i in (rng.integers(0, n, n) for _ in range(BOOTS))]
    return float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))


def context_agreement(res, rng):
    """Do an immune line and a core line rank the same compounds the same way?

    The immune lines share only 20-49 compounds with OP3 (measured 2026-11-06), too few to calibrate, but they
    share hundreds with the core lines. If a monocytic line and a melanoma line agree on those, the cell context
    costs little for this signature and the core expansion is on firmer ground than its OP3 interval alone shows;
    if they disagree, the context matters and a passing OP3 calibration on a core line deserves suspicion.
    This is an internal control between proxies, not evidence about real monocytes.
    """
    out = []
    if res.empty:
        return pd.DataFrame()
    for pop, g in res.groupby("population", sort=False):
        W = g.pivot_table(index="sm_name", columns="cell_line", values="specific_reversal", aggfunc="first")
        tier_of = g.drop_duplicates("cell_line").set_index("cell_line").tier
        for imm in [c for c in W if tier_of.get(c) == "immune"]:
            for core in [c for c in W if tier_of.get(c) == "core"]:
                k = W[[imm, core]].dropna()
                if len(k) >= 30:
                    x, y = k[imm].to_numpy(), k[core].to_numpy()
                    lo, hi = interval(x, y, rng)
                    out.append(dict(population=pop, immune_line=imm, core_line=core, shared_compounds=len(k),
                                    spearman=round(float(spearmanr(x, y).statistic), 3),
                                    ci_lo=round(lo, 3), ci_hi=round(hi, 3),
                                    reading="context costs little" if lo > 0.5 else
                                            "partial agreement" if lo > 0 else "lines disagree; context matters"))
    return pd.DataFrame(out)


def pert_blocks(meta):
    out = {}
    for f in ("GSE92742_Broad_LINCS_pert_info.txt.gz", "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz"):
        if not (meta / f).exists():
            continue
        pi = read(meta / f)
        smi, key = col(pi, "canonical_smiles", "smiles"), col(pi, "inchi_key_prefix", "inchi_key")
        for pid, s, k in zip(pi.pert_id, pi[smi] if smi else [None] * len(pi), pi[key] if key else [None] * len(pi)):
            b = block(s) if s is not None and str(s) not in ("nan", "-666", "restricted") else None
            if b is None and k is not None and str(k) != "nan":
                b = str(k)[:14]
            if b:
                out.setdefault(str(pid), b)
    return out


def self_check():
    d = np.array([2.0, -1.0, 0.5])
    m = np.zeros(3)
    assert np.isclose(spec(-d, m, d), 1) and np.isclose(spec(d, m, d), -1) and np.isnan(spec(m, m, d)).all()
    import tempfile
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "t.gctx"
        M = np.array([[1.0, 2, 3], [4, 5, 6]], dtype=np.float32)  # (cols, rows), as gctx stores it
        with h5py.File(p, "w") as f:
            f.create_dataset("0/DATA/0/matrix", data=M)
            f.create_dataset("0/META/ROW/id", data=np.array(["7", "8", "9"], dtype="S"))
            f.create_dataset("0/META/COL/id", data=np.array(["s1", "s2"], dtype="S"))
        got = gctx_slice(p, ["s2", "s1"], ["9", "7"])
        assert list(got.index) == ["s2", "s1"] and got.loc["s1", "9"] == 3 and got.loc["s2", "7"] == 4, got


def main():
    p = argparse.ArgumentParser()
    p.add_argument("gctx", type=Path, nargs="?"); p.add_argument("meta", type=Path, nargs="?")
    p.add_argument("out", type=Path, nargs="?")
    p.add_argument("--phase", default="GSE92742", choices=list(SIG_FILES))
    p.add_argument("--de-dir", default=None, help="study results dir holding DE_human_*.csv (strongly preferred)")
    p.add_argument("--lines", nargs="+", default=["THP1", "U937", "NOMO1", "SKM1", "PL21", "JURKAT", "WSUDLCL2",
                                                  "MCF7", "PC3", "A375", "HA1E"])
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
    if "pert_type" in sig:
        sig = sig[sig.pert_type == "trt_cp"]
    cell_c, sid = col(sig, "cell_id", "cell_iname"), col(sig, "sig_id", "distil_id")
    sig = sig.assign(cell_key=sig[cell_c].map(norm))

    gi = read(a.meta / gene_f)
    sym, gid = col(gi, "pr_gene_symbol", "gene_symbol"), col(gi, "pr_gene_id", "gene_id")
    lm = col(gi, "pr_is_lm", "feature_space")
    gi = gi.assign(is_lm=gi[lm].astype(str).isin(["1", "landmark", "True"]))
    landmarks = set(gi[gi.is_lm][sym])
    id_of = dict(zip(gi[sym], gi[gid].astype(str)))

    vectors = disease_vectors(a.de_dir, landmarks)
    op3, op3_basis = op3_scores(vectors)
    blocks = pert_blocks(a.meta)
    pops_of = {norm(l): pops.split(", ") for pops, lines in PROXY.values() for l in lines}
    pops_of.update({norm(l): ALL_POPS.split(", ") for l in CORE})

    rows, calib = [], []
    for line in a.lines:
        s = sig[sig.cell_key == norm(line)]
        if s.empty:
            calib.append(dict(cell_line=line, population=None, note=f"FAIL no trt_cp signatures in {a.phase}"))
            continue
        tier = "immune" if any(norm(line) == norm(x) for _, ls in PROXY.values() for x in ls) else "core"
        for pop in pops_of.get(norm(line), []):
            v = vectors[pop]
            d = v["d"][[g for g in v["d"].index if g in id_of]]
            if len(d) < MIN_SIG:
                calib.append(dict(cell_line=line, population=pop, tier=tier, genes_basis=v["basis"],
                                  note=f"FAIL only {len(d)} basis genes are in L1000; not scored"))
                continue
            Z = gctx_slice(a.gctx, s[sid].to_list(), [id_of[g] for g in d.index])
            sym_of = {id_of[g]: g for g in d.index}
            Z.columns = [sym_of[c] for c in Z.columns]
            d = d[Z.columns]
            name = s.set_index(sid).pert_iname.reindex(Z.index)
            pid = s.set_index(sid).pert_id.reindex(Z.index) if "pert_id" in s else name
            prof = Z.groupby(name.to_numpy()).median()          # over doses and times
            n_sig = Z.groupby(name.to_numpy()).size()
            blk = pd.Series(pid.to_numpy(), name.to_numpy()).groupby(level=0).first().map(blocks)
            X = prof.to_numpy()
            res = pd.DataFrame(dict(population=pop, cell_line=line, tier=tier, phase=a.phase, sm_name=prof.index,
                                    inchikey_block=blk.reindex(prof.index).to_numpy(),
                                    specific_reversal=spec(X, np.median(X, 0), d.to_numpy()),
                                    signatures=n_sig.to_numpy(), genes_scored=len(d), genes_basis=v["basis"],
                                    evidence=f"measured_{line}"))
            res["specific_rank"] = res.specific_reversal.rank(ascending=False, method="min")
            rows.append(res)
            o = op3[pop]
            k = res.dropna(subset=["inchikey_block"]).drop_duplicates("inchikey_block").set_index("inchikey_block")
            both = k.index.intersection(o.index)
            if len(both) >= MIN_OVERLAP:
                x, y = k.specific_reversal[both].to_numpy(), o[both].to_numpy()
                lo, hi = interval(x, y, rng)
                calib.append(dict(cell_line=line, population=pop, tier=tier, op3_overlap=len(both),
                                  spearman_vs_op3=round(float(spearmanr(x, y).statistic), 3),
                                  ci_lo=round(lo, 3), ci_hi=round(hi, 3), compounds_added=len(res),
                                  genes_basis=v["basis"],
                                  note="transfers" if lo > 0 else "FAIL CI includes zero; do not rank beside OP3"))
            else:
                calib.append(dict(cell_line=line, population=pop, tier=tier, op3_overlap=len(both),
                                  compounds_added=len(res), genes_basis=v["basis"],
                                  note=f"FAIL fewer than {MIN_OVERLAP} shared OP3 compounds; no calibration"))
    res = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    cal = pd.DataFrame(calib)
    res.to_csv(a.out / "lincs_rescue.tsv", sep="\t", index=False, float_format="%.4f")
    cal.to_csv(a.out / "calibration.tsv", sep="\t", index=False)
    ctx = context_agreement(res, rng)
    ctx.to_csv(a.out / "context_agreement.tsv", sep="\t", index=False)
    passed = cal[cal.note.astype(str).eq("transfers")] if len(cal) else cal
    out = [f"LINCS {a.phase}: scored {0 if res.empty else res.sm_name.nunique()} distinct compounds across "
           f"{0 if res.empty else res.cell_line.nunique()} cell line(s); OP3 contributes 140 in real PBMCs.",
           "Disease-vector basis: " + "; ".join(f"{k}: {v['basis']}" for k, v in vectors.items()),
           "OP3 side: " + "; ".join(f"{k}: {v}" for k, v in op3_basis.items()), "", cal.to_string(index=False), ""]
    for _, r in cal.iterrows():
        out.append(f"  [{'ok' if r.get('note') == 'transfers' else 'FAIL'}] {r.get('tier')} {r.cell_line} -> "
                   f"{r.population}: {r.get('note')}")
    if len(ctx):
        out += ["Immune line vs core line on the compounds they share (internal control, not evidence about real "
                "monocytes):", ctx.to_string(index=False), ""]
    usable = 0 if res.empty or passed.empty else res.merge(
        passed[["cell_line", "population"]], on=["cell_line", "population"]).sm_name.nunique()
    out += ["", f"Usable expansion: {len(passed)} line-population pair(s), {usable} distinct compounds whose ranking "
                "is calibrated against real PBMCs. Compounds from a FAILed pair stay in lincs_rescue.tsv with their "
                "evidence tag and must not be pooled into a ranking with OP3.",
            "L1000 z-scores are not log2FC, so only the rank of specific_reversal is comparable to OP3; the "
            "match-to-healthy value is not computed for these rows.",
            "A 'signature_FALLBACK' basis means the DE table was missing and the score rests on genes L1000 "
            "inferred from its 978 landmarks; rerun with --de-dir before using those numbers."]
    (a.out / "run_summary.txt").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
