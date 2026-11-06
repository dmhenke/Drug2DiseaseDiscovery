#!/usr/bin/env python
"""Can LINCS compounds measured in cell lines be ranked for SLE reversal in PBMC populations?

For each line and population, on one gene set: L1000 landmarks (measured, not inferred) that OP3 also holds and
that reach padj < 0.05 in the full SLE DE, ribosomal and mitochondrial genes removed.
  spec = -((t - m) . d) / (|t - m| |d|)   t = compound profile, m = median compound in that context, d = SLE log2FC
Scale-free, because L1000 z-scores and OP3 log2FC are different units. OP3 is rescored here from its h5ad on the
same genes, so the two sides of the calibration are like for like.

GATE: one-sided Spearman of line score vs OP3 score over shared compounds (InChIKey block), BH across every pair
with at least MIN_ACT compounds; q < 0.05 passes. 20 to 49 shared compounds is report-only, below 20 is not tested.
A pass means this line orders those shared compounds like PBMCs do. It does not validate any one unshared compound.
No DE table means no score: there is no fallback to the signature genes, which L1000 mostly infers.

Expectation, written before the first valid run: one OP3 cell type predicts another at about 0.30, so that is a
ceiling. Monocytic lines should sit below it and non-immune lines at or near zero. Passes confined to one population
or a negative correlation would point to artefact. Known mismatches: OP3 is 1 uM, most L1000 signatures are higher
and doses are pooled by median; the scoring assumes disease and drug effects add.

Usage: python scripts/lincs_rescue.py <level5.gctx> shared/lincs/meta outdir --de studies/<study>/results/DE_all_human.csv \
           --op3 shared/perturbation/op3 --phase GSE92742 [--time 24] [--lines ...] [--seed 0]
Test:  python scripts/test_lincs_rescue.py
"""
import argparse
import re
import socket
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
from compare_cohorts import POPULATIONS  # noqa: E402
from lincs_triage import CORE, PROXY, col, norm, pert_blocks, read  # noqa: E402
from lupus_use import block  # noqa: E402
from op3_rescue import bh, load  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
BOOTS, CHUNK, MIN_SIG, MIN_REPORT, MIN_ACT = 2000, 2000, 20, 20, 50  # MIN_REPORT, MIN_ACT as in lincs_triage
RIBO = re.compile(r"^(RPL|RPS|MT-)")
CONTROLS = re.compile("ruxolitinib|filgotinib|tofacitinib|baricitinib|prednisolone|dexamethasone|"
                      "hydroxychloroquine|methotrexate|mycophenol", re.I)
SIG_FILES = {"GSE92742": ("GSE92742_Broad_LINCS_sig_info.txt.gz", "GSE92742_Broad_LINCS_gene_info.txt.gz"),
             "GSE70138": ("GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz",
                          "GSE70138_Broad_LINCS_gene_info_2017-03-06.txt.gz")}
PERT_FILES = ["GSE92742_Broad_LINCS_pert_info.txt.gz", "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz"]
DEFAULT_LINES = ["THP1", "U937", "NOMO1", "SKM1", "PL21", "JURKAT", "WSUDLCL2", "MCF7", "PC3", "A375", "A549",
                 "HA1E", "HT29", "VCAP", "HCC515", "HEPG2", "NPC"]


def spec(T, m, d):
    """Cosine between each profile's compound-specific part and the reversed disease vector."""
    R = np.atleast_2d(T) - m
    n = np.linalg.norm(R, axis=1) * np.linalg.norm(d)
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.where(n > 0, -(R @ d) / n, np.nan)


def gctx_slice(path, col_ids, row_ids):
    """Signatures by genes. A gctx stores (cols, rows); orientation is checked against the ids."""
    with h5py.File(path, "r") as f:
        decode = lambda k: np.array([x.decode() if isinstance(x, bytes) else str(x) for x in f[k][:]])  # noqa: E731
        rid, cid, M = decode("0/META/ROW/id"), decode("0/META/COL/id"), f["0/DATA/0/matrix"]
        assert M.shape == (len(cid), len(rid)), f"matrix {M.shape} is not ({len(cid)} cols, {len(rid)} rows)"
        rpos = pd.Series(np.arange(len(rid)), rid).reindex(row_ids).dropna().astype(int)
        cpos = pd.Series(np.arange(len(cid)), cid).reindex(col_ids).dropna().astype(int)
        assert len(rpos) and len(cpos), "no requested gene or signature id is in this matrix"
        gi, ci = np.sort(rpos.to_numpy()), cpos.to_numpy()
        out = np.empty((len(ci), len(gi)), dtype=np.float32)
        for s in range(0, len(ci), CHUNK):
            blk = ci[s:s + CHUNK]
            order = np.argsort(blk)
            out[np.arange(s, s + len(blk))[order]] = M[blk[order], :][:, gi]
        return pd.DataFrame(out, index=cpos.index, columns=rid[gi])


def disease_vectors(path, allowed):
    """{population: Series of log2FC by gene} from the full DE, and one review row per population."""
    de = pd.read_csv(path).dropna(subset=["log2FoldChange", "padj"])
    out, rows = {}, []
    for cell, (pop, op3_cell) in POPULATIONS.items():
        g = de[(de.cell_type == cell) & de.gene.isin(allowed) & ~de.gene.str.match(RIBO)].drop_duplicates("gene")
        s = g[g.padj < 0.05].set_index("gene").log2FoldChange
        rows.append(dict(population=pop, op3_cell_type=op3_cell, landmark_genes_tested=len(g), genes_padj05=len(s)))
        if len(s) >= MIN_SIG:
            out[pop] = s
    return out, pd.DataFrame(rows)


def op3_scores(obs, Y, genes, vectors):
    """OP3 measured scores per population on that population's gene set, keyed by InChIKey block."""
    col_of = {g: i for i, g in enumerate(genes)}
    blk = {s: block(s) for s in obs.SMILES.unique()}
    op3_cell = dict(POPULATIONS.values())
    out = {}
    for pop, d in vectors.items():
        rows = np.flatnonzero(obs.cell_type == op3_cell[pop])
        T = Y[rows][:, d.index.map(col_of).to_numpy()]
        sc = spec(T, np.median(T, 0), d.to_numpy())
        s = pd.Series(sc, obs.SMILES.to_numpy()[rows].astype(object)).rename(index=blk).dropna()
        out[pop] = s.groupby(level=0).mean()
    return out


def interval(x, y, rng):
    n = len(x)
    v = [spearmanr(x[i], y[i]).statistic for i in (rng.integers(0, n, n) for _ in range(BOOTS))]
    return float(np.nanpercentile(v, 2.5)), float(np.nanpercentile(v, 97.5))


def context_agreement(res, rng):
    """Do an immune line and a non-immune line rank shared compounds alike? Compares two proxies, not real cells."""
    out = []
    if res.empty:
        return pd.DataFrame()
    for pop, g in res.groupby("population", sort=False):
        W = g.pivot_table(index="sm_name", columns="cell_line", values="specific_reversal", aggfunc="first")
        tier_of = g.drop_duplicates("cell_line").set_index("cell_line").tier
        for imm in [c for c in W if tier_of[c] == "immune"]:
            for core in [c for c in W if tier_of[c] == "core"]:
                k = W[[imm, core]].dropna()
                if len(k) >= 30:
                    lo, hi = interval(k[imm].to_numpy(), k[core].to_numpy(), rng)
                    out.append(dict(population=pop, immune_line=imm, core_line=core, shared_compounds=len(k),
                                    spearman=round(float(spearmanr(k[imm], k[core]).statistic), 3),
                                    ci_lo=round(lo, 3), ci_hi=round(hi, 3)))
    return pd.DataFrame(out)


def gate(cal):
    """BH over pairs that can pass. Status comes from q and the sign; the interval is reported, not used."""
    cal = cal.assign(q_bh=np.nan, status=cal.get("status", "not_tested"))
    act = cal.index[cal.n_op3 >= MIN_ACT]
    if len(act):
        cal.loc[act, "q_bh"] = bh(cal.loc[act, "p_one_sided"].to_numpy())
        cal.loc[act, "status"] = np.where(cal.loc[act, "q_bh"] < 0.05, "transfers",
                                          np.where(cal.loc[act, "ci_hi"] < 0, "inverted", "no_transfer"))
    cal.loc[(cal.n_op3 >= MIN_REPORT) & (cal.n_op3 < MIN_ACT), "status"] = "report_only"
    return cal


def main():
    p = argparse.ArgumentParser()
    p.add_argument("gctx", type=Path); p.add_argument("meta", type=Path); p.add_argument("out", type=Path)
    p.add_argument("--de", type=Path, required=True, help="DE_all_human.csv computed on this machine")
    p.add_argument("--op3", type=Path, required=True)
    p.add_argument("--phase", default="GSE92742", choices=list(SIG_FILES))
    p.add_argument("--lines", nargs="+", default=DEFAULT_LINES)
    p.add_argument("--time", type=int, default=None, help="keep only signatures at this many hours")
    p.add_argument("--seed", type=int, default=0)
    a = p.parse_args()
    for path in (a.gctx, a.de, a.op3 / "de_train.h5ad"):
        if not path.exists():
            sys.exit(f"ERROR: {path} is missing; nothing was scored")
    a.out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(a.seed)

    sig_f, gene_f = SIG_FILES[a.phase]
    sig = read(a.meta / sig_f)
    sig = sig[sig.pert_type == "trt_cp"].assign(cell_key=lambda x: x[col(x, "cell_id", "cell_iname")].map(norm))
    if a.time is not None:
        sig = sig[sig[col(sig, "pert_time")] == a.time]
    gi = read(a.meta / gene_f)
    is_lm = gi[col(gi, "pr_is_lm", "feature_space")].astype(str).isin(["1", "landmark", "True"])
    sym, gid = col(gi, "pr_gene_symbol", "gene_symbol"), col(gi, "pr_gene_id", "gene_id")
    id_of = dict(zip(gi[sym][is_lm], gi[gid][is_lm].astype(str)))

    obs, Y, op3_genes = load(a.op3)
    vectors, basis = disease_vectors(a.de, set(id_of) & set(op3_genes))
    basis.to_csv(a.out / "disease_basis.tsv", sep="\t", index=False)
    if not vectors:
        sys.exit(f"ERROR: no population has {MIN_SIG} landmark genes at padj<0.05 in {a.de}\n{basis.to_string(index=False)}")
    pd.concat([d.rename("log2fc").reset_index().assign(population=pop) for pop, d in vectors.items()]).to_csv(
        a.out / "disease_vectors.tsv", sep="\t", index=False)
    op3 = op3_scores(obs, Y, op3_genes, vectors)
    blocks = pert_blocks(a.meta, PERT_FILES)
    union = sorted({g for d in vectors.values() for g in d.index})
    pops_of = {norm(x): pops.split(", ") for pops, xs in PROXY.values() for x in xs}
    pops_of.update({norm(x): [pop for pop, _ in POPULATIONS.values()] for x in CORE})
    sym_of = {id_of[g]: g for g in union}

    rows, calib = [], []
    for line in a.lines:
        s = sig[sig.cell_key == norm(line)]
        pops = pops_of.get(norm(line), [])
        tier = "core" if norm(line) in {norm(x) for x in CORE} else "immune"
        if s.empty:
            calib += [dict(cell_line=line, population=q, tier=tier, n_op3=0, status="no_signatures") for q in pops]
            continue
        sid = col(s, "sig_id")
        Z = gctx_slice(a.gctx, s[sid].to_list(), [id_of[g] for g in union]).rename(columns=sym_of)
        meta = s.set_index(sid).reindex(Z.index)
        prof, n_sig = Z.groupby(meta.pert_iname.to_numpy()).median(), Z.groupby(meta.pert_iname.to_numpy()).size()
        blk = meta.groupby("pert_iname").pert_id.first().astype(str).map(blocks).reindex(prof.index)
        dose = float(meta[col(meta, "pert_dose")].median()) if col(meta, "pert_dose") else np.nan
        for pop in pops:
            if pop not in vectors:
                calib.append(dict(cell_line=line, population=pop, tier=tier, n_op3=0, status="too_few_de_genes"))
                continue
            d = vectors[pop]
            X = prof[d.index].to_numpy()
            res = pd.DataFrame(dict(population=pop, cell_line=line, tier=tier, phase=a.phase, sm_name=prof.index,
                                    inchikey_block=blk.to_numpy(), specific_reversal=spec(X, np.median(X, 0), d.to_numpy()),
                                    signatures=n_sig.to_numpy(), genes_scored=len(d), evidence=f"measured_{line}"))
            res["specific_rank"] = res.specific_reversal.rank(ascending=False, method="min")
            res["compounds_in_context"] = len(res)
            rows.append(res)
            k = res.dropna(subset=["inchikey_block", "specific_reversal"]).groupby("inchikey_block").specific_reversal.mean()
            both = k.index.intersection(op3[pop].index)
            row = dict(cell_line=line, population=pop, tier=tier, n_op3=len(both), compounds_scored=len(res),
                       genes=len(d), median_dose=dose, signatures=len(s))
            if len(both) >= MIN_REPORT:
                x, y = k[both].to_numpy(), op3[pop][both].to_numpy()
                row["spearman_vs_op3"] = round(float(spearmanr(x, y).statistic), 3)
                row["ci_lo"], row["ci_hi"] = (round(v, 3) for v in interval(x, y, rng))
                row["p_one_sided"] = float(spearmanr(x, y, alternative="greater").pvalue)
            calib.append(row)
    res = pd.concat(rows, ignore_index=True) if rows else pd.DataFrame()
    cal = gate(pd.DataFrame(calib))
    ctx = context_agreement(res, rng)
    passed = cal[cal.status == "transfers"]
    keep = res.merge(passed[["cell_line", "population"]], on=["cell_line", "population"]) if len(passed) and len(res) else res.iloc[:0]
    ctrl = res[res.sm_name.str.contains(CONTROLS)] if len(res) else res
    ctrl = ctrl.assign(percentile=(100 * ctrl.specific_rank / ctrl.compounds_in_context).round(1))
    res.to_csv(a.out / "lincs_rescue.tsv", sep="\t", index=False, float_format="%.4f")
    keep.to_csv(a.out / "passing_pairs.tsv", sep="\t", index=False, float_format="%.4f")
    cal.to_csv(a.out / "calibration.tsv", sep="\t", index=False)
    ctx.to_csv(a.out / "context_agreement.tsv", sep="\t", index=False)
    ctrl.drop(columns=["inchikey_block", "evidence", "phase"]).to_csv(a.out / "controls.tsv", sep="\t", index=False, float_format="%.4f")

    git = lambda *c: subprocess.run(["git", "-C", str(ROOT), *c], capture_output=True, text=True).stdout.strip()  # noqa: E731
    head = (f"REAL RUN {a.phase}, time {a.time or 'all'} h, DE {a.de}; commit {git('rev-parse', '--short', 'HEAD')}; "
            f"{len(git('status', '--porcelain').splitlines())} changed files; host {socket.gethostname()}; seed {a.seed}")
    show = ["cell_line", "population", "tier", "n_op3", "compounds_scored", "genes", "median_dose", "spearman_vs_op3",
            "ci_lo", "ci_hi", "q_bh", "status"]
    tested = int((cal.n_op3 >= MIN_ACT).sum())
    out = [head, "", "Disease basis (landmarks shared with OP3, padj<0.05 needs at least " f"{MIN_SIG}):", basis.to_string(index=False), "",
           cal.reindex(columns=show).round(4).to_string(index=False), "",
           f"BH family: {tested} pairs with at least {MIN_ACT} shared compounds. Status: " + ", ".join(
               f"{k} {v}" for k, v in cal.status.value_counts().items()),
           f"Pairs passing: {len(passed)}; compounds in those lines: {keep.sm_name.nunique() if len(keep) else 0}. "
           "These are candidates in a line that passed, not validated PBMC measurements.",
           "Controls are recorded, not a gate, in controls.tsv. Only ranks are comparable with OP3.", head]
    (a.out / "run_summary.txt").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
