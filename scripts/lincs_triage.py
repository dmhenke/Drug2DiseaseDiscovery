#!/usr/bin/env python
"""Which LINCS cell lines are worth using as proxies for each PBMC population, and how many compounds do they add?

Reads only LINCS metadata (sig_info, cell_info, gene_info). Answers three questions before any multi-GB
download:
  1. How many distinct compounds (trt_cp) were profiled in each haematopoietic cell line?
  2. How many of the 140 OP3 compounds are in that line? Those are the calibration set: OP3 measured them in
     real PBMCs, so the overlap measures what the cell-line context costs. No overlap means no calibration
     and the expansion is not interpretable.
  3. How many SLE signature genes exist in the L1000 12,328-gene space, split by whether L1000 measured the
     gene directly (landmark) or inferred it from the 978 landmarks. Inferred genes are a model output, not a
     measurement, and are tagged as such downstream.

Expectation, written before the run: the monocytic lines (THP1, HL60, U937, NOMO1, PL21) are the closest
proxies for the myeloid populations and JURKAT for CD4 T; no line is a proxy for NK. Phase I (GSE92742) should
hold more compounds per line than Phase II. If a line has fewer than ~50 OP3 compounds, its calibration
correlation will have an interval too wide to act on and it should be dropped.

Usage: python scripts/lincs_triage.py shared/lincs/meta analysis/lincs_triage [--op3 shared/perturbation/op3]
       python scripts/lincs_triage.py --self-check
"""
import argparse
import gzip
import re
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
# Haematopoietic lines, by lineage, with the PBMC population each is nearest to. Lineage is from the cell line's
# origin (checkable), not from a similarity score: a monocytic leukaemia line is the nearest available thing to a
# monocyte, which is not the same as being one.
PROXY = {
    "myeloid": ("Classical mono, Non-classical mono",
                ["THP1", "HL60", "U937", "NOMO1", "PL21", "SKM1", "MONOMAC6", "OCIAML2", "OCIAML3", "EOL1",
                 "MV411", "KASUMI1", "HEL", "K562", "HAP1"]),
    "T": ("CD4 T", ["JURKAT", "CCRFCEM", "MOLT4", "HUT78", "SUPT1"]),
    "B": ("B cells", ["RAJI", "RAMOS", "DAUDI", "GRANTA519", "MINO", "JEKO1", "SUDHL4", "TOLEDO", "WSUDLCL2",
                      "RPMI8226", "MM1S", "U266", "NCIH929", "KMS11", "L363", "OPM2"]),
    "NK": ("NK", ["NK92", "KHYG1", "YT"]),
}
MIN_OP3 = 50  # below this the calibration interval is too wide to act on; set here, not after seeing results


def read(path):
    with gzip.open(path, "rt", errors="replace") as fh:
        return pd.read_csv(fh, sep="\t", low_memory=False)


def norm(s):
    """Cell ids and compound names differ only in punctuation and case across the two phases and OP3."""
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def op3_compounds(d):
    """OP3 compound names and their LINCS BRD ids, from the DE tables' obs. Returns (names, brd_ids)."""
    import anndata as ad
    obs = pd.concat([ad.read_h5ad(d / f).obs for f in ("de_train.h5ad", "de_test.h5ad")], ignore_index=True)
    names = {norm(x) for x in obs.sm_name.unique()}
    brd = {str(x) for x in obs.get("sm_lincs_id", pd.Series(dtype=str)).unique() if str(x).startswith("BRD")}
    return names, brd


def phase(meta, label, sig_file, op3_names, op3_brd):
    sig = read(meta / sig_file)
    cell_col = "cell_id" if "cell_id" in sig else "cell_iname"
    sig = sig[sig.pert_type == "trt_cp"] if "pert_type" in sig else sig
    sig["cell_key"] = sig[cell_col].map(norm)
    sig["name_key"] = sig.pert_iname.map(norm)
    brd_col = "pert_id" if "pert_id" in sig else None
    rows = []
    for lineage, (pops, lines) in PROXY.items():
        for line in lines:
            s = sig[sig.cell_key == norm(line)]
            if s.empty:
                continue
            overlap = s[s.name_key.isin(op3_names) | (s[brd_col].isin(op3_brd) if brd_col else False)]
            rows.append(dict(phase=label, lineage=lineage, pbmc_populations=pops, cell_line=line,
                             signatures=len(s), compounds=s.name_key.nunique(),
                             op3_compounds_present=overlap.name_key.nunique(),
                             doses=s.pert_dose.nunique() if "pert_dose" in s else pd.NA,
                             times=s.pert_time.nunique() if "pert_time" in s else pd.NA))
    all_lines = sig.groupby("cell_key").pert_iname.nunique().sort_values(ascending=False)
    return pd.DataFrame(rows), all_lines, sig.name_key.nunique()


def genes(meta, sig_path):
    g = read(meta / sig_path)
    sym = "pr_gene_symbol" if "pr_gene_symbol" in g else "gene_symbol"
    typ = "pr_is_lm" if "pr_is_lm" in g else "feature_space"
    return g[[sym, typ]].rename(columns={sym: "gene", typ: "is_landmark"})


def self_check():
    assert norm("THP-1") == norm("THP1") == "thp1"
    assert norm("Mometasone Furoate") == "mometasonefuroate"
    s = pd.DataFrame(dict(pert_type="trt_cp", cell_id=["THP1"] * 3 + ["A375"], pert_iname=["a", "b", "a", "c"],
                          pert_id=["BRD-1", "BRD-2", "BRD-1", "BRD-3"], pert_dose=1.0, pert_time=24))
    assert s[s.cell_id.map(norm) == "thp1"].pert_iname.nunique() == 2


def main():
    p = argparse.ArgumentParser()
    p.add_argument("meta", type=Path, nargs="?"); p.add_argument("out", type=Path, nargs="?")
    p.add_argument("--op3", type=Path, default=ROOT / "shared/perturbation/op3")
    p.add_argument("--self-check", action="store_true")
    a = p.parse_args()
    self_check()
    if a.self_check:
        print("self-check ok")
        return
    assert a.meta and a.out, "need meta and out directories"
    a.out.mkdir(parents=True, exist_ok=True)

    if (a.op3 / "de_train.h5ad").exists():
        op3_names, op3_brd = op3_compounds(a.op3)
        src = f"{len(op3_names)} OP3 compounds, {len(op3_brd)} with a LINCS BRD id"
    else:
        tsv = ROOT / "analysis/op3_stage_p/compounds.tsv"
        op3_names = {norm(x) for x in pd.read_csv(tsv, sep="\t").sm_name}
        op3_brd, src = set(), f"{len(op3_names)} OP3 compound names from {tsv.name} (no BRD ids: OP3 h5ad absent)"

    tabs, lines, totals = [], {}, {}
    for label, sig_f in (("GSE92742", "GSE92742_Broad_LINCS_sig_info.txt.gz"),
                         ("GSE70138", "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz")):
        t, ln, tot = phase(a.meta, label, sig_f, op3_names, op3_brd)
        tabs.append(t); lines[label] = ln; totals[label] = tot
    tab = pd.concat(tabs, ignore_index=True).sort_values(["lineage", "op3_compounds_present"], ascending=[True, False])
    tab.to_csv(a.out / "cell_line_triage.tsv", sep="\t", index=False)

    gi = genes(a.meta, "GSE92742_Broad_LINCS_gene_info.txt.gz")
    sig = pd.read_csv(ROOT / "analysis/sle_signatures/signatures.tsv", sep="\t")
    cov = []
    for pop, s in sig.groupby("population", sort=False):
        m = gi[gi.gene.isin(s.gene)]
        lm = m.is_landmark.astype(str).isin(["1", "landmark", "True"]).sum()
        cov.append(dict(population=pop, signature_genes=len(s), in_l1000=len(m), landmark=int(lm),
                        inferred=len(m) - int(lm), fraction=round(len(m) / len(s), 3)))
    cov = pd.DataFrame(cov)
    cov.to_csv(a.out / "signature_gene_coverage.tsv", sep="\t", index=False)

    lines_tab = pd.concat([v.rename(k) for k, v in lines.items()], axis=1).fillna(0).astype(int)
    lines_tab.to_csv(a.out / "compounds_per_cell_line.tsv", sep="\t")

    keep = tab[tab.op3_compounds_present >= MIN_OP3]
    out = [f"LINCS triage. {src}.",
           f"Distinct compounds per phase: " + ", ".join(f"{k} {v}" for k, v in totals.items()), "",
           "Haematopoietic lines (proxies for the PBMC populations):", tab.to_string(index=False), "",
           "SLE signature gene coverage in the L1000 12,328-gene space:", cov.to_string(index=False), "",
           f"Lines with at least {MIN_OP3} OP3 compounds for calibration: "
           f"{len(keep)} ({', '.join(keep.phase + ':' + keep.cell_line) if len(keep) else 'none'})"]
    for pop, s in cov.set_index("population").iterrows():
        out.append(f"  [{'ok' if s.fraction >= 0.5 else 'FAIL'}] {pop}: {s.fraction:.0%} of signature genes in L1000 "
                   f"({s.landmark} measured, {s.inferred} inferred)")
    for lineage, (pops, _) in PROXY.items():
        k = keep[keep.lineage == lineage]
        out.append(f"  [{'ok' if len(k) else 'FAIL'}] {lineage} -> {pops}: "
                   f"{len(k)} usable line(s); best adds {int(k.compounds.max()) if len(k) else 0} compounds")
    best = keep.groupby("phase").compounds.max()
    out += ["", "Level 5 matrix to download next (gzipped sizes from the GEO FTP listing, 2026-11-06):",
            "  GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz   20G  (~23 GB unzipped)",
            "  GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz  5.0G  (~5.8 GB unzipped)",
            f"  Compounds available in usable lines: {best.to_dict() if len(best) else 'none'}",
            "A FAIL above means the expansion cannot be calibrated for that population; report it, do not pool it."]
    (a.out / "triage.txt").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
