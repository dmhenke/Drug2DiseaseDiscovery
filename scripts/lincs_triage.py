#!/usr/bin/env python
"""Which LINCS cell lines can be used for each PBMC population, and how many compounds does each one add?

Reads only LINCS metadata (sig_info, pert_info, gene_info). Answers, before any multi-GB download:
  1. How many distinct compounds (trt_cp) each cell line holds, in two tiers:
       immune  haematopoietic lines, the nearest available context to a PBMC population. Measured 2026-11-06:
               a few hundred compounds each, not thousands.
       core    the lines LINCS profiled densely (A375, A549, ... VCAP). Far context, roughly 20k compounds.
               Admissible only with a measured transfer calibration, never on the assumption that a drug acts
               the same in a melanoma line as in a monocyte.
  2. How many of the 140 OP3 compounds each line holds. Matching is by InChIKey connectivity block of the
     largest fragment, as in lupus_use.py: OP3's sm_lincs_id values do not start with BRD, so id matching found
     nothing, and name matching misses aliases (GLPG0634 / filgotinib). Those compounds are the calibration set.
  3. How many SLE signature genes L1000 measured directly (landmark, 978 of 12,328) rather than inferred from
     the landmarks. An inferred gene is a linear model's output; scoring an interferon signature on inferred
     genes scores the model. Measured 2026-11-06: 4 of 66 classical-monocyte signature genes are landmarks,
     which is why scripts/lincs_rescue.py builds its disease vector from the landmarks of the full DE instead.

Expectation, written before the first run: monocytic lines are the closest proxies for the myeloid populations
and JURKAT for CD4 T; no line is a proxy for NK. Confirmed, with the unexpected part being how thin the immune
lines are.

Usage: python scripts/lincs_triage.py shared/lincs/meta analysis/lincs_triage [--op3 shared/perturbation/op3]
       python scripts/lincs_triage.py --self-check
"""
import argparse
import gzip
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from lupus_use import block  # noqa: E402  InChIKey connectivity block of the largest fragment

ROOT = Path(__file__).resolve().parent.parent
# Haematopoietic lines by cell origin, with the PBMC population each is nearest to. Lineage is checkable;
# a similarity score is not. A monocytic leukaemia line is the nearest available thing to a monocyte, not one.
PROXY = {
    "myeloid": ("Classical mono, Non-classical mono",
                ["THP1", "HL60", "U937", "NOMO1", "PL21", "SKM1", "MONOMAC6", "OCIAML2", "OCIAML3", "EOL1",
                 "MV411", "KASUMI1", "HEL", "K562", "HAP1"]),
    "T": ("CD4 T, CD8 T", ["JURKAT", "CCRFCEM", "MOLT4", "HUT78", "SUPT1"]),
    "B": ("B", ["RAJI", "RAMOS", "DAUDI", "GRANTA519", "MINO", "JEKO1", "SUDHL4", "TOLEDO", "WSUDLCL2",
                      "RPMI8226", "MM1S", "U266", "NCIH929", "KMS11", "L363", "OPM2"]),
    "NK": ("NK", ["NK92", "KHYG1", "YT"]),
}
CORE = ["A375", "A549", "HA1E", "HCC515", "HEPG2", "HT29", "MCF7", "PC3", "VCAP", "HEK293T", "NPC", "MCF10A"]
ALL_POPS = "Classical mono, Non-classical mono, NK, CD4 T, CD8 T, B"
# Pre-registered bands, fixed before any count was seen: at n = 50 a true Spearman near 0.3 is distinguishable
# from zero, at n = 20 only a strong one is, below that nothing is.
MIN_OP3_ACT, MIN_OP3_REPORT = 50, 20


def read(path):
    with gzip.open(path, "rt", errors="replace") as fh:
        return pd.read_csv(fh, sep="\t", low_memory=False)


def norm(s):
    """Cell ids and compound names differ only in punctuation and case across phases and OP3."""
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def col(df, *names):
    for n in names:
        if n in df:
            return n
    return None


def pert_blocks(meta, files):
    """pert_id -> InChIKey block, from the SMILES where there is one, so salts resolve as on the OP3 side."""
    out = {}
    for f in files:
        if not (meta / f).exists():
            continue
        p = read(meta / f)
        smi, key = col(p, "canonical_smiles", "smiles"), col(p, "inchi_key_prefix", "inchi_key")
        for pid, s, k in zip(p.pert_id, p[smi] if smi else [None] * len(p), p[key] if key else [None] * len(p)):
            b = block(s) if s is not None and str(s) not in ("nan", "-666", "restricted") else None
            if b is None and k is not None and str(k) != "nan":
                b = str(k)[:14]
            if b:
                out[str(pid)] = b
    return out


def op3_compounds(d):
    """OP3 compound names and InChIKey blocks. The h5ad ships SMILES, which is what the blocks come from."""
    import anndata as ad
    obs = pd.concat([ad.read_h5ad(d / f).obs for f in ("de_train.h5ad", "de_test.h5ad")], ignore_index=True)
    comp = obs.drop_duplicates("sm_name")
    return ({norm(x) for x in comp.sm_name}, {b for b in (block(s) for s in comp.SMILES) if b})


def phase(meta, label, sig_file, op3_names, op3_blocks, blocks):
    sig = read(meta / sig_file)
    if "pert_type" in sig:
        sig = sig[sig.pert_type == "trt_cp"]
    cell_c = col(sig, "cell_id", "cell_iname")
    dose_c, time_c = col(sig, "pert_dose", "pert_idose"), col(sig, "pert_time", "pert_itime")
    sig = sig.assign(cell_key=sig[cell_c].map(norm), name_key=sig.pert_iname.map(norm),
                     blk=sig.pert_id.map(blocks) if "pert_id" in sig else None)
    tiers = [("immune", lineage, pops, line) for lineage, (pops, lines) in PROXY.items() for line in lines]
    tiers += [("core", "non-immune", ALL_POPS, line) for line in CORE]
    rows = []
    for tier, lineage, pops, line in tiers:
        s = sig[sig.cell_key == norm(line)]
        if s.empty:
            continue
        by_block = s[s.blk.isin(op3_blocks)] if s.blk.notna().any() else s.iloc[:0]
        by_name = s[s.name_key.isin(op3_names)]
        rows.append(dict(phase=label, tier=tier, lineage=lineage, pbmc_populations=pops, cell_line=line,
                         signatures=len(s), compounds=s.name_key.nunique(),
                         op3_by_inchikey=by_block.name_key.nunique(), op3_by_name=by_name.name_key.nunique(),
                         op3_compounds_present=pd.concat([by_block, by_name]).name_key.nunique(),
                         doses=s[dose_c].nunique() if dose_c else pd.NA,
                         times=s[time_c].nunique() if time_c else pd.NA))
    lines = sig.groupby("cell_key").pert_iname.nunique().sort_values(ascending=False)
    return pd.DataFrame(rows), lines, sig.name_key.nunique()


def genes(meta, f):
    g = read(meta / f)
    sym, lm = col(g, "pr_gene_symbol", "gene_symbol"), col(g, "pr_is_lm", "feature_space")
    return pd.DataFrame(dict(gene=g[sym], is_landmark=g[lm].astype(str).isin(["1", "landmark", "True"])))


def self_check():
    assert norm("THP-1") == norm("THP1") == "thp1"
    assert block("CCS(=O)(=O)N1CC(C1)(CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3") is not None
    # a salt resolves to its parent's block, which is why matching is not done on names or vendor ids
    assert block("CC(=O)OC1=CC=CC=C1C(=O)O") == block("CC(=O)OC1=CC=CC=C1C(=O)[O-].[Na+]")


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
        op3_names, op3_blocks = op3_compounds(a.op3)
        src = f"{len(op3_names)} OP3 compounds, {len(op3_blocks)} with an InChIKey block"
    else:
        c = pd.read_csv(ROOT / "analysis/op3_stage_p/compounds.tsv", sep="\t")
        op3_names = {norm(x) for x in c.sm_name}
        op3_blocks = {b for b in (block(s) for s in c.SMILES) if b}
        src = f"{len(op3_names)} OP3 compounds from compounds.tsv, {len(op3_blocks)} with an InChIKey block"

    blocks = pert_blocks(a.meta, ["GSE92742_Broad_LINCS_pert_info.txt.gz",
                                  "GSE70138_Broad_LINCS_pert_info_2017-03-06.txt.gz",
                                  "GSE70138_Broad_LINCS_pert_info.txt.gz"])
    tabs, lines, totals = [], {}, {}
    for label, sig_f in (("GSE92742", "GSE92742_Broad_LINCS_sig_info.txt.gz"),
                         ("GSE70138", "GSE70138_Broad_LINCS_sig_info_2017-03-06.txt.gz")):
        t, ln, tot = phase(a.meta, label, sig_f, op3_names, op3_blocks, blocks)
        tabs.append(t); lines[label] = ln; totals[label] = tot
    tab = pd.concat(tabs, ignore_index=True).sort_values(["tier", "op3_compounds_present"], ascending=[True, False])
    tab.to_csv(a.out / "cell_line_triage.tsv", sep="\t", index=False)
    pd.concat([v.rename(k) for k, v in lines.items()], axis=1).fillna(0).astype(int).to_csv(
        a.out / "compounds_per_cell_line.tsv", sep="\t")

    gi = genes(a.meta, "GSE92742_Broad_LINCS_gene_info.txt.gz")
    gi.to_csv(a.out / "l1000_genes.tsv", sep="\t", index=False)
    sig = pd.read_csv(ROOT / "analysis/sle_signatures/signatures.tsv", sep="\t")
    cov = []
    for pop, s in sig.groupby("population", sort=False):
        m = gi[gi.gene.isin(s.gene)]
        n_lm = int(m.is_landmark.sum())
        cov.append(dict(population=pop, signature_genes=len(s), in_l1000=len(m), landmark=n_lm,
                        inferred=len(m) - n_lm, landmark_fraction=round(n_lm / len(s), 3)))
    cov = pd.DataFrame(cov)
    cov.to_csv(a.out / "signature_gene_coverage.tsv", sep="\t", index=False)

    act = tab[tab.op3_compounds_present >= MIN_OP3_ACT]
    rep = tab[(tab.op3_compounds_present >= MIN_OP3_REPORT) & (tab.op3_compounds_present < MIN_OP3_ACT)]
    out = [f"LINCS triage. {src}.",
           "Distinct compounds per phase: " + ", ".join(f"{k} {v}" for k, v in totals.items()), "",
           tab.to_string(index=False), "",
           "SLE signature genes in the L1000 12,328-gene space (landmark = measured, inferred = modelled from "
           "the 978 landmarks):", cov.to_string(index=False), "",
           f"Calibratable lines (>= {MIN_OP3_ACT} OP3 compounds): "
           f"{', '.join(act.phase + ':' + act.cell_line) if len(act) else 'none'}",
           f"Report-only lines ({MIN_OP3_REPORT}-{MIN_OP3_ACT - 1} OP3 compounds, interval too wide to act on): "
           f"{', '.join(rep.phase + ':' + rep.cell_line) if len(rep) else 'none'}", ""]
    for _, s in cov.iterrows():
        out.append(f"  [{'ok' if s.landmark_fraction >= 0.5 else 'FAIL'}] {s.population}: {s.landmark} of "
                   f"{s.signature_genes} signature genes are L1000 landmarks, {s.inferred} inferred -> score on "
                   f"landmarks of the full DE, not on this signature")
    for tier in ("immune", "core"):
        t, k = tab[tab.tier == tier], act[act.tier == tier]
        out.append(f"  [{'ok' if len(k) else 'FAIL'}] {tier}: {len(t)} line(s) present, {len(k)} calibratable; "
                   f"largest adds {int(t.compounds.max()) if len(t) else 0} compounds")
    out += ["", "Level 5 matrix (gzipped sizes from the GEO FTP listing, 2026-11-06):",
            "  GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx.gz   20G  (~23 GB unzipped)",
            "  GSE70138_Broad_LINCS_Level5_COMPZ_n118050x12328_2017-03-06.gctx.gz  5.0G  (~5.8 GB unzipped)",
            "A FAIL above is a result to report, not a line to pool anyway."]
    (a.out / "triage.txt").write_text("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
