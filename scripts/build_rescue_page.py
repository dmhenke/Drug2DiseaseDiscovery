#!/usr/bin/env python
"""Build the chemical rescue lookup page from the output directory of scripts/op3_rescue.py.

Usage: python scripts/build_rescue_page.py analysis/op3_rescue [--out docs/chemical-rescue.html] [--synthetic]
--synthetic stamps a banner on the page; it is for previewing the interface, never for publishing.
"""
import argparse
import json
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from mol_fingerprints import fingerprints  # noqa: E402

ROOT = Path(__file__).parent.parent


def seed_flags(smiles):
    """Which compounds are also in our seed table, matched on InChIKey connectivity."""
    seed = pd.read_csv(ROOT / "shared/compounds/seed_compounds.tsv", sep="\t").dropna(subset=["smiles"])
    _, seed_keys, _ = fingerprints(seed.smiles)
    labels = {"positive_control": "lupus positive control", "l1000_hit": "L1000 reversal hit"}
    block = {}
    for k, src in zip(seed_keys, seed.source):
        if k:
            block.setdefault(k[:14], set()).add(labels.get(src, src))
    _, keys, _ = fingerprints(smiles)
    return [", ".join(sorted(block.get(k[:14], []))) if k else "" for k in keys]


def main():
    p = argparse.ArgumentParser()
    p.add_argument("src", type=Path)
    p.add_argument("--out", type=Path, default=ROOT / "docs/chemical-rescue.html")
    p.add_argument("--synthetic", action="store_true")
    a = p.parse_args()
    res = pd.read_csv(a.src / "rescue.tsv", sep="\t")
    genes = pd.read_csv(a.src / "rescue_genes.tsv", sep="\t")
    cov = pd.read_csv(a.src / "coverage.tsv", sep="\t")
    summary = (a.src / "run_summary.txt").read_text().splitlines()[0]
    commit = re.search(r"commit (\w+)", summary).group(1)
    assert len(res) == genes.groupby(["population", "sm_name"]).ngroups, "score and gene tables disagree"

    comp = res.drop_duplicates("sm_name")[["sm_name", "SMILES"]]
    flags = dict(zip(comp.sm_name, seed_flags(comp.SMILES)))
    pops = {}
    for pop, g in genes.groupby("population", sort=False):
        first = g[g.sm_name == g.sm_name.iat[0]]
        c = cov[cov.population == pop].iloc[0]
        pops[pop] = dict(genes=first.gene.tolist(), d=first.disease.round(2).tolist(), cell=c.op3_cell_type,
                         total=int(c.signature_genes), compounds=int(c.compounds))
    drug = {k: v.round(2).tolist() for k, v in genes.groupby(["population", "sm_name"], sort=False).drug}
    rows = [dict(p=r.population, c=r.sm_name, s=r.SMILES, e=r.evidence, g=round(r.match_gain, 3), lo=round(r.ci_lo, 3),
                 hi=round(r.ci_hi, 3), q=round(r.q_bh, 4), n=int(r.genes_toward_healthy), m=int(r.genes_scored),
                 cos=round(r.reversal_cosine, 3), off=round(r.off_target_rms, 3), f=flags[r.sm_name],
                 t=drug[(r.population, r.sm_name)]) for r in res.itertuples()]
    payload = json.dumps(dict(pops=pops, rows=rows), separators=(",", ":")).replace("</", "<\\/")
    css = (ROOT / "docs/report.css").read_text().strip()
    banner = ("<div class='note' style='border-left-color:var(--fail)'><p><b>SYNTHETIC DATA.</b> Interface preview only. "
              "These responses are invented.</p></div>") if a.synthetic else ""
    html = TEMPLATE.replace("/*CSS*/", css).replace("/*DATA*/", payload).replace("<!--BANNER-->", banner) \
                   .replace("<!--COMMIT-->", commit).replace("<!--NCOMP-->", str(res.sm_name.nunique()))
    a.out.write_text(html)
    print(f"wrote {a.out}: {len(rows)} rows, {res.sm_name.nunique()} compounds, {len(html) // 1024} KB"
          + (" (SYNTHETIC)" if a.synthetic else ""))


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Chemical Rescue Lookup</title>
<meta name="description" content="Filter measured compounds by how far each one would move an SLE cell population toward its healthy state.">
<style>
/*CSS*/
.controls { display: grid; grid-template-columns: repeat(auto-fit, minmax(170px, 1fr)); gap: 12px; margin: 20px 0 10px; }
.controls label { display: block; font-size: 0.8rem; font-weight: 600; color: var(--ink-2); }
.controls input:not([type=checkbox]), .controls select { display: block; width: 100%; margin-top: 4px; padding: 8px; font: inherit; color: var(--ink);
  background: var(--surface); border: 1px solid var(--ink-3); border-radius: 6px; }
.check { display: flex !important; gap: 8px; align-items: center; align-self: end; padding-bottom: 8px; }
.status { display: flex; flex-wrap: wrap; justify-content: space-between; gap: 8px; align-items: center; color: var(--ink-2); font-size: 0.88rem; margin: 8px 0; }
button { font: inherit; color: var(--ink); }
.dl { padding: 6px 12px; background: var(--surface); border: 1px solid var(--ink-3); border-radius: 6px; cursor: pointer; }
th button { all: unset; cursor: pointer; font-weight: 600; font-size: 0.82rem; color: var(--ink-2); }
th button:focus-visible, .open:focus-visible, .dl:focus-visible { outline: 2px solid var(--known); outline-offset: 2px; }
.open { all: unset; cursor: pointer; font-weight: 600; color: var(--known); }
#tbl td:first-child { min-width: 150px; }
.bar { position: relative; width: 120px; height: 8px; margin-top: 5px; background: var(--chip); border-radius: 4px; }
.bar i { position: absolute; top: 0; height: 8px; }
.bar::after { content: ""; position: absolute; left: 50%; top: -2px; height: 12px; width: 1px; background: var(--ink-3); }
.detail td { background: var(--surface); }
.genes { display: grid; grid-template-columns: repeat(auto-fill, minmax(200px, 1fr)); gap: 2px 14px; font-size: 0.82rem; margin: 8px 0; }
.smiles { font: 0.78rem ui-monospace, Menlo, monospace; overflow-wrap: anywhere; color: var(--ink-2); }
.flag { color: var(--ink-3); font-size: 0.8rem; display: block; font-weight: 400; }
.warn { color: var(--guess); }
</style>
</head>
<body>
<main>
<!--BANNER-->
<h1>Chemical rescue lookup</h1>
<p class="sub">How far would each compound move an SLE cell population toward healthy? <!--NCOMP--> compounds measured in OP3 &middot; commit <code><!--COMMIT--></code> &middot; <a href="index.html">All reports</a></p>
<div class="summary">
<p style="margin-top:0"><b>Match to healthy</b> is 1 &minus; &Vert;disease + drug&Vert; / &Vert;disease&Vert; across the population's SLE signature genes. <b>1</b> means the signature is fully closed, <b>0</b> means no change, and <b>below 0</b> means further from healthy.</p>
<p style="margin-bottom:0"><b>This is a calculation, not an observation.</b> No lupus cell here was treated. Drug effects were measured in cultured healthy PBMCs and added to the disease signature, which assumes the two effects add. OP3's coarse cell types stand in for finer populations: both monocyte groups use Myeloid. Use it to choose compounds to test, not as evidence that a drug works in lupus.</p>
</div>
<div class="controls">
<div><label for="pop">Cell population</label><select id="pop"></select></div>
<div><label for="name">Compound</label><input id="name" type="search" placeholder="e.g. ruxolitinib"></div>
<div><label for="gene">Moves this gene toward healthy</label><input id="gene" type="search" list="genelist" placeholder="e.g. IFI44L"><datalist id="genelist"></datalist></div>
<div><label for="min">Minimum match to healthy</label><input id="min" type="number" step="0.05" min="-5" max="1" placeholder="any"></div>
<div><label for="off">Maximum off-target effect</label><input id="off" type="number" step="0.05" min="0" placeholder="any"></div>
<label class="check"><input id="sig" type="checkbox"> Better than chance only (q &lt; 0.05, above 0)</label>
</div>
<p id="thin" class="warn" hidden></p>
<div class="status"><span id="count" role="status"></span><button class="dl" id="dl">Download filtered CSV</button></div>
<div class="scroll"><table class="plain" id="tbl"><thead><tr id="head"></tr></thead><tbody id="body"></tbody></table></div>
<h2>Reading the columns</h2>
<ul>
<li><b>Match to healthy</b> has a 95&thinsp;% interval from resampling signature genes. It reflects which genes are in the signature, not donor or experimental variation.</li>
<li><b>q</b> compares the compound with the same disease values placed on random non-signature genes, using Benjamini-Hochberg across compounds within the population.</li>
<li><b>Off-target</b> is the root-mean-square log<sub>2</sub> change on genes outside the signature. A high value next to a good match can mean a broadly disruptive or toxic compound rather than a specific one.</li>
<li><b>Evidence</b> is &ldquo;measured&rdquo; for every row: OP3 profiled the compound. Structure-based predictions are not shown, because the <a href="op3-stage-p.html">Stage P model</a> failed its test.</li>
</ul>
<footer>Diagnostic tool, not a result for publication and not treatment guidance. Drug data: OP3, Szalata et al., NeurIPS 2024 Datasets and Benchmarks (GSE279945), CC BY 4.0. Disease data: Perez 2022 cohort 4, 20 vs 20 donors.</footer>
</main>
<script>
const D = /*DATA*/;
const $ = id => document.getElementById(id);
const esc = s => String(s).replace(/[&<>"']/g, c => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[c]));
const num = (x, n = 3) => (x < 0 ? "−" : "") + Math.abs(x).toFixed(n);
const COLS = [["c","Compound"],["g","Match to healthy"],["n","Genes toward healthy"],["cos","Reversal cosine"],["q","q"],["off","Off-target"],["e","Evidence"]];
let sort = {k: "g", dir: -1}, open = null, shown = [];
for (const p of Object.keys(D.pops)) $("pop").add(new Option(`${p} (${D.pops[p].genes.length} genes)`, p));
function geneIndex() { return D.pops[$("pop").value].genes.findIndex(g => g.toUpperCase() === $("gene").value.trim().toUpperCase()); }
function render() {
  const pop = $("pop").value, P = D.pops[pop], name = $("name").value.trim().toLowerCase(), gi = geneIndex();
  const typed = $("gene").value.trim(), min = parseFloat($("min").value), off = parseFloat($("off").value);
  $("genelist").innerHTML = P.genes.map(g => `<option value="${esc(g)}">`).join("");
  const thin = P.genes.length < 15 ? `${pop} has only ${P.genes.length} signature genes, so its scores are unstable. ` : "";
  const missing = typed && gi < 0 ? `${typed} is not in the ${pop} signature, so the gene filter is not applied. ` : "";
  $("thin").textContent = thin + missing; $("thin").hidden = !(thin + missing);
  shown = D.rows.filter(r => r.p === pop && (!name || r.c.toLowerCase().includes(name))
    && (isNaN(min) || r.g >= min) && (isNaN(off) || r.off <= off) && (!$("sig").checked || (r.q < 0.05 && r.g > 0))
    && (gi < 0 || Math.abs(P.d[gi] + r.t[gi]) < Math.abs(P.d[gi])));
  shown.sort((a, b) => (typeof a[sort.k] === "string" ? a[sort.k].localeCompare(b[sort.k]) : a[sort.k] - b[sort.k]) * sort.dir);
  $("head").innerHTML = COLS.map(([k, l]) => `<th aria-sort="${sort.k === k ? (sort.dir > 0 ? "ascending" : "descending") : "none"}" class="${k === "c" || k === "e" ? "" : "num"}"><button data-k="${k}">${l}${sort.k === k ? (sort.dir > 0 ? " ▲" : " ▼") : ""}</button></th>`).join("");
  $("count").textContent = `${shown.length} of ${P.compounds} compounds. ${pop} uses OP3 ${P.cell}; ${P.genes.length} of ${P.total} signature genes measured.`;
  $("body").innerHTML = shown.length ? shown.map(row).join("") : `<tr><td colspan="7">No compound meets these filters. Loosen one to see more.</td></tr>`;
}
function row(r) {
  const c = Math.max(-1, Math.min(1, r.g)), w = Math.abs(c) * 50;
  const bar = `<div class="bar" aria-hidden="true"><i style="left:${c < 0 ? 50 - w : 50}%;width:${w}%;background:var(${c < 0 ? "--bar-2" : "--bar"})"></i></div>`;
  const key = r.p + "|" + r.c;
  let h = `<tr><td><button class="open" data-key="${esc(key)}" aria-expanded="${open === key}">${esc(r.c)}</button>${r.f ? `<span class="flag">${esc(r.f)}</span>` : ""}</td>
  <td class="num">${num(r.g)} <span class="muted">(${num(r.lo, 2)} to ${num(r.hi, 2)})</span>${bar}</td><td class="num">${r.n} / ${r.m}</td><td class="num">${num(r.cos, 2)}</td><td class="num">${r.q < 0.001 ? "&lt;0.001" : r.q.toFixed(3)}</td><td class="num">${r.off.toFixed(2)}</td><td>${esc(r.e)}</td></tr>`;
  if (open === key) {
    const P = D.pops[r.p], items = P.genes.map((g, i) => ({g, d: P.d[i], a: P.d[i] + r.t[i]})).sort((x, y) => (Math.abs(x.a) - Math.abs(x.d)) - (Math.abs(y.a) - Math.abs(y.d)));
    h += `<tr class="detail"><td colspan="7"><b>SLE vs healthy log<sub>2</sub> change, before &rarr; after adding the drug.</b> Closer to 0 is closer to healthy. Most improved first.
    <div class="genes">${items.map(x => `<span>${esc(x.g)}: ${num(x.d, 2)} &rarr; ${num(x.a, 2)}${Math.abs(x.a) < Math.abs(x.d) ? "" : " <span class='warn'>(worse)</span>"}</span>`).join("")}</div>
    <div class="smiles">SMILES: ${esc(r.s)}</div></td></tr>`;
  }
  return h;
}
$("head").addEventListener("click", e => { const k = e.target.closest("button")?.dataset.k; if (!k) return;
  sort = sort.k === k ? {k, dir: -sort.dir} : {k, dir: k === "c" || k === "e" || k === "q" || k === "off" ? 1 : -1}; render(); });
$("body").addEventListener("click", e => { const b = e.target.closest(".open"); if (!b) return; open = open === b.dataset.key ? null : b.dataset.key; render(); });
for (const id of ["pop","name","gene","min","off","sig"]) { $(id).addEventListener("input", () => { if (id === "pop") { $("gene").value = ""; open = null; } render(); }); }
$("dl").addEventListener("click", () => {
  const head = "population,compound,match_to_healthy,ci_lo,ci_hi,genes_toward_healthy,genes_scored,reversal_cosine,q,off_target,evidence,smiles";
  const q = v => { const s = String(v); return /^[=+\-@\t\r]/.test(s) && isNaN(Number(s)) ? `"'${s.replace(/"/g, '""')}"` : `"${s.replace(/"/g, '""')}"`; };
  const body = shown.map(r => [r.p, r.c, r.g, r.lo, r.hi, r.n, r.m, r.cos, r.q, r.off, r.e, r.s].map(q).join(","));
  const a = document.createElement("a"); a.href = URL.createObjectURL(new Blob([[head, ...body].join("\n")], {type: "text/csv"}));
  a.download = "chemical_rescue_filtered.csv"; a.click(); URL.revokeObjectURL(a.href);
});
render();
</script>
</body>
</html>
"""

if __name__ == "__main__":
    main()
