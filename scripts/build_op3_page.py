#!/usr/bin/env python
"""Build docs/op3-stage-p.html from analysis/op3_stage_p/. Every number is read or recomputed
from those files, never typed here; recomputed values are asserted against the run's own summary.tsv.

Usage: python scripts/build_op3_page.py
"""
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from mol_fingerprints import fingerprints  # noqa: E402

ROOT = Path(__file__).parent.parent
SRC = ROOT / "analysis/op3_stage_p"
OUT = ROOT / "docs/op3-stage-p.html"
LAYERS = {"logFC": "logFC", "clipped_sign_log10_pval": "clipped signed log<sub>10</sub> p"}
BINS = [0, 0.25, 0.4, 1.0]
rng = np.random.default_rng(0)


def f(x, n=3):
    return "&minus;" + f"{-x:.{n}f}" if x < 0 else f"{x:.{n}f}"


def ci(x):
    """Bootstrap mean CI over compounds. Rows of one compound are not independent, so the unit is the compound."""
    x = np.asarray(x, float)
    x = x[~np.isnan(x)]
    m = rng.choice(x, (4000, len(x))).mean(1)
    return x.mean(), *np.percentile(m, [2.5, 97.5])


def interval(est, lo, hi, n=3):
    return f"{f(est, n)} <span class='muted'>({f(lo, n)} to {f(hi, n)})</span>"


def load():
    out = {}
    for layer in LAYERS:
        d = SRC / f"stage_p_{layer}"
        s = pd.read_csv(d / "summary.tsv", sep="\t")
        per = pd.read_csv(d / "per_row.tsv", sep="\t")
        pred = pd.read_csv(d / "prednisolone.tsv", sep="\t")
        verdict = (d / "verdict.txt").read_text()
        w = per[per.evaluation == "whole_compound"]
        byc = w.groupby("sm_name").mean(numeric_only=True)
        # the page must not drift from the run: recompute the headline means and compare
        for metric in ("r.ridge", "r.cell_type_mean", "r.nearest_compound", "r_resid.ridge", "r_resid.nearest_compound"):
            m, c = metric.split(".", 1)
            said = s.query("evaluation=='whole_compound' and metric==@m and comparison==@c").estimate.iat[0]
            assert abs(byc[metric].mean() - said) < 2e-3, f"{layer} {metric}: page {byc[metric].mean():.4f} vs run {said:.4f}"
        out[layer] = dict(summary=s, per=per, byc=byc, pred=pred, verdict=verdict,
                          criteria=dict(re.findall(r"\[(ok|FAIL)\] (\d [^\n]+)", verdict)[::1] and
                                        [(m.group(2), m.group(1) == "ok") for m in re.finditer(r"\[(ok|FAIL)\] (\d [^\n]+)", verdict)]))
    return out


def get(s, evaluation, metric, comparison):
    r = s.query("evaluation==@evaluation and metric==@metric and comparison==@comparison")
    assert len(r) == 1, (evaluation, metric, comparison)
    return r.iloc[0]


def dotplot(bins, colours, values, title_id):
    """Horizontal dot-and-interval: estimate with 95 % CI, one row per bin per series.
    Chosen over bars because these are estimates with uncertainty, not magnitudes from zero."""
    W, left, right, row_h, grp_gap, pad_t, pad_b = 640, 150, 54, 26, 14, 34, 42
    lo = min(v[1] for s in values.values() for v in s)
    hi = max(v[2] for s in values.values() for v in s)
    span = hi - lo or 1
    x0, x1 = lo - span * 0.08, hi + span * 0.08
    plot_w = W - left - right
    sx = lambda v: left + (v - x0) / (x1 - x0) * plot_w  # noqa: E731
    n = len(bins)
    H = pad_t + n * (len(colours) * row_h + grp_gap) + pad_b
    g = [f'<svg viewBox="0 0 {W} {H}" role="img" aria-labelledby="{title_id}-t">',
         f'<title id="{title_id}-t">Compound-specific signal by chemical similarity, with 95 % intervals</title>']
    ticks = np.linspace(0, max(hi, 0), 4) if lo >= -0.02 else np.linspace(lo, hi, 5)
    for t in ticks:
        g.append(f'<line x1="{sx(t):.1f}" y1="{pad_t - 10}" x2="{sx(t):.1f}" y2="{H - pad_b + 4}" '
                 f'stroke="var(--rule)" stroke-width="1"/>')
        g.append(f'<text x="{sx(t):.1f}" y="{H - pad_b + 20}" text-anchor="middle" font-size="11" '
                 f'fill="var(--ink-3)">{t:.2f}</text>')
    g.append(f'<line x1="{sx(0):.1f}" y1="{pad_t - 10}" x2="{sx(0):.1f}" y2="{H - pad_b + 4}" '
             f'stroke="var(--ink-3)" stroke-width="1"/>')
    y = pad_t
    for i, (label, count) in enumerate(bins):
        g.append(f'<text x="0" y="{y + row_h - 8}" font-size="12" fill="var(--ink-2)">{label}</text>')
        g.append(f'<text x="0" y="{y + row_h + 8}" font-size="11" fill="var(--ink-3)">{count} compounds</text>')
        for k, (name, colour) in enumerate(colours.items()):
            est, lo_, hi_ = values[name][i]
            cy = y + row_h / 2 + k * row_h
            g.append(f'<g><title>{name}, {label}: {est:.3f} (95 % CI {lo_:.3f} to {hi_:.3f}), {count} compounds</title>')
            g.append(f'<line x1="{sx(lo_):.1f}" y1="{cy:.1f}" x2="{sx(hi_):.1f}" y2="{cy:.1f}" '
                     f'stroke="{colour}" stroke-width="2" stroke-linecap="round"/>')
            g.append(f'<circle cx="{sx(est):.1f}" cy="{cy:.1f}" r="5" fill="{colour}" '
                     f'stroke="var(--surface)" stroke-width="2"/>')
            g.append(f'<text x="{sx(hi_) + 8:.1f}" y="{cy + 4:.1f}" font-size="11" fill="var(--ink-2)">{est:.3f}</text>')
            g.append("</g>")
        y += len(colours) * row_h + grp_gap
    g.append(f'<text x="{left + plot_w / 2:.0f}" y="{H - 6}" text-anchor="middle" font-size="11" '
             f'fill="var(--ink-3)">correlation after the cell-type mean is removed</text>')
    g.append("</svg>")
    return "\n".join(g)


def main():
    data = load()
    L = data["logFC"]
    run = (SRC / "run_log.txt").read_text()
    head = L["verdict"].splitlines()
    rows, ncomp, ngenes = (int(x) for x in re.search(r"(\d+) rows, (\d+) compounds, (\d+) genes", head[0]).groups())
    commit = re.search(r"commit (\w+)", head[1]).group(1)
    host = re.search(r"host (\S+);", head[1]).group(1)
    seed = re.search(r"seed (\d+)", head[1]).group(1)
    date = re.search(r"date:\s+(.+)", run).group(1).strip()
    folds, perms, boots = (int(x) for x in re.search(r"folds (\d+), permutations (\d+), bootstraps (\d+)", head[2]).groups())
    alphas = re.search(r"alphas chosen: \[(.+)\]", L["verdict"]).group(1)

    comps = pd.read_csv(SRC / "compounds.tsv", sep="\t")
    assert len(comps) == ncomp
    fp, keys, bad = fingerprints(comps.SMILES)
    assert not bad
    inter = fp @ fp.T
    k = fp.sum(1)
    tan = inter / (k[:, None] + k[None, :] - inter)
    np.fill_diagonal(tan, -1)
    nn = tan.max(1)
    seed_tbl = pd.read_csv(ROOT / "shared/compounds/seed_compounds.tsv", sep="\t").dropna(subset=["smiles"])
    _, seed_keys, _ = fingerprints(seed_tbl.smiles)
    seed_blocks = {s[:14] for s in seed_keys if s}
    overlap = sorted(comps.sm_name[[i for i, s in enumerate(keys) if s[:14] in seed_blocks]])

    split = pd.read_csv(SRC / "split_by_cell_type.tsv", sep="\t", index_col=0)
    gc = pd.read_csv(SRC / "glucocorticoid_targets.tsv", sep="\t")

    # --- criteria table -------------------------------------------------------
    crit_names = {"global_mean": "global mean", "cell_type_mean": "cell-type mean",
                  "permuted": "permuted null", "nearest_compound": "nearest training compound"}
    crit_rows = []
    for comp, label in crit_names.items():
        r = get(L["summary"], "whole_compound", "r", f"ridge_minus_{comp}")
        passed = bool(r.ci_lo > 0)
        crit_rows.append(f"<tr><td>Beats the {label}</td><td class='num'>{interval(r.estimate, r.ci_lo, r.ci_hi)}</td>"
                         f"<td class='{'pass' if passed else 'fail'}'>{'pass' if passed else 'FAIL'}</td></tr>")
    pr = L["pred"].query("training == 'mometasone_in_training'").set_index("gene")
    pred_pass = bool(((pr.ridge > 0) & (pr.ridge > pr.cell_type_mean)).all())
    crit_rows.append(
        "<tr><td>Unseen prednisolone in myeloid cells induces TSC22D3 and FKBP5, above the cell-type mean</td>"
        f"<td class='num'>TSC22D3 {f(pr.ridge['TSC22D3'], 2)}, FKBP5 {f(pr.ridge['FKBP5'], 2)}<br>"
        f"<span class='muted'>measured {f(pr.measured['TSC22D3'], 2)} and {f(pr.measured['FKBP5'], 2)}</span></td>"
        f"<td class='{'pass' if pred_pass else 'fail'}'>{'pass' if pred_pass else 'FAIL'}</td></tr>")

    # --- raw r vs compound-specific signal ------------------------------------
    sig_rows = []
    for comp, label in [("ridge", "Ridge on the fingerprint"), ("nearest_compound", "Nearest training compound"),
                        ("cell_type_mean", "Cell-type mean"), ("permuted", "Permuted null")]:
        r = get(L["summary"], "whole_compound", "r", comp)
        rr = L["summary"].query("evaluation=='whole_compound' and metric=='r_resid' and comparison==@comp")
        resid = "not defined" if rr.empty or np.isnan(rr.estimate.iat[0]) else interval(*rr.iloc[0][["estimate", "ci_lo", "ci_hi"]])
        sig_rows.append(f"<tr><td>{label}</td><td class='num'>{interval(r.estimate, r.ci_lo, r.ci_hi)}</td>"
                        f"<td class='num'>{resid}</td></tr>")
    ridge_vs_near = get(L["summary"], "whole_compound", "r_resid", "ridge_minus_nearest_compound")

    # --- similarity bins ------------------------------------------------------
    byc = L["byc"].copy()
    byc["bin"] = pd.cut(byc.tanimoto, BINS)
    bins, series = [], {"Ridge on the fingerprint": "var(--bar)", "Nearest training compound": "var(--bar-2)"}
    series_vals = {k: [] for k in series}
    bin_rows = []
    for b, sub in byc.groupby("bin", observed=True):
        label = f"up to {b.right:g}" if b.left == 0 else (f"above {b.left:g}" if b.right == 1 else f"{b.left:g} to {b.right:g}")
        bins.append((label, len(sub)))
        cells = []
        for name, col in [("Ridge on the fingerprint", "r_resid.ridge"), ("Nearest training compound", "r_resid.nearest_compound")]:
            est, lo, hi = ci(sub[col])
            series_vals[name].append((est, lo, hi))
            cells.append(f"<td class='num'>{interval(est, lo, hi)}</td>")
        bin_rows.append(f"<tr><td>{label}</td><td class='num'>{len(sub)}</td>{''.join(cells)}</tr>")
    chart = dotplot(bins, series, series_vals, "simchart")
    legend = " ".join(f"<span style='color:var(--ink-2)'><span style='display:inline-block;width:10px;height:10px;"
                      f"border-radius:50%;background:{c};vertical-align:-1px;margin-right:5px'></span>{n}</span>"
                      for n, c in series.items())

    # --- prednisolone ---------------------------------------------------------
    pred_rows = []
    for _, r in L["pred"].iterrows():
        label = ("mometasone kept in training" if r.training == "mometasone_in_training" else "no glucocorticoid in training")
        pred_rows.append(f"<tr><td>{label}</td><td>{r.gene}</td><td class='num'>{f(r.measured, 2)}</td>"
                         f"<td class='num'>{f(r.ridge, 2)}</td><td class='num'>{f(r.cell_type_mean, 2)}</td>"
                         f"<td class='num'>{f(r.nearest_compound, 2)}</td>"
                         f"<td>{r.nearest_name} <span class='muted'>({r.tanimoto:.2f})</span></td></tr>")

    # --- official split -------------------------------------------------------
    off_rows = []
    for comp, label in [("ridge", "Ridge on the fingerprint"), ("same_compound", "The same compound, measured in T and NK cells"),
                        ("nearest_compound", "Nearest training compound")]:
        r = get(L["summary"], "official", "r_resid", comp)
        off_rows.append(f"<tr><td>{label}</td><td class='num'>{interval(r.estimate, r.ci_lo, r.ci_hi)}</td></tr>")
    same = get(L["summary"], "official", "r_resid", "same_compound")
    off_ridge = get(L["summary"], "official", "r_resid", "ridge")

    # --- the other layer agrees ----------------------------------------------
    other = data["clipped_sign_log10_pval"]
    other_lines = []
    for comp, label in crit_names.items():
        r = get(other["summary"], "whole_compound", "r", f"ridge_minus_{comp}")
        other_lines.append(f"<tr><td>Beats the {label}</td><td class='num'>{interval(r.estimate, r.ci_lo, r.ci_hi)}</td>"
                           f"<td class='{'pass' if r.ci_lo > 0 else 'fail'}'>{'pass' if r.ci_lo > 0 else 'FAIL'}</td></tr>")
    other_resid = get(other["summary"], "whole_compound", "r_resid", "ridge_minus_nearest_compound")

    split_rows = "".join(f"<tr><td>{i}</td><td class='num'>{int(r.train)}</td><td class='num'>{int(r.test)}</td>"
                         f"<td class='num'>{int(r.train + r.test)}</td></tr>" for i, r in split.iterrows())
    verdicts = {k: ("PASS" if "PASS" in v["verdict"].splitlines()[-1] else "FAIL") for k, v in data.items()}
    assert set(verdicts.values()) == {"FAIL"}, verdicts

    html = f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>OP3 Stage P</title>
<meta name="description" content="Pilot test of whether a compound's chemical structure predicts its transcriptional effect on PBMCs, on the Open Problems OP3 perturbation dataset.">
<style>
{(ROOT / "docs/report.css").read_text().strip()}
</style>
</head>
<body data-palette="#2a78d6,#eb6834">
<main>
<h1>Can a drug's structure predict its effect?</h1>
<p class="sub">Stage P pilot on OP3, {ncomp} compounds in human PBMCs &middot; run {date} on {host} &middot;
commit <code>{commit}</code></p>

<div class="summary">
<p><b>No, not from a fingerprint.</b> A ridge model on 2048-bit Morgan fingerprints does not predict a
held-out compound's effect any better than simply using the average effect across all compounds in that
cell type: a difference of {interval(*get(L['summary'], 'whole_compound', 'r', 'ridge_minus_cell_type_mean')[['estimate', 'ci_lo', 'ci_hi']])}
in correlation. <b>Stage P fails</b>, on both target measures.</p>
<p>The data is not at fault. A compound's measured effect in T and NK cells predicts its own effect in B and
myeloid cells far better than the model does ({f(same.estimate)} against {f(off_ridge.estimate)}), so
response predicts response. What fails is the link from <i>structure</i> to response.</p>
<p>Looking up the most chemically similar compound OP3 already measured beats the model on compound-specific
signal, by {interval(-ridge_vs_near.estimate, -ridge_vs_near.ci_hi, -ridge_vs_near.ci_lo)}. The baseline
built to be beaten is the method that works.</p>
</div>

<div class="tiles">
<div class="tile"><b>{ncomp}</b><span>compounds, 4 cell types, {rows} measured profiles</span></div>
<div class="tile"><b>{np.median(nn):.2f}</b><span>median similarity to the nearest other compound</span></div>
<div class="tile"><b>{(nn >= 0.4).sum()}</b><span>compounds within 0.4 of another <span class="muted">({(nn >= 0.4).mean() * 100:.0f}&thinsp;%)</span></span></div>
<div class="tile"><b>{len(overlap)}</b><span>also in our own compound list</span></div>
</div>

<h2>What was asked, and what came back</h2>
<p>The question was whether structure predicts response well enough to build on. The criteria below were
fixed in writing before the run, so they could fail. <span class="muted">Target: {LAYERS['logFC']}.
Held out whole compounds, {folds}-fold, so a test compound appears in no cell type during training.
&ldquo;Beats&rdquo; means the 95&thinsp;% interval of the paired difference, resampling compounds, excludes zero.</span></p>

<div class="scroll"><table class="plain">
<thead><tr><th>Criterion</th><th class="num">Difference in correlation</th><th>Result</th></tr></thead>
<tbody>{''.join(crit_rows)}</tbody>
</table></div>

<div class="note">
<p><b>The one criterion that passed was badly specified, and I wrote it. <span class="tag k">known</span></b>
Ridge chose heavy regularisation (alpha {alphas}) and so sits almost exactly on the cell-type mean. Raw
correlation rewards that, which means the cell-type mean on its own would also have &ldquo;beaten&rdquo; the
nearest-compound baseline. It is not evidence that the model learned any chemistry.</p>
<p>The criterion is left as written rather than quietly rewritten. The column below that it should have used
is the second one.</p>
</div>

<h2>Raw correlation hides the thing we care about</h2>
<p>Most of any correlation here comes from the cell type, not the compound: all compounds in myeloid cells
look broadly alike. Removing the cell-type mean first leaves only the compound-specific part, which is what a
drug-discovery claim rests on.</p>

<div class="scroll"><table class="plain">
<thead><tr><th>Method</th><th class="num">Correlation</th><th class="num">After removing the cell-type mean</th></tr></thead>
<tbody>{''.join(sig_rows)}</tbody>
</table></div>
<p class="muted">The cell-type mean has no compound-specific part by construction, and the permuted null keeps
every structure but shuffles which compound it belongs to.</p>

<h2>Chemical similarity is what carries the signal <span class="tag g">exploratory</span></h2>
<p>This split was chosen after seeing the results, so treat it as a lead rather than a test. Binned by how
close each compound is to its nearest training neighbour, compound-specific signal appears only when
something chemically similar was already measured &mdash; for both methods.</p>

<figure>
<p class="muted" style="margin:0 0 4px">{legend}</p>
{chart}
<figcaption>Correlation after the cell-type mean is removed, with 95&thinsp;% intervals over compounds.
Chemical similarity is Tanimoto on the same fingerprints the model used. Intervals, not point estimates,
because {len(byc)} compounds is not many; the top bin holds only {bins[-1][1]}.</figcaption>
</figure>

<div class="scroll"><table class="plain">
<thead><tr><th>Similarity to nearest training compound</th><th class="num">Compounds</th>
<th class="num">Ridge</th><th class="num">Nearest compound</th></tr></thead>
<tbody>{''.join(bin_rows)}</tbody>
</table></div>

<h2>The steroid test, in detail</h2>
<p>Prednisolone induces TSC22D3 and FKBP5 through the glucocorticoid receptor &mdash; the most reliable
pharmacology in this dataset, and the easiest thing a working model should recover. With every prednisolone
row removed from training, mometasone furoate remains, a steroid {L['pred'].tanimoto.max():.2f} similar.</p>

<div class="scroll"><table class="plain">
<thead><tr><th>Training set</th><th>Gene</th><th class="num">Measured</th><th class="num">Ridge</th>
<th class="num">Cell-type mean</th><th class="num">Nearest compound</th><th>Which compound</th></tr></thead>
<tbody>{''.join(pred_rows)}</tbody>
</table></div>
<p>The model predicts almost nothing. Looking up mometasone lands close to the measured value. With no
steroid at all in training the lookup fails too, which is the honest limit of that method: it can only
transfer from chemistry it has seen.</p>

<h2>The dataset, and why the published split is not a test of new chemistry</h2>
<p>OP3's own split holds out B and myeloid cells, not compounds: all {int(split.test.max())} test compounds
also appear in training in T and NK cells. It measures transfer between cell types. That is why the verdict
above uses a whole-compound holdout instead, and why the published split can be used to show that compound
signal is present.</p>

<div class="scroll"><table class="plain">
<thead><tr><th>Cell type</th><th class="num">Compounds in train</th><th class="num">In test</th><th class="num">Total</th></tr></thead>
<tbody>{split_rows}</tbody>
</table></div>

<h3>On the published split, compound-specific signal</h3>
<div class="scroll"><table class="plain">
<thead><tr><th>Method</th><th class="num">After removing the cell-type mean</th></tr></thead>
<tbody>{''.join(off_rows)}</tbody>
</table></div>
<p>Carrying a compound's own measured profile across cell types beats the structure model. The signal is
there to be predicted; fingerprints are not the way to reach it.</p>

<h2>Both target measures agree</h2>
<p>The run was repeated against {LAYERS['clipped_sign_log10_pval']}, the measure the published benchmark
scores on. {LAYERS['logFC']} decides the verdict here because p-value-based measures scale with how many
cells a cell type contributed, while an effect size does not.</p>
<div class="scroll"><table class="plain">
<thead><tr><th>Criterion, on {LAYERS['clipped_sign_log10_pval']}</th><th class="num">Difference in correlation</th><th>Result</th></tr></thead>
<tbody>{''.join(other_lines)}</tbody>
</table></div>
<p>Nearest compound again wins on compound-specific signal, by
{interval(-other_resid.estimate, -other_resid.ci_hi, -other_resid.ci_lo)}.</p>

<h2>What this changes</h2>
<ul>
<li><b>Do not scale the same model.</b> The rule set before the run was that a bigger model of the same kind
does not follow a failure like this. <span class="tag k">known</span></li>
<li><b>Test the molecule representation next, holding everything else fixed.</b> Same folds, same seed, same
baselines, with a pretrained chemical language model in place of the fingerprint. If that also cannot beat
nearest compound, the limit is {ncomp} compounds of sparse chemistry rather than the encoder.
<span class="tag i">inferred</span></li>
<li><b>Denser chemistry matters more than a bigger model.</b> Pretraining on LINCS L1000, roughly 30&thinsp;000
compounds, is the only route to enough chemical coverage for a structure-to-response map. This result is
evidence for that, not against it. <span class="tag i">inferred</span></li>
<li><b>For the lupus application today, use the lookup.</b> For any candidate within about 0.4 of something
OP3 measured, the defensible answer is that neighbour's measured profile. {len(overlap)} of our own
candidate compounds were measured in OP3 directly, so for those there is no prediction to make.
<span class="tag k">known</span></li>
</ul>

<h2>Honest limits</h2>
<ul>
<li><b>One seed, one set of folds.</b> Seed {seed}; stability across seeds was not tested.</li>
<li><b>No noise ceiling.</b> These summarised tables carry no replicates, so how much of the unexplained
variation is simply measurement noise cannot be estimated from them.</li>
<li><b>The similarity bins were chosen after seeing results</b> and the intervals are not adjusted for having
looked. Treat that section as a lead.</li>
<li><b>Healthy cells, in culture.</b> OP3 is PBMCs from 3 healthy donors at a single dose for 24 hours.
Nothing here tests prediction in diseased cells, which is where this pipeline is aimed.</li>
<li><b>A linear model on a binary fingerprint is a floor, not a fair test of chemistry.</b> It cannot
represent a compound acting differently in different cell types. That is deliberate for a baseline, and it
is why the next step is the representation.</li>
</ul>

<h2>Reproducing this</h2>
<p>On the compute node, one command runs the content checks, then the pilot, and packs the outputs:</p>
<p><code>scripts/op3_report.sh</code></p>
<p>Inputs are the Open Problems build of OP3 (<code>de_train.h5ad</code>, <code>de_test.h5ad</code>), fetched by
<code>scripts/fetch_op3.sh</code>. {folds} folds, {perms} permutations for the null, {boots} bootstrap
resamples, seed {seed}. The model, baselines and criteria are in <code>scripts/op3_ridge.py</code>, with a
self-test in <code>scripts/test_op3_ridge.py</code> that requires a planted signal to be recovered and pure
noise not to be. The outputs this page is built from are in <code>analysis/op3_stage_p/</code>, and every
number here is read or recomputed from them by <code>scripts/build_op3_page.py</code>, which fails if a
recomputed value disagrees with the run's own summary.</p>

<footer>
<p>Diagnostic report, not a result for publication. OP3 is Szalata et&nbsp;al., NeurIPS 2024 Datasets and
Benchmarks (GEO GSE279945), CC&nbsp;BY&nbsp;4.0. {rows} profiles, {ngenes} genes.
Glucocorticoid direction confirmed in all {len(gc)} steroid measurements before any modelling.
&middot; <a href="index.html">All reports</a></p>
</footer>
</main>
</body>
</html>
"""
    OUT.write_text(html)
    print(f"wrote {OUT} ({len(html)} bytes); verdicts {verdicts}; overlap {len(overlap)}: {', '.join(overlap)}")
    print(f"median nn tanimoto {np.median(nn):.3f}; >=0.4 {(nn >= 0.4).sum()}")


if __name__ == "__main__":
    main()
