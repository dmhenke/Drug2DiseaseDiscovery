# OP3 pilot and full run

**Decided 2026-11-02:** train on OP3; start with Morgan fingerprints.
Design context: [`chemical_perturbation_model.md`](chemical_perturbation_model.md).

Two stages. **Stage P** is a pilot that can be run in an afternoon on a laptop and
can fail cheaply. **Stage F** is the full run: pretrained molecule encoders, L1000
pretraining, single-cell resolution, on the node's GPU.

**Stage P failed on 2026-11-04** (result below). Under the rule set before running, Stage F does not
start on a bigger model of the same kind. What the failure changes is recorded under Stage F.

## What OP3 actually contains (checked 2026-11-02)

From Szalata et al., NeurIPS 2024 D&B, and the OP3 data files (**Known** unless marked):

- **146 compounds** in PBMCs from **3 healthy donors** (2 male, 1 female), 24 h, **1 µM**.
- **4 cell types** after their QC: T, B, NK, Myeloid. 301,785 cells, 21,265 genes; the benchmark build filters to **5,317 genes** and **140 compounds** (3 dropped for toxicity).
- **Controls:** DMSO vehicle; belinostat (HDAC, 0.1 µM) and dabrafenib (BRAF) as positive controls, in every multiplexed library.
- **SMILES ship with the data** (`obs["SMILES"]`, `sm_name`, `sm_lincs_id`), so compounds map straight to LINCS. No structure lookup needed for training.
- **Licence CC BY 4.0.** Compounds were chosen for being in LINCS CMap and for giving a robust response in CD34+ cells.

**Two things this settles from the earlier design doc:**

1. **The interferon/steroid axis is present after all.** The 76-compound test split alone contains **Ruxolitinib** (JAK1/2) and **Prednisolone**, plus Mometasone Furoate. Both ruxolitinib and prednisolone are in our positive-control list. The earlier worry that OP3 could not teach "JAK inhibitor → interferon genes down" was too pessimistic, though cultured healthy PBMCs still have low baseline interferon, so the effect may be small (**Inferred**).
2. **The 14.2 GB GEO file is not the right download.** The Open Problems benchmark build is the same experiment, QC'd and pseudobulked, at **1.1 GB in single-cell form** and **292 MB** for the DE tables alone (**Known**, from file sizes).

**Confirmed overlap with our own work (Known, matched on InChIKey connectivity, 2026-11-04):** 15 of the 140 OP3 compounds are in our seed table: BMS-387032, CHIR-99021, crizotinib, dactolisib, IMD-0354, LDN 193189, methotrexate, penfluridol, prednisolone, ruxolitinib, C646 (STK219801), TL_HRAS26, TPCA-1, tivantinib, vorinostat.

**Content check on the node, 2026-11-04 (commit 8790252): PASS.** 553 rows: T 140, NK 139, B 139, Myeloid 135 compounds. In all 16 steroid rows (prednisolone and mometasone, 4 cell types, 2 genes) TSC22D3 and FKBP5 are up, logFC 2.3 to 5.0, so the DE sign is as expected.

## Stage P — pilot: fingerprints → pseudobulk response

**Question:** given a compound's structure and a cell type, can we predict its
transcriptional effect well enough to be worth building on?

**Data:** `de_train.h5ad` / `de_test.h5ad` (292 MB). These are limma DE statistics per
(compound, cell type) against DMSO, 5,317 genes. No single cells, no GPU.

**Model:** cell-type mean + ridge on the fingerprint (2048) → 5,317 DE values. The cell-type term is
unpenalised, so strong regularisation falls back to the cell-type mean. Alpha is chosen by inner
compound-grouped CV. This is deliberately the dumbest thing that uses structure.

**Splits, corrected 2026-11-04 from the real data:** OP3's official split does **not** hold out
chemistry. All 76 test compounds have T and NK rows in train, so it measures transfer between cell
types. The verdict therefore uses our own **5-fold whole-compound holdout** over all 553 rows (**primary**);
the official split is reported alongside with an added `same_compound` baseline, which should be hard to beat.

**Target (fixed before running):** `logFC`, because p-based values scale with the cells per cell type
(prednisolone's TSC22D3 logFC is 2.35 in myeloid and 2.97 in T, but its signed log10 p is 39 and 83).
`clipped_sign_log10_pval`, the published benchmark's target, is run as a secondary.

**Baselines it must beat.** Without these the result is meaningless:

| Baseline | What it tests |
|---|---|
| Global mean DE across training compounds | Is there any compound-specific signal at all? |
| Cell-type mean DE | Is the structure doing anything the cell type doesn't? |
| **Nearest training compound** by Tanimoto | Is the model learning chemistry, or just copying the closest neighbour? |
| Random fingerprint permutation | The null. |

**Pass criteria, fixed before running:**

"Beats" means the 95 % bootstrap CI of the paired difference, resampling compounds, is above zero.

1. The ridge model beats the global mean, the cell-type mean and the permuted null on held-out compounds, by mean correlation across compounds.
2. It beats the nearest-training-compound baseline. *(If it does not, say so and use the nearest-neighbour method; it is simpler and would be the honest answer.)*
3. **Biology check:** with all prednisolone rows removed from training, predicted TSC22D3 and FKBP5 in myeloid cells are up **and above the cell-type mean**. Mometasone stays in training (Tanimoto 0.29, prednisolone's nearest neighbour), so this tests transfer from one steroid to another. A run with no glucocorticoid in training is also reported, for information.
4. **Bootstrap CIs** over held-out compounds, not point estimates.

**Run:** `scripts/op3_report.sh` on the node runs the checks, then this pilot, and packs the results.
Tested on synthetic values only: planted signal is recovered, pure noise falls back to the cell-type mean.

### Stage P result: FAIL (node, 2026-11-04, commit 20a5847, seed 0)

Published as [`docs/op3-stage-p.html`](../docs/op3-stage-p.html); outputs in `analysis/op3_stage_p/`.

Whole-compound holdout, 140 compounds, CIs over compounds. Both target layers give the same verdict; logFC shown.

| Criterion | Result | Verdict |
|---|---|---|
| 1. Beats cell-type mean, r | difference 0.001 (-0.002 to 0.005) | **FAIL** |
| 1. Beats global mean, r | difference 0.071 (0.064 to 0.077) | ok |
| 1. Beats permuted null, r | difference 0.005 (0.001 to 0.008) | ok |
| 2. Beats nearest compound, r | difference 0.144 (0.116 to 0.172) | ok, but see below |
| 3. Unseen prednisolone, myeloid | TSC22D3 0.21, FKBP5 -0.04; measured 2.35 and 2.38 | **FAIL** |

**Criterion 2 was badly specified, and its "ok" should not be read as a win.** Ridge chose heavy regularisation
(alpha 1,000 to 10,000) and sits almost on the cell-type mean (r 0.336 vs 0.335), which raw r rewards; the
cell-type mean alone would pass the same test. On compound-specific signal (`r_resid`, correlation after the
cell-type mean is removed) nearest compound wins: 0.095 against ridge 0.048, difference -0.047 (-0.077 to -0.018).
The criterion is left as written rather than rewritten after the fact.

**Exploratory, decided after seeing the results.** Compound-specific signal tracks chemical proximity to a
training compound. `r_resid` binned by nearest-neighbour Tanimoto:

| Tanimoto | Compounds | Ridge | Nearest compound |
|---|---|---|---|
| up to 0.25 | 67 | 0.017 (-0.018 to 0.051) | 0.037 (0.003 to 0.072) |
| 0.25 to 0.4 | 60 | 0.053 (0.006 to 0.100) | 0.138 (0.083 to 0.195) |
| above 0.4 | 13 | 0.183 (0.107 to 0.252) | 0.198 (0.131 to 0.262) |

Nearest compound recovers known mechanism pairs: belinostat and resminostat (HDAC, 0.77 and 0.76), AT13387 and
ganetespib (HSP90, 0.70). For unseen prednisolone it picks mometasone and predicts TSC22D3 3.15, FKBP5 2.47,
against measured 2.35 and 2.38.

**Not checked:** other seeds or fold assignments, and the noise ceiling, which these DE tables cannot give
because they carry no replicates.

**What this does and does not say.** It says a 2048-bit fingerprint with a linear model cannot predict a held-out
compound's effect on PBMCs beyond the cell-type average. It does not say the data lack compound signal: the
official split's `same_compound` baseline reaches `r_resid` 0.269, well above ridge's 0.189, so a compound's
effect in T and NK cells does predict its effect in B and myeloid cells. The failure is structure to response,
not response to response.

## Stage F — full run

Three changes, each testable on its own against Stage P:

| Change | What it buys | Cost |
|---|---|---|
| **F1. Pretrained molecule encoders** instead of fingerprints: ChemBERTa (the NeurIPS competition winner used it), plus RDKit descriptors and a GNN | Chemistry the fingerprint cannot express. Compare on the same splits | Encoder inference only; minutes on GPU |
| **F2. L1000 pretraining, then fine-tune on OP3** — chemCPA's central claim | ~30k compounds of bulk signal to generalize from. Needs the L1000 level-5 matrix (several GB) and gene-space alignment with OP3's 5,317 | Hours on the P4000 (**Guess**) |
| **F3. Single-cell chemCPA** on `sc_train.h5ad` (893 MB, ~300k cells): latent `z_treated = z_basal + g(m, dose)`, negative-binomial decoder | Per-cell predictions, which is what the SLE application needs | MLP autoencoder; should fit in 8 GB (**Guess**, untested) |

**Each must beat Stage P on the same held-out compounds**, on `r_resid` rather than raw r, since Stage P
showed raw r is passed by the cell-type mean alone. Stage P's nearest-compound baseline, not its ridge
model, is the bar: it is the stronger of the two on compound-specific signal.

**Revised after the Stage P failure (2026-11-04).** The order changes, because the failure points at the
molecule representation, not at model size:

- **F1 is now the next step, and is a test of representation, not of capacity.** Swap the fingerprint for
  ChemBERTa, keeping everything else identical: same folds, same seed, same baselines, same ridge. If a
  pretrained encoder also cannot beat nearest compound on `r_resid`, the limit is 140 compounds of sparse
  chemistry, and F2 and F3 should not be attempted on OP3 alone.
- **F2 gains in importance.** L1000's roughly 30k compounds are the only route to dense enough chemistry for
  a structure-to-response map, and Stage P is evidence for that claim rather than against it.
- **F3 is unaffected in principle** but has nothing to condition on until F1 or F2 produces a usable molecule
  term. A single-cell model with no working chemistry term would reduce to a cell-type average per compound.

**Dose:** OP3 is single-dose (1 µM), so the dose term `g(m, d)` cannot be fitted from
OP3 alone. Either fix dose, or fit it from L1000/sciPlex, which are multi-dose
(**Inferred**).

## The SLE application (after F3)

1. Embed SLE and healthy classical monocytes from `studies/lupus_test` with the trained encoder.
2. For each compound in `shared/compounds/`, generate treated SLE cells.
3. Score: fraction of the 68-gene SLE monocyte signature moved toward healthy, and latent distance to the healthy centroid.
4. **Null:** 100 random library compounds. **Positive controls:** ruxolitinib and prednisolone, which OP3 measured directly, so their true effect in myeloid cells is known.
5. **The honest caveat throughout:** OP3 cells are healthy PBMCs in culture, not SLE patient cells. Predicting how a drug acts on an interferon-high SLE monocyte is extrapolation beyond the training distribution.

## Risks, in the order they are likely to bite

1. **Compound-specific signal may be weak.** OP3 itself reports that many compounds give subtle effects and that DE is sensitive to low-count genes. Stage P criterion 1 is exactly this test, which is why it comes first.
2. **Nearest-neighbour wins (confirmed 2026-11-04, not a prediction any more).** On compound-specific signal
   it beats ridge by 0.047 (0.018 to 0.077). For any compound within Tanimoto 0.4 of something OP3 measured,
   the honest method today is to look up that neighbour's measured profile.
3. **Chemistry is sparse (Known, from `compounds.tsv`).** Median nearest-neighbour Tanimoto among the 140 compounds is 0.26; only 11 % reach 0.4. Ruxolitinib's nearest neighbour is crizotinib at 0.18, and filgotinib (GLPG0634) is not close to it despite sharing a target. Fingerprints may not connect compounds that act alike. Expect higher numbers on OP3's own split than on the whole-compound holdout, because the compound is already seen in other cell types.
4. **Healthy PBMC → SLE monocyte transfer is unvalidated** and probably cannot be validated without treated SLE cells.
5. **Three donors** means donor effects cannot be separated from compound effects well. Keep donor as a covariate and never as a split.

## Commands

```bash
scripts/fetch_op3.sh                         # ~1.1 GB, resumable, size-checked
scripts/op3_report.sh                        # checks, Stage P, one tarball in ~/
```

`scripts/mol_fingerprints.py` already passes its analogue check (prednisolone–methyl­prednisolone
0.70 and baricitinib–ruxolitinib 0.48, against 0.09 and 0.12 for unrelated pairs) and has
been run on the 129-compound seed table: 2048-bit vectors, InChIKeys round-tripping 100 %
against PubChem's (**Known**).

The Stage P script is `scripts/op3_ridge.py`, with a self-test in `scripts/test_op3_ridge.py`.

## Sources

- OP3: [NeurIPS 2024 D&B paper](https://papers.neurips.cc/paper_files/paper/2024/file/24c4d51f3ef48dd2dbab78243ecb26a1-Paper-Datasets_and_Benchmarks_Track.pdf) · [GSE279945](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945) · [task repo](https://github.com/openproblems-bio/task_perturbation_prediction)
- chemCPA: [Hetzel et al. 2022](https://www.arxiv.org/pdf/2204.13545) · github.com/theislab/chemCPA
- ChemBERTa, used by the competition's winning method (ref. 30 of the OP3 paper)
