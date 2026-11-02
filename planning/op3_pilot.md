# OP3 pilot and full run

**Decided 2026-11-02:** train on OP3; start with Morgan fingerprints.
Design context: [`chemical_perturbation_model.md`](chemical_perturbation_model.md).

Two stages. **Stage P** is a pilot that can be run in an afternoon on a laptop and
can fail cheaply. **Stage F** is the full run: pretrained molecule encoders, L1000
pretraining, single-cell resolution, on the node's GPU.

**Run Stage F only if Stage P passes.** If a fingerprint model cannot beat the
baselines on OP3's own held-out compounds, a larger model on the same data will
not rescue it; the problem is the data or the task (**Inferred**).

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

**Already confirmed overlap with our own work:** 5 of our 129 seed compounds appear in the OP3 test split — prednisolone, ruxolitinib, vorinostat, IMD-0354, TPCA-1 — and three of those (vorinostat, IMD-0354, TPCA-1) came out of the L1000 monocyte reversal independently. The full 140-compound list will add more; it needs `de_train.h5ad`.

## Stage P — pilot: fingerprints → pseudobulk response

**Question:** given a compound's structure and a cell type, can we predict its
transcriptional effect well enough to be worth building on?

**Data:** `de_train.h5ad` / `de_test.h5ad` (292 MB). These are limma DE statistics per
(compound, cell type) against DMSO, 5,317 genes. No single cells, no GPU.

**Model:** ridge regression, fingerprint (2048) + cell type (one-hot) → 5,317 DE values.
One `sklearn` call. This is deliberately the dumbest thing that uses structure.

**Splits:** hold out **whole compounds**, never rows. OP3's own split also holds out
whole cell types (B and Myeloid), which is the harder benchmark task; report both.

**Baselines it must beat.** Without these the result is meaningless:

| Baseline | What it tests |
|---|---|
| Global mean DE across training compounds | Is there any compound-specific signal at all? |
| Cell-type mean DE | Is the structure doing anything the cell type doesn't? |
| **Nearest training compound** by Tanimoto | Is the model learning chemistry, or just copying the closest neighbour? |
| Random fingerprint permutation | The null. |

**Pass criteria, fixed before running:**

1. The ridge model beats the global and cell-type means on held-out compounds, by mean correlation across compounds.
2. It beats the nearest-training-compound baseline. *(If it does not, say so and use the nearest-neighbour method; it is simpler and would be the honest answer.)*
3. **Biology check:** for held-out prednisolone in myeloid cells, predicted TSC22D3 and FKBP5 are up. These are canonical glucocorticoid targets, so a model that misses them has not learned the strongest signal in the data.
4. **Bootstrap CIs** over held-out compounds, not point estimates.

**Cost:** about 30 min of compute on the Mac (**Guess**; ridge on 560×2048 is seconds, the h5ad reads dominate).

## Stage F — full run

Three changes, each testable on its own against Stage P:

| Change | What it buys | Cost |
|---|---|---|
| **F1. Pretrained molecule encoders** instead of fingerprints: ChemBERTa (the NeurIPS competition winner used it), plus RDKit descriptors and a GNN | Chemistry the fingerprint cannot express. Compare on the same splits | Encoder inference only; minutes on GPU |
| **F2. L1000 pretraining, then fine-tune on OP3** — chemCPA's central claim | ~30k compounds of bulk signal to generalize from. Needs the L1000 level-5 matrix (several GB) and gene-space alignment with OP3's 5,317 | Hours on the P4000 (**Guess**) |
| **F3. Single-cell chemCPA** on `sc_train.h5ad` (893 MB, ~300k cells): latent `z_treated = z_basal + g(m, dose)`, negative-binomial decoder | Per-cell predictions, which is what the SLE application needs | MLP autoencoder; should fit in 8 GB (**Guess**, untested) |

**Each must beat Stage P on the same held-out compounds** to be kept. Stage P's
ridge model is the control for the whole of Stage F, and the paper to beat for
F2 is chemCPA's own claim that bulk pretraining helps on unseen drugs.

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
2. **Nearest-neighbour may win.** Common in this field. It would be a real result and should be reported, not buried.
3. **Cross-cell-type generalization is the published hard task.** Expect worse numbers on OP3's own split than on a plain compound holdout.
4. **Healthy PBMC → SLE monocyte transfer is unvalidated** and probably cannot be validated without treated SLE cells.
5. **Three donors** means donor effects cannot be separated from compound effects well. Keep donor as a covariate and never as a split.

## Commands

```bash
scripts/fetch_op3.sh                         # ~1.1 GB, resumable, size-checked
python scripts/check_op3.py shared/perturbation/op3
python scripts/mol_fingerprints.py --self-check
```

`scripts/mol_fingerprints.py` already passes its analogue check (prednisolone–methyl­prednisolone
0.70 and baricitinib–ruxolitinib 0.48, against 0.09 and 0.12 for unrelated pairs) and has
been run on the 129-compound seed table: 2048-bit vectors, InChIKeys round-tripping 100 %
against PubChem's (**Known**).

Not yet written: the Stage P training script. It needs the downloaded files first,
so that the splits and baselines are written against the real columns.

## Sources

- OP3: [NeurIPS 2024 D&B paper](https://papers.neurips.cc/paper_files/paper/2024/file/24c4d51f3ef48dd2dbab78243ecb26a1-Paper-Datasets_and_Benchmarks_Track.pdf) · [GSE279945](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE279945) · [task repo](https://github.com/openproblems-bio/task_perturbation_prediction)
- chemCPA: [Hetzel et al. 2022](https://www.arxiv.org/pdf/2204.13545) · github.com/theislab/chemCPA
- ChemBERTa, used by the competition's winning method (ref. 30 of the OP3 paper)
