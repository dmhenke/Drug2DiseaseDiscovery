# Chemical × single-cell perturbation model: design

**Goal (2026-11-02):** take a chemical's embedding, combine it with the embedded
expression profile of a single cell, run the pair through a model, and get a
generated, treated version of that cell. Then compare it with healthy cells.

Tags: **Known** (checked, with source) · **Inferred** · **Guess**.

## Recommendation

- **Chemical database: ChEMBL** (release 37, May 2026, about 2.92M compounds)
  as the master source. It provides structures, max clinical phase, mechanisms
  and targets.
  - Add the **Broad Drug Repurposing Hub** as the curated subset of about 5–7k
    compounds with mechanism-of-action annotations, linked to L1000/CMap by
    BRD ID.
  - Add **COCONUT** later, if natural products are wanted.
- **Join everything on InChIKey.** ChEMBL, PubChem, the Repurposing Hub, L1000,
  OP3 and Tahoe all carry it or can derive it.
- **Model: a chemCPA-style compositional autoencoder.** The cell latent is
  combined with a learned function of the molecule embedding and dose, then
  decoded to expression. It is trained on single-cell drug perturbations of
  **immune cells** (OP3: 146 compounds in PBMCs), optionally pretrained on L1000.

**The database is not the bottleneck.** A model can only give trustworthy outputs
for compounds that are in, or chemically close to, its **perturbation training
data**. A 2.9M-compound library makes 2.9M predictions, but almost all of them
are extrapolations. The perturbation dataset sets the ceiling (**Inferred**).

## Chemical sources

| Source | Size | What it adds | License | Role |
|---|---|---|---|---|
| **ChEMBL 37** | 2,921,148 compounds (**Known**, release notes) | SMILES/InChIKey, `max_phase`, mechanism and target tables | CC BY-SA 3.0 (**Known** for ChEMBL generally; check the release-37 notes) | Master structures and annotation |
| **Drug Repurposing Hub** (Broad) | 4,707 compounds in 2017; later 5,691 unique (**Known**) | Curated MoA, target and clinical phase; BRD IDs match L1000 | Not confirmed; check the portal | Interpretable subset; links to CMap/L1000 |
| PubChem | >100M | Name/CID resolution | Public domain | Lookup only, not a library |
| COCONUT | ~600k natural products (team notes; not checked) | Natural products | Not checked | Later. Mostly outside any training set |
| DrugBank | ~17k drugs (**Guess**) | Rich drug annotation | Restrictive licence for redistribution | Avoid in a public repo |

**Seed table built:** [`shared/compounds/seed_compounds.tsv`](../shared/compounds/seed_compounds.tsv),
from [`scripts/build_compound_table.py`](../scripts/build_compound_table.py).
- 129 compounds: 15 small-molecule SLE positive controls, plus the 114 PubChem-linked compounds from the L1000 reversal lookup.
- All 129 have a SMILES that RDKit parses; 99 are in ChEMBL; 28 have `max_phase` 4 (approved for some indication) (**Known**).
- The InChIKey join confirmed that L1000's TG101348 is fedratinib (**Known**).
- Note that `max_phase` is for any indication. Deucravacitinib shows 4 because it's approved for psoriasis, not SLE.

## How the generated cell would be produced

```
SMILES ──RDKit──▶ molecule embedding m (Morgan fingerprint 2048-bit, or a pretrained encoder)
dose   ─────────▶ d
cell x ──encoder─▶ z_basal            (expression → latent; UCE 1,280-d can be an input feature)
z_treated = z_basal + g(m, d)          (chemCPA: drug effect added in latent space)
x̂_treated = decoder(z_treated)        (negative-binomial expression for every gene)
```

**Comparison** for each compound c, applied to SLE classical monocytes:
- **Reversal score:** the fraction of the SLE signature (genes at padj < 0.05)
  that x̂ moves toward healthy.
- **Distance:** in latent space, from treated SLE cells to the healthy-monocyte
  centroid, against untreated SLE cells.
- **Controls:**
  - vehicle (DMSO) treated
  - positive-control drugs from the seed table
  - 100 random library compounds (the null distribution)

## Design constraints

1. **UCE has no decoder** (**Inferred**: UCE is an encoder; nothing in its repo
   maps embeddings back to genes). There are three ways to get a generated cell:
   - (a) **Stay in UCE space.** Predict the treated embedding and compare it with
     healthy embeddings. This is cheap, but the output has no genes, so "is IFI27
     down?" can't be answered.
   - (b) **Train a decoder from UCE to expression** on our own data, then use (a).
   - (c) **Use chemCPA's own expression autoencoder** and feed UCE in as an extra
     cell feature.

   **Recommend (c).** It gives gene-level outputs and reuses tested code.
2. **The training perturbations must be in immune cells.**

   | Dataset | Cells | Compounds | Fit for SLE |
   |---|---|---|---|
   | **OP3** (GSE279945) | PBMC from healthy donors: T, B, NK, myeloid (**Known**) | 146 (**Known**) | Best available. Same tissue type as our SLE data |
   | Tahoe-100M | ~100M cells, 50 cancer cell lines (**Known**, secondary sources) | ~1,100 (**Known**) | Cancer lines; too big for the node |
   | sciPlex3 | 3 cancer lines (**Guess**, from memory) | 188 (**Guess**) | chemCPA's benchmark; useful for checking our implementation |
   | LINCS L1000 | Bulk, ~978 measured genes, mostly epithelial lines | ~30k | Pretraining only (as chemCPA does) |
   | Kang 2018 (GSE96583) | Lupus PBMC ± IFN-β | 1 cytokine | Teaches the interferon axis |

3. **The interferon axis may be invisible in the training data.** OP3 PBMCs are
   cultured from healthy donors without interferon stimulation, so interferon
   genes are low at baseline. A JAK inhibitor then has little to suppress, and
   the model can't learn "JAK inhibitor → interferon genes down". Yet that is the
   main SLE axis (**Inferred**).
   - Mitigation: add IFN-stimulated PBMC data (Kang 2018) so the model sees the axis.
   - Check whether OP3 includes any JAK/TYK2 inhibitor or glucocorticoid. Its compound list hasn't been checked (**Guess** that it does).
4. **Generalizing to unseen compounds is the weak point.** chemCPA's main claim is
   better performance on unseen drugs, helped by bulk pretraining (**Known**,
   abstract). Arc's State model is described as "the first model to consistently
   beat simple linear baselines" (**Known**, Arc's own claim). So:
   - **Every evaluation needs a baseline:** the mean shift of the nearest training compound, or linear regression on fingerprints.
   - **Hold out whole compounds,** not cells, in the test split.
5. **Hardware.** chemCPA-scale models (MLP encoders and decoders) should train on
   the node's 8 GB GPU (**Inferred**, not tested). State and Tahoe-scale training
   would not fit. UCE-embedding OP3 at 0.93 cells/s costs about 30 h per 100k
   cells; OP3's cell count hasn't been checked.

## Stages (proposed; none built)

| Stage | Input → output | Check that can fail on its own |
|---|---|---|
| C0 Compound table | Names/IDs → `shared/compounds/*.tsv` with SMILES, InChIKey, ChEMBL ID, phase; RDKit-standardized | Every row parses; InChIKey round-trips; no duplicate InChIKey with conflicting annotation. **Seed done.** |
| C1 Molecule embedding | SMILES → fingerprints (`.npy`, keyed by InChIKey) | Tanimoto similarity of known analogue pairs > unrelated pairs (e.g. prednisolone vs methylprednisolone > prednisolone vs methotrexate). **Passed on the seed table** (Morgan r=2, 2048-bit): prednisolone–methylprednisolone 0.70, baricitinib–ruxolitinib 0.48, prednisolone–dexamethasone 0.42, against 0.09 (prednisolone–methotrexate) and 0.12 (baricitinib–hydroxychloroquine) |
| C2 Training data | OP3 (+ Kang 2018) → h5ad with `compound_inchikey`, `dose`, `donor`, `cell_type`, raw counts | Each compound's DMSO-normalized response in OP3 recovers known biology (glucocorticoids induce TSC22D3/FKBP5 in myeloid cells) |
| C3 Model | chemCPA, trained with whole compounds held out | Beats the linear and nearest-compound baselines on held-out compounds |
| C4 In-silico treatment | SLE monocytes × compound → generated cells | Positive controls rank above the random-compound null on reversal score |
| C5 Report | Scores → `docs/` page, extending the drug lookup | — |

## Open decisions

1. Accept OP3 as the training set? It's the only immune-cell, multi-compound single-cell set we've found.
2. Molecule encoder: start with Morgan fingerprints (no extra model, RDKit already installed) before trying pretrained encoders.
3. Run on the node, or on a larger GPU for Tahoe-scale pretraining?

## Sources

- chemCPA: Hetzel et al. 2022, [NeurIPS](https://proceedings.neurips.cc/paper_files/paper/2022/hash/aa933b5abc1be30baece1d230ec575a7-Abstract.html) · [arXiv 2204.13545](https://www.arxiv.org/pdf/2204.13545) · code github.com/theislab/chemCPA
- PRnet: Qi et al. 2024, [Nat Commun](https://ideas.repec.org/a/nat/natcom/v15y2024i1d10.1038_s41467-024-53457-1.html) (SMILES input not confirmed)
- TranSiGen: Tong et al. 2024, [PMC11199551](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC11199551/)
- Arc State / Tahoe-100M: [github.com/arcinstitute/state](https://github.com/arcinstitute/state) · [BiopharmaTrend](https://www.biopharmatrend.com/news/arc-institute-releases-its-first-virtual-cell-model-trained-on-270-million-cells-1302/)
- OP3: [NeurIPS 2024 D&B](https://proceedings.neurips.cc/paper_files/paper/2024/hash/24c4d51f3ef48dd2dbab78243ecb26a1-Abstract.html) · [GSE279945](https://www.omicsdi.org/dataset/geo/GSE279945)
- ChEMBL 37: [ChEMBL blog](https://chembl.blogspot.com/2026/05/chembl-37-is-here-bigger-cleaner-and.html) · [Cambridge MedChem](https://cambridgemedchemconsulting.com/2026/06/07/chembl-37-is-out/)
- Repurposing Hub: [Corsello et al. 2017](https://pmc.ncbi.nlm.nih.gov/articles/PMC5568558) · [NCI processed files with SMILES](https://computational.cancer.gov/dataset/drug-mechanism-action-information)
