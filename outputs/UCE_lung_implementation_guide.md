# Implementing UCE for cross-study, cross-species integration and DE of lung scRNA-seq (R workflow)

**Paper:** Rosen, Roohani, Agrawal *et al.* "Universal cell embedding provides a foundation model for cell biology." *Nature* (2026). DOI 10.1038/s41586-026-10689-z.
**Code:** https://github.com/snap-stanford/UCE
**Model weights (33-layer, the paper model):** https://figshare.com/articles/dataset/Universal_Cell_Embedding_Model_Files/24320806

---

## 1. What UCE does — and the two constraints that shape this workflow

UCE is a **zero-shot transformer foundation model**. Feed it raw scRNA-seq counts for a cell and it returns a fixed **1,280-dimensional embedding**. It puts data from different studies, platforms, gene panels, and **species** onto one "playing field" because it tokenizes genes with **ESM2 protein-language embeddings** (the amino-acid sequence of each gene's protein), not by gene identity. No training, fine-tuning, or gene harmonization is required, and — critically for you — mouse and human cells land in the *same* universal space.

Two hard constraints determine how we build the pipeline:

**Constraint 1 — the embedding step is Python-only.** UCE ships as a PyTorch package; there is no R port. So Stage 1 (producing the embedding) is run as the authors' command-line script. It is a single self-contained command per dataset and needs no R interop. **Everything else — assembling the input, joint clustering/annotation, and the differential expression — is done in R**, which is what the rest of this guide covers.

**Constraint 2 — UCE homogenizes the cell *grouping*, not the counts.** The output is a per-cell embedding (`obsm["X_uce"]`), not a corrected gene-by-cell matrix, and the paper explicitly says UCE **does not** replace batch correction. Its 1,280 dimensions are latent features, not genes, so **gene-level DE cannot be run on the embedding**. UCE's job is to define cross-study/cross-species-consistent cell types; the DE then runs on expression counts with **`study` as a covariate in the DESeq2 design** — the batch-correction approach you asked to keep.

---

## 2. How mouse + human changes the structure

UCE embeds each dataset under a single `--species` flag, so **you cannot mix species inside one embedding run.** The workflow therefore splits by species for embedding, then re-joins in the shared space:

```
 human h5ads ─┐                              ┌─ embed (--species human) ─┐
              ├─ R: build 2 species files ───┤                           ├─ R: stack X_uce ─→ joint
 mouse h5ads ─┘   (counts, symbols, meta)    └─ embed (--species mouse) ─┘   embeddings         clustering
                                                                                                    │
                                                                              per-species pseudobulk DE
                                                                              DESeq2  ~ study + condition
```

- **Joint clustering is cross-species.** Because both species occupy the same UCE space, you stack the two embedding matrices and cluster once. Homologous cell types (human alveolar macrophage ↔ mouse alveolar macrophage) co-cluster, giving you one consistent `cell_type` label set across the whole master collection. This is the payoff of UCE's universality.
- **DE stays within species.** Genes are not identical across species, and dispersion differs, so the negative-binomial test is run **separately for mouse and for human**, each with `~ study + condition`. You then compare the two hit lists on one-to-one orthologs (script provides the mapping hook). Pooling species into one DESeq2 model is only advisable after collapsing to shared orthologs with `species` added as a further covariate — offered as an option, but per-species-then-intersect is the cleaner default.

---

## 3. Stage-by-stage

**Stage 0 — Build one AnnData per species (`01_build_species_anndata.R`).**
For each lung dataset, ensure raw integer counts, gene **symbols** as feature names (human symbols for human, MGI symbols for mouse — *not* Ensembl IDs), and per-cell metadata `study`, `sample_id` (donor/replicate), `condition` (`healthy`/`diseased`), `species`. The script assembles all human datasets into `human_lung.h5ad` and all mouse datasets into `mouse_lung.h5ad` using `zellkonverter::writeH5AD` (writes the counts assay into `.X`, which is what UCE reads).

**Stage 1 — Embed each species with UCE (`run_uce.sh`).** Two runs, one per species file, each with the matching `--species`. Output: two h5ads with `obsm["X_uce"]` added. Needs a GPU; if you have none locally, run just this step on a cloud GPU and bring the embedded files back — all R stages are CPU-only.

**Stage 2 — Joint integration & annotation (`02_integrate_and_cluster.R`).** Read both embedded h5ads with `zellkonverter::readH5AD` (the `X_uce` embedding arrives in `reducedDims`). Stack the two embedding matrices, build a shared-nearest-neighbour graph on the combined embedding (`bluster`), cluster (Louvain/Leiden via `igraph`), and — recommended — transfer labels from an annotated reference lung atlas embedded with the same UCE model. Attach the resulting `cell_type` back onto each per-species object for DE. Optional `scran`/silhouette diagnostics confirm studies and species are genuinely mixed.

**Stage 3 — Pseudobulk DE with DESeq2 (`03_pseudobulk_de.R`).** For each species, sum raw counts per `(cell_type, sample_id)` with `scuttle::aggregateAcrossCells`, then run `DESeq2` per cell type with **`design = ~ study + condition`**. Pseudobulk-per-sample respects the true replication unit (the donor, not the cell) and avoids the false-positive inflation of cell-level tests; the `study` term is the batch correction. Results are written per cell type and combined, with guards for confounded designs.

---

## 4. Why this design

- **Why `~ study + condition` rather than trusting UCE to remove batch?** The authors state UCE aligns cell *types* but does **not** strip known experimental confounders from expression. Modeling `study` in DESeq2 is what yields valid p-values. This is the batch-correction step you asked to retain.
- **Why pseudobulk, not `Seurat::FindMarkers`/Wilcoxon per cell?** Cells from one donor are correlated; treating each as an independent replicate massively inflates false positives on multi-sample designs. Pseudobulk collapses to one profile per sample and tests across samples — the biologically correct level.
- **Why per-species DE?** Orthology is not identity; mouse and human genes differ in sequence, expression scale, and dispersion. UCE unifies the cell *grouping* across species, but the gene-level test belongs within a species (then compared on orthologs).
- **Design requirement:** ideally multiple studies contribute **both** conditions within each species. If `study` and `condition` are perfectly confounded (all diseased from one study, all healthy from another), no covariate can separate disease from batch — the script flags and, where unavoidable, drops the `study` term with a warning.

---

## 5. Environment

```bash
# --- Embedding (Python, GPU box) ---
git clone https://github.com/snap-stanford/UCE.git && cd UCE
pip install -r requirements.txt
# download 33-layer weights from figshare (README link) into model_files/
```

```r
# --- Analysis (R, CPU) ---  Bioconductor + CRAN
install.packages("BiocManager")
BiocManager::install(c(
  "zellkonverter",   # AnnData <-> SingleCellExperiment (h5ad I/O)
  "SingleCellExperiment", "scuttle", "scran", "bluster", "DESeq2",
  "SummarizedExperiment", "S4Vectors"
))
install.packages("igraph")   # graph clustering
# optional, for cross-species ortholog mapping in Stage 3:
BiocManager::install("babelgene")   # or biomaRt
```

R ≥ 4.3, Bioconductor ≥ 3.18. `zellkonverter` uses a managed Python (basilisk) under the hood purely for h5ad reading/writing — you do not interact with it.

---

## 6. Scripts in this bundle

| File | Stage | Language | What it does |
|------|-------|----------|--------------|
| `01_build_species_anndata.R` | 0 | R | Merge lung datasets into `human_lung.h5ad` + `mouse_lung.h5ad`; enforce counts, symbols, metadata. |
| `run_uce.sh` | 1 | shell→Python | Two UCE runs (`--species human`, `--species mouse`). The only Python step. |
| `02_integrate_and_cluster.R` | 2 | R | Stack `X_uce` from both species, SNN graph + Louvain, optional reference label transfer, mixing diagnostics. |
| `03_pseudobulk_de.R` | 3 | R | `aggregateAcrossCells` pseudobulk + DESeq2 `~ study + condition`, per species/cell type; ortholog comparison hook. |

Edit the CONFIG block at the top of each script (paths + the metadata column names), then run in order. Note: the earlier all-Python scripts are superseded by these R versions; only `run_uce.sh` (the embedding step) remains shared.
