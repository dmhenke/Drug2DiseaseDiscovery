# SLE bulk whole-blood meta-analysis

Imported 2026-11-02 from `github.com/lkmklsmn/sle_analysis`, commit `ab8fc77`
(2025-04-18), received as a zip. `code/` and `outputs/merged_sle_diff_expr.csv`
are copied unchanged. The `.git`, `.Rproj.user`, `.Rhistory` and the MSigDB `.gmt`
files were left out; re-download the `.gmt` files from MSigDB (v2024.1, Hallmark
`h.all` and immunologic `c7.all`).

## What it does

`code/run_de_analysis_across_cohorts.R`:
1. Downloads four microarray studies from GEO: GSE49454 (Chiche), GSE65391 (Banchereau), GSE45291 and GSE121239 (Petri).
2. Collapses probes to genes.
3. Fits limma `~ study + condition`, SLE vs control, after dropping the rheumatoid arthritis samples.

`code/GSEA_test.R` is draft fgsea code for the Hallmark and C7 gene sets.

**Result as delivered:** 5,794 genes, 2,857 at adj.P < 0.05. The top hit is IFI27 at log2FC 4.59 (adj.P 6e-79).

## Problems found in review (2026-11-02)

The code was read but **not rerun**: that needs the GEO downloads and Bioconductor.
Problems 1 and 2 can be seen in the delivered output.

1. **The gene intersection counts probes, not studies.**
   - The line `allgenes <- names(which(table(c(chiche_gene, banchereau_gene, petri_1_gene, petri_2_gene)) == 4))` counts probe rows. A gene with two probes on one array is excluded, even if it's on all four arrays.
   - Consequence: **most of the interferon signature is missing** from the output. IFI44L, IFI6, IFIT1, IFIT3, RSAD2, OAS1, STAT1, SIGLEC1 and TNFSF13B are all absent; ISG15, MX1 and LY6E survive.
   - The same line lets in genes absent from some studies: four probes on one array also add up to 4.
   - Fix: `allgenes <- Reduce(intersect, list(rownames(chiche_expr), rownames(banchereau_expr), rownames(petri_1_expr), rownames(petri_2_expr)))`.
2. **700 output rows are named `NA..*`.** These are genes missing from the first study (Chiche): `match()` returns NA rows, and the row names come from that study. The second-strongest "gene" in the output (`NA..638`, log2FC −2.85) is one of these. This is a consequence of problem 1.
3. **Repeated samples from the same patients** (pseudoreplication). GSE65391 and GSE121239 are, to my recollection, longitudinal: several visits per patient (**Guess**, check the GEO sample sheets). limma treats every visit as an independent person, which inflates significance. Fix: keep one visit per patient, or use `duplicateCorrelation(block = patient)`.
4. **`merged_expr[, -other]` deletes every column if `other` is empty**, i.e. if no rheumatoid arthritis samples are present. Fix: `if (length(other)) …`.
5. **Probe ties are broken with an unseeded `sample()`.** Add `set.seed()`.
6. **The condition labels have not been checked.** Only the `merged_treatSLE` coefficient is reported. Any label other than exactly `SLE` or `CTRL` becomes its own level and is silently excluded from the contrast. Check with `table(study, merged_treat)`.
7. **`GSEA_test.R` doesn't run:** there's a stray `}`, and `read.gmt` needs clusterProfiler. It also ranks genes by logFC; the moderated t (`top_table$t`) is the usual choice, because it accounts for variance.
8. **Whole blood vs PBMC.** These arrays include neutrophils; the single-cell data (PBMC) doesn't. Differences in which cells are present affect the two datasets differently.

## Agreement with the single-cell signature

`compare_with_scrna.py` compares this table with the pseudobulk DE from the lupus UCE test run (`studies/lupus_test`, commit `bc1cb88`). The output is `outputs/scrna_concordance_lupus_test_bc1cb88.tsv`.

| Cell type (scRNA) | Significant genes | …of which on the bulk table | Replicated in bulk (same sign, adj.P < 0.05) |
|---|---|---|---|
| Classical monocyte | 68 | 20 | 17 |
| Non-classical monocyte | 35 | 12 | 10 |
| CD4 T | 46 | 9 | 8 |
| NK | 38 | 12 | 4 |

- **Where a single-cell gene is on the bulk table, it mostly replicates in monocytes and CD4 T cells** (**Known**). Most single-cell hits aren't on the bulk table at all, largely because of problem 1.
- **Ribosomal genes:** 26 of the 31 RPL/RPS genes on the bulk table are significantly down in SLE whole blood. So the single-cell ribosomal signal is not peculiar to one dataset. This fits a change in which cells are present, such as fewer lymphocytes (**Inferred**). It does not separate that from sample handling.
- **FOS** is not changed in bulk (log2FC 0.07, adj.P 0.32), unlike in single-cell T and NK cells.

Rerun the comparison after fixing problem 1. The interferon genes should then appear in the bulk table and replicate.
