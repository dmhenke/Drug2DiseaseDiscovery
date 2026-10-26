#!/usr/bin/env Rscript
## 03_pseudobulk_de.R
## ---------------------------------------------------------------------------
## Stage 3 of the UCE lung pipeline (R): differential expression on the
## homogenized master set, per species.
##
## Design rationale (see the guide):
##   * UCE homogenized the CELL GROUPING, not the counts -> DE on expression.
##   * DE PER CELL TYPE (defined in the shared UCE space) so we compare like
##     with like across studies (and, via joint clustering, across species).
##   * PSEUDOBULK: sum raw counts per (cell_type, sample_id). Unit of replication
##     is the donor/sample -> valid, non-inflated p-values.
##   * BATCH CORRECTION via the DESeq2 design: ~ study + condition. The study
##     term absorbs the batch effect UCE does NOT remove from expression; the
##     tested coefficient is condition (diseased vs healthy).
##   * Run separately per species; compare hits on 1:1 orthologs (hook below).
##
##   Rscript 03_pseudobulk_de.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(scuttle)
  library(DESeq2)
  library(S4Vectors)
})

## ===========================================================================
## CONFIG
## ===========================================================================
INPUTS <- c(human = "human_lung_annotated.rds",
            mouse = "mouse_lung_annotated.rds")   # from Stage 2
OUTDIR <- "de_results"
MIN_CELLS_PER_SAMPLE   <- 10   # a (cell_type,sample) pseudobulk needs >= this many cells
MIN_SAMPLES_PER_GROUP  <- 2    # need >=2 healthy AND >=2 diseased samples to test
MIN_COUNT              <- 10   # gene kept if >=MIN_COUNT in >=MIN_TOTAL_SAMPLES samples
MIN_TOTAL_SAMPLES      <- 2
ALPHA                  <- 0.05
## Optional cross-species ortholog comparison (needs 'babelgene').
COMPARE_ORTHOLOGS      <- TRUE
## ===========================================================================

## pseudobulk one species: sum counts per (sample_id, cell_type)
make_pseudobulk <- function(sce) {
  assay(sce, "counts") <- as.matrix(assay(sce, "counts"))
  pb <- aggregateAcrossCells(
    sce, ids = colData(sce)[, c("sample_id", "cell_type")],
    statistics = "sum", use.assay.type = "counts")
  pb   # SCE: columns = sample x cell_type; colData carries study/condition/species
}

de_one_celltype <- function(pb, ct, species) {
  sub <- pb[, pb$cell_type == ct & pb$ncells >= MIN_CELLS_PER_SAMPLE]
  cond_tab <- table(sub$condition)
  if ((cond_tab["healthy"] %||% 0) < MIN_SAMPLES_PER_GROUP ||
      (cond_tab["diseased"] %||% 0) < MIN_SAMPLES_PER_GROUP)
    return(list(res = NULL, note = sprintf(
      "skip [%s|%s]: too few samples/group (%s)", species, ct,
      paste(names(cond_tab), cond_tab, sep = "=", collapse = ", "))))

  ## batch-correction design + confound guard
  design <- ~ study + condition
  note <- NULL
  n_study <- length(unique(sub$study))
  both <- if (n_study > 1)
    sum(apply(table(sub$study, sub$condition) > 0, 1, all)) else 0
  if (n_study < 2) { design <- ~ condition }
  else if (both == 0) {
    design <- ~ condition
    note <- sprintf("WARNING [%s|%s]: study and condition fully confounded; ",
                    species, ct)
    note <- paste0(note, "study effect not separable -- interpret with caution.")
  }

  cts <- round(as.matrix(assay(sub, "counts")))
  mode(cts) <- "integer"
  cd <- as.data.frame(colData(sub))
  cd$condition <- factor(cd$condition, levels = c("healthy", "diseased"))
  cd$study <- factor(cd$study)

  dds <- DESeqDataSetFromMatrix(countData = cts, colData = cd, design = design)
  keep <- rowSums(counts(dds) >= MIN_COUNT) >= MIN_TOTAL_SAMPLES
  dds <- dds[keep, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("condition", "diseased", "healthy"),
                               alpha = ALPHA))
  res$gene <- rownames(res)
  res$cell_type <- ct; res$species <- species
  res$design <- paste(deparse(design), collapse = "")
  res <- res[order(res$padj), ]
  list(res = res, note = note)
}
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

run_species <- function(sce, species) {
  message("\n===== ", toupper(species), " : ", ncol(sce), " cells, ",
          length(unique(sce$cell_type)), " cell types =====")
  pb <- make_pseudobulk(sce)
  out <- list(); notes <- character()
  for (ct in sort(unique(pb$cell_type))) {
    r <- tryCatch(de_one_celltype(pb, ct, species),
                  error = function(e) list(res = NULL,
                    note = sprintf("ERROR [%s|%s]: %s", species, ct, conditionMessage(e))))
    if (!is.null(r$note)) { message(r$note); notes <- c(notes, r$note) }
    if (!is.null(r$res)) {
      nsig <- sum(r$res$padj < ALPHA, na.rm = TRUE)
      message(sprintf("  %s: %d DE genes at padj<%.2f", ct, nsig, ALPHA))
      write.csv(r$res, file.path(OUTDIR,
        sprintf("DE_%s_%s.csv", species, gsub("[/ ]", "_", ct))), row.names = FALSE)
      out[[ct]] <- r$res
    }
  }
  list(res = if (length(out)) do.call(rbind, out) else NULL, notes = notes)
}

## optional: intersect human vs mouse hits on 1:1 orthologs
compare_orthologs <- function(human_res, mouse_res) {
  if (!requireNamespace("babelgene", quietly = TRUE)) {
    message("(babelgene not installed; skipping ortholog comparison)"); return(invisible())
  }
  h <- human_res[which(human_res$padj < ALPHA), ]
  m <- mouse_res[which(mouse_res$padj < ALPHA), ]
  ort <- babelgene::orthologs(genes = unique(h$gene), species = "mouse", human = TRUE)
  h$mouse_gene <- ort$symbol[match(h$gene, ort$human_symbol)]
  merged <- merge(h, m, by.x = c("mouse_gene", "cell_type"),
                  by.y = c("gene", "cell_type"), suffixes = c("_human", "_mouse"))
  merged$concordant_sign <-
    sign(merged$log2FoldChange_human) == sign(merged$log2FoldChange_mouse)
  write.csv(merged, file.path(OUTDIR, "DE_conserved_human_mouse.csv"), row.names = FALSE)
  message(sprintf("\nConserved (ortholog) DE genes: %d (%d sign-concordant)",
                  nrow(merged), sum(merged$concordant_sign, na.rm = TRUE)))
}

main <- function() {
  dir.create(OUTDIR, showWarnings = FALSE)
  INPUTS <- INPUTS[file.exists(INPUTS)]
  results <- list(); all_notes <- character()
  for (sp in names(INPUTS)) {
    r <- run_species(readRDS(INPUTS[[sp]]), sp)
    results[[sp]] <- r$res; all_notes <- c(all_notes, r$notes)
    if (!is.null(r$res)) write.csv(r$res,
      file.path(OUTDIR, sprintf("DE_all_%s.csv", sp)), row.names = FALSE)
  }
  if (COMPARE_ORTHOLOGS && all(c("human", "mouse") %in% names(results)) &&
      !is.null(results$human) && !is.null(results$mouse))
    compare_orthologs(results$human, results$mouse)
  writeLines(all_notes, file.path(OUTDIR, "warnings.txt"))
  message("\nResults + warnings written to ", OUTDIR, "/")
}

main()
