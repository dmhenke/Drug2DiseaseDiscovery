#!/usr/bin/env Rscript
## 03_pseudobulk_de.R  --  Stage 3 (organ-agnostic)
## ---------------------------------------------------------------------------
## Pseudobulk differential expression per species, per cell type, with the
## DESeq2 batch-correction design  ~ study + condition  (diseased vs healthy).
## UCE homogenized the cell GROUPING; DE runs on expression with study as the
## batch covariate. Cross-species hits are compared on 1:1 orthologs.
##
## Usage:  Rscript pipeline/03_pseudobulk_de.R studies/lung
## Reads   <STUDY_DIR>/interim/<STUDY_NAME>_<species>_annotated.rds
## Writes  <STUDY_DIR>/results/DE_*.csv
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(scuttle); library(DESeq2); library(S4Vectors)
})
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
STUDY_DIR <- if (length(args) >= 1) args[1] else Sys.getenv("STUDY_DIR", unset = NA)
if (is.na(STUDY_DIR)) stop("Usage: Rscript 03_pseudobulk_de.R <study_dir>")
STUDY_DIR <- normalizePath(STUDY_DIR, mustWork = TRUE)
source(file.path(STUDY_DIR, "config.R"))
INTERIM <- file.path(STUDY_DIR, "interim")
RESULTS <- file.path(STUDY_DIR, "results"); dir.create(RESULTS, showWarnings = FALSE, recursive = TRUE)
ann_file <- function(sp) file.path(INTERIM, sprintf("%s_%s_annotated.rds", STUDY_NAME, sp))

make_pseudobulk <- function(sce) {
  ## No as.matrix(): keep the sparse counts (B2). aggregateAcrossCells accepts dgCMatrix.
  aggregateAcrossCells(sce, ids = colData(sce)[, c("sample_id", "cell_type")],
                       statistics = "sum", use.assay.type = "counts")
}

de_one <- function(pb, ct, sp) {
  sub <- pb[, pb$cell_type == ct & pb$ncells >= MIN_CELLS_PER_SAMPLE]
  ctab <- table(sub$condition)
  if ((ctab["healthy"] %||% 0) < MIN_SAMPLES_PER_GROUP ||
      (ctab["diseased"] %||% 0) < MIN_SAMPLES_PER_GROUP)
    return(list(res = NULL, note = sprintf("skip [%s|%s]: too few samples/group (%s)",
      sp, ct, paste(names(ctab), ctab, sep = "=", collapse = ", "))))

  design <- ~ study + condition; note <- NULL
  n_study <- length(unique(sub$study))
  both <- if (n_study > 1) sum(apply(table(sub$study, sub$condition) > 0, 1, all)) else 0
  if (n_study < 2) design <- ~ condition
  else if (both == 0) { design <- ~ condition
    note <- sprintf("WARNING [%s|%s]: study & condition fully confounded; study effect not separable.", sp, ct) }

  cts <- round(as.matrix(assay(sub, "counts"))); mode(cts) <- "integer"
  cd <- as.data.frame(colData(sub))
  cd$condition <- factor(cd$condition, levels = c("healthy", "diseased"))
  cd$study <- factor(cd$study)
  dds <- DESeqDataSetFromMatrix(cts, cd, design)
  dds <- dds[rowSums(counts(dds) >= MIN_COUNT) >= MIN_TOTAL_SAMPLES, ]
  dds <- DESeq(dds, quiet = TRUE)
  res <- as.data.frame(results(dds, contrast = c("condition", "diseased", "healthy"), alpha = ALPHA))
  res$gene <- rownames(res); res$cell_type <- ct; res$species <- sp
  res$design <- paste(deparse(design), collapse = "")
  list(res = res[order(res$padj), ], note = note)
}

run_species <- function(sce, sp) {
  message("\n===== ", toupper(sp), " : ", ncol(sce), " cells, ",
          length(unique(sce$cell_type)), " cell types =====")
  pb <- make_pseudobulk(sce); out <- list(); notes <- character()
  for (ct in sort(unique(pb$cell_type))) {
    r <- tryCatch(de_one(pb, ct, sp), error = function(e)
      list(res = NULL, note = sprintf("ERROR [%s|%s]: %s", sp, ct, conditionMessage(e))))
    if (!is.null(r$note)) { message(r$note); notes <- c(notes, r$note) }
    if (!is.null(r$res)) {
      message(sprintf("  %s: %d DE genes at padj<%.2f", ct, sum(r$res$padj < ALPHA, na.rm = TRUE), ALPHA))
      write.csv(r$res, file.path(RESULTS, sprintf("DE_%s_%s.csv", sp, gsub("[/ ]", "_", ct))), row.names = FALSE)
      out[[ct]] <- r$res
    }
  }
  list(res = if (length(out)) do.call(rbind, out) else NULL, notes = notes)
}

compare_orthologs <- function(h, m) {
  if (!requireNamespace("babelgene", quietly = TRUE)) {
    message("(babelgene not installed; skipping ortholog comparison)"); return(invisible()) }
  h <- h[which(h$padj < ALPHA), ]; m <- m[which(m$padj < ALPHA), ]
  ort <- babelgene::orthologs(genes = unique(h$gene), species = "mouse", human = TRUE)
  h$mouse_gene <- ort$symbol[match(h$gene, ort$human_symbol)]
  merged <- merge(h, m, by.x = c("mouse_gene", "cell_type"), by.y = c("gene", "cell_type"),
                  suffixes = c("_human", "_mouse"))
  merged$concordant_sign <- sign(merged$log2FoldChange_human) == sign(merged$log2FoldChange_mouse)
  write.csv(merged, file.path(RESULTS, "DE_conserved_human_mouse.csv"), row.names = FALSE)
  message(sprintf("\nConserved (ortholog) DE genes: %d (%d sign-concordant)",
                  nrow(merged), sum(merged$concordant_sign, na.rm = TRUE)))
}

species <- c("human", "mouse")[file.exists(c(ann_file("human"), ann_file("mouse")))]
results <- list(); all_notes <- character()
for (sp in species) {
  r <- run_species(readRDS(ann_file(sp)), sp); results[[sp]] <- r$res; all_notes <- c(all_notes, r$notes)
  if (!is.null(r$res)) write.csv(r$res, file.path(RESULTS, sprintf("DE_all_%s.csv", sp)), row.names = FALSE)
}
if (isTRUE(COMPARE_ORTHOLOGS) && all(c("human", "mouse") %in% names(results)) &&
    !is.null(results$human) && !is.null(results$mouse))
  compare_orthologs(results$human, results$mouse)
writeLines(all_notes, file.path(RESULTS, "warnings.txt"))
message("\nResults written to ", RESULTS, "/")
