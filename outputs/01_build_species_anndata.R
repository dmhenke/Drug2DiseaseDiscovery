#!/usr/bin/env Rscript
## 01_build_species_anndata.R
## ---------------------------------------------------------------------------
## Stage 0 of the UCE lung pipeline (R).
##
## Assemble lung scRNA-seq datasets into ONE AnnData PER SPECIES, each meeting
## UCE's input contract:
##   * assay "counts"  = raw integer counts   (UCE samples genes by count)
##   * rownames        = GENE SYMBOLS         (human or MGI; NOT Ensembl IDs)
##   * colData$study     study/dataset of origin   (DESeq2 covariate)
##   * colData$sample_id biological replicate/donor (pseudobulk unit)
##   * colData$condition 'healthy' or 'diseased'    (the contrast)
##   * colData$species   'human' or 'mouse'
##
## UCE embeds one species per run, so we emit human_lung.h5ad and mouse_lung.h5ad
## separately. Written with zellkonverter (counts assay -> AnnData .X).
##
##   Rscript 01_build_species_anndata.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter)
  library(S4Vectors)
})

## ===========================================================================
## CONFIG -- edit for your data
## Each entry describes one input .h5ad. `obs_map` renames colData columns to the
## canonical names; `symbols_from` names a rowData column holding gene symbols if
## the rownames are Ensembl IDs (set NULL if rownames are already symbols).
## ===========================================================================
INPUTS <- list(
  list(path = "data/studyA_human_lung.h5ad", study = "studyA", species = "human",
       obs_map = c(donor = "sample_id", disease_state = "condition"),
       symbols_from = NULL,          counts_assay = "counts"),
  list(path = "data/studyB_human_lung.h5ad", study = "studyB", species = "human",
       obs_map = c(patient = "sample_id", status = "condition"),
       symbols_from = "feature_name", counts_assay = "counts"),
  list(path = "data/studyC_mouse_lung.h5ad", study = "studyC", species = "mouse",
       obs_map = c(mouse_id = "sample_id", group = "condition"),
       symbols_from = "symbol",      counts_assay = "counts")
  ## ... add every lung study (human and mouse) here ...
)

## Map raw condition labels onto exactly {'healthy','diseased'}. Unlisted -> dropped.
CONDITION_MAP <- c(healthy = "healthy", normal = "healthy", control = "healthy",
                   ctrl = "healthy",
                   diseased = "diseased", disease = "diseased", ipf = "diseased",
                   copd = "diseased", fibrosis = "diseased", tumor = "diseased")

OUT_HUMAN <- "human_lung.h5ad"
OUT_MOUSE <- "mouse_lung.h5ad"
MIN_GENES_PER_CELL <- 200   # paper CxG preprocessing
MIN_CELLS_PER_GENE <- 10
## ===========================================================================

looks_like_counts <- function(m) {
  s <- m[, seq_len(min(200, ncol(m))), drop = FALSE]
  s <- as.matrix(s)
  if (any(s < 0)) return(FALSE)
  isTRUE(all.equal(s, round(s), check.attributes = FALSE))
}

load_one <- function(cfg) {
  message("== ", cfg$path, " (study=", cfg$study, ", species=", cfg$species, ") ==")
  sce <- readH5AD(cfg$path, use_hdf5 = FALSE)

  ## 1) counts assay
  an <- cfg$counts_assay
  if (!is.null(an) && an %in% assayNames(sce)) assay(sce, "counts") <- assay(sce, an)
  if (!"counts" %in% assayNames(sce)) assay(sce, "counts") <- assay(sce, 1L)
  if (!looks_like_counts(assay(sce, "counts")))
    stop(cfg$path, ": 'counts' assay is not raw integer counts. UCE needs counts; ",
         "point counts_assay at the raw layer. Do NOT feed log-normalized data.")

  ## 2) gene symbols as rownames
  if (!is.null(cfg$symbols_from)) rownames(sce) <- as.character(rowData(sce)[[cfg$symbols_from]])
  if (mean(grepl("^ENS", rownames(sce), ignore.case = TRUE)) > 0.5)
    stop(cfg$path, ": rownames look like Ensembl IDs. UCE needs gene SYMBOLS; ",
         "set symbols_from to the rowData column with symbols.")
  rownames(sce) <- make.unique(rownames(sce))

  ## 3) canonical colData
  cd <- colData(sce)
  for (i in seq_along(cfg$obs_map)) {
    old <- names(cfg$obs_map)[i]; new <- cfg$obs_map[[i]]
    if (old %in% colnames(cd)) colnames(cd)[colnames(cd) == old] <- new
  }
  colData(sce) <- cd
  sce$study   <- cfg$study
  sce$species <- cfg$species
  for (col in c("sample_id", "condition"))
    if (!col %in% colnames(colData(sce)))
      stop(cfg$path, ": missing '", col, "' after obs_map. Have: ",
           paste(colnames(colData(sce)), collapse = ", "))

  ## normalize condition, drop unmapped
  cond <- CONDITION_MAP[tolower(trimws(as.character(sce$condition)))]
  keep <- !is.na(cond)
  if (any(!keep))
    message("   WARNING dropping ", sum(!keep), " cells with unmapped condition: ",
            paste(unique(as.character(sce$condition)[!keep]), collapse = ", "))
  sce <- sce[, keep]
  sce$condition <- factor(cond[keep], levels = c("healthy", "diseased"))
  ## globally-unique sample IDs
  sce$sample_id <- paste0(cfg$study, "::", as.character(sce$sample_id))

  ## keep only counts + minimal colData to avoid merge conflicts
  keep_cols <- c("study", "species", "sample_id", "condition")
  sce2 <- SingleCellExperiment(
    assays = list(counts = assay(sce, "counts")),
    colData = colData(sce)[, keep_cols, drop = FALSE]
  )
  rownames(sce2) <- rownames(sce)
  message("   -> ", ncol(sce2), " cells x ", nrow(sce2), " genes")
  sce2
}

## union-of-genes merge for a set of SCEs (genes differ across studies; missing = 0)
merge_sces <- function(sces) {
  all_genes <- Reduce(union, lapply(sces, rownames))
  pad <- function(s) {
    m <- matrix(0, nrow = length(all_genes), ncol = ncol(s),
                dimnames = list(all_genes, colnames(s)))
    m[rownames(s), ] <- as.matrix(assay(s, "counts"))
    SingleCellExperiment(assays = list(counts = as(m, "dgCMatrix")), colData = colData(s))
  }
  do.call(cbind, lapply(sces, pad))
}

qc_and_write <- function(sces, outfile, label) {
  if (length(sces) == 0) { message("No ", label, " datasets; skipping ", outfile); return(invisible()) }
  merged <- merge_sces(sces)
  ## light QC matching the paper's CxG preprocessing (no HVG subsetting!)
  ngenes <- Matrix::colSums(assay(merged, "counts") > 0)
  ncells <- Matrix::rowSums(assay(merged, "counts") > 0)
  merged <- merged[ncells >= MIN_CELLS_PER_GENE, ngenes >= MIN_GENES_PER_CELL]
  message("\n", label, " cells per study x condition:")
  print(table(merged$study, merged$condition))
  writeH5AD(merged, outfile, X_name = "counts")   # counts assay -> AnnData .X
  message("Wrote ", outfile, ": ", ncol(merged), " cells x ", nrow(merged), " genes")
}

main <- function() {
  loaded <- lapply(INPUTS, function(cfg) list(sp = cfg$species, sce = load_one(cfg)))
  human <- lapply(Filter(function(x) x$sp == "human", loaded), `[[`, "sce")
  mouse <- lapply(Filter(function(x) x$sp == "mouse", loaded), `[[`, "sce")
  qc_and_write(human, OUT_HUMAN, "HUMAN")
  qc_and_write(mouse, OUT_MOUSE, "MOUSE")
  message("\nNext: embed each species file with UCE (see run_uce.sh).")
}

main()
