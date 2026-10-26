#!/usr/bin/env Rscript
## 01_build_species_anndata.R  --  Stage 0 (organ-agnostic)
## ---------------------------------------------------------------------------
## Assemble a study's lung/kidney/... datasets into ONE AnnData PER SPECIES that
## meets UCE's input contract (raw counts, gene SYMBOLS, canonical metadata).
##
## Usage:
##   Rscript pipeline/01_build_species_anndata.R studies/lung
##   # or: STUDY_DIR=studies/lung Rscript pipeline/01_build_species_anndata.R
##
## Reads  <STUDY_DIR>/config.R and <STUDY_DIR>/data/*.h5ad
## Writes <STUDY_DIR>/interim/<STUDY_NAME>_<species>.h5ad
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(zellkonverter); library(S4Vectors)
})

## ---- resolve study dir + config ----
args <- commandArgs(trailingOnly = TRUE)
STUDY_DIR <- if (length(args) >= 1) args[1] else Sys.getenv("STUDY_DIR", unset = NA)
if (is.na(STUDY_DIR)) stop("Usage: Rscript 01_build_species_anndata.R <study_dir>")
STUDY_DIR <- normalizePath(STUDY_DIR, mustWork = TRUE)
source(file.path(STUDY_DIR, "config.R"))
INTERIM <- file.path(STUDY_DIR, "interim"); dir.create(INTERIM, showWarnings = FALSE, recursive = TRUE)
out_file <- function(sp) file.path(INTERIM, sprintf("%s_%s.h5ad", STUDY_NAME, sp))
abspath  <- function(p) if (grepl("^(/|~)", p)) p else file.path(STUDY_DIR, p)

looks_like_counts <- function(m) {
  s <- as.matrix(m[, seq_len(min(200, ncol(m))), drop = FALSE])
  if (any(s < 0)) return(FALSE)
  isTRUE(all.equal(s, round(s), check.attributes = FALSE))
}

load_one <- function(cfg) {
  message("== ", cfg$path, " (study=", cfg$study, ", species=", cfg$species, ") ==")
  sce <- readH5AD(abspath(cfg$path), use_hdf5 = FALSE)

  an <- cfg$counts_assay
  if (!is.null(an) && an %in% assayNames(sce)) assay(sce, "counts") <- assay(sce, an)
  if (!"counts" %in% assayNames(sce)) assay(sce, "counts") <- assay(sce, 1L)
  if (!looks_like_counts(assay(sce, "counts")))
    stop(cfg$path, ": 'counts' assay is not raw integer counts. UCE needs counts.")

  if (!is.null(cfg$symbols_from)) rownames(sce) <- as.character(rowData(sce)[[cfg$symbols_from]])
  if (mean(grepl("^ENS", rownames(sce), ignore.case = TRUE)) > 0.5)
    stop(cfg$path, ": rownames look like Ensembl IDs. UCE needs gene SYMBOLS.")
  rownames(sce) <- make.unique(rownames(sce))

  cd <- colData(sce)
  for (i in seq_along(cfg$obs_map)) {
    old <- names(cfg$obs_map)[i]; new <- cfg$obs_map[[i]]
    if (old %in% colnames(cd)) colnames(cd)[colnames(cd) == old] <- new
  }
  colData(sce) <- cd
  sce$study <- cfg$study; sce$species <- cfg$species
  for (col in c("sample_id", "condition"))
    if (!col %in% colnames(colData(sce)))
      stop(cfg$path, ": missing '", col, "' after obs_map. Have: ",
           paste(colnames(colData(sce)), collapse = ", "))

  cond <- CONDITION_MAP[tolower(trimws(as.character(sce$condition)))]
  keep <- !is.na(cond)
  if (any(!keep))
    message("   WARNING dropping ", sum(!keep), " cells with unmapped condition: ",
            paste(unique(as.character(sce$condition)[!keep]), collapse = ", "))
  sce <- sce[, keep]
  sce$condition <- factor(cond[keep], levels = c("healthy", "diseased"))
  sce$sample_id <- paste0(cfg$study, "::", as.character(sce$sample_id))

  keep_cols <- c("study", "species", "sample_id", "condition")
  sce2 <- SingleCellExperiment(
    assays = list(counts = assay(sce, "counts")),
    colData = colData(sce)[, keep_cols, drop = FALSE])
  rownames(sce2) <- rownames(sce)
  colnames(sce2) <- paste0(cfg$study, "::", colnames(sce2))  # barcodes collide across studies (B3)
  message("   -> ", ncol(sce2), " cells x ", nrow(sce2), " genes")
  sce2
}

merge_sces <- function(sces) {
  ## Intersect genes, not union: a gene a study never measured is missing, not zero (C1).
  ## Sparse throughout: no dense gene x cell matrix (B2).
  shared <- Reduce(intersect, lapply(sces, rownames))
  message("   genes shared across inputs: ", length(shared))
  do.call(cbind, lapply(sces, function(s) s[shared, ]))
}

qc_and_write <- function(sces, sp) {
  if (length(sces) == 0) return(invisible())
  merged <- merge_sces(sces)
  stopifnot(!anyDuplicated(colnames(merged)))  # unique cell IDs across studies (B3)
  ng <- Matrix::colSums(assay(merged, "counts") > 0)
  nc <- Matrix::rowSums(assay(merged, "counts") > 0)
  merged <- merged[nc >= MIN_CELLS_PER_GENE, ng >= MIN_GENES_PER_CELL]
  message("\n", toupper(sp), " cells per study x condition:")
  print(table(merged$study, merged$condition))
  writeH5AD(merged, out_file(sp), X_name = "counts")
  message("Wrote ", out_file(sp), ": ", ncol(merged), " cells x ", nrow(merged), " genes")
}

loaded <- lapply(INPUTS, function(cfg) list(sp = cfg$species, sce = load_one(cfg)))
for (sp in unique(vapply(loaded, `[[`, "", "sp")))
  qc_and_write(lapply(Filter(function(x) x$sp == sp, loaded), `[[`, "sce"), sp)
message("\nNext: embed each species file with UCE -> pipeline/run_uce.sh ", STUDY_DIR)
