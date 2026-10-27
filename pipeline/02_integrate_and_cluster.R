#!/usr/bin/env Rscript
## 02_integrate_and_cluster.R  --  Stage 2 (organ-agnostic)
## ---------------------------------------------------------------------------
## Jointly cluster a study's UCE-embedded species files in the shared UCE space
## (mouse + human co-embed there), assign one consistent cell_type label set,
## and write per-species SCEs (with cell_type) for DE.
##
## Usage:  Rscript pipeline/02_integrate_and_cluster.R studies/lung
## Reads   <STUDY_DIR>/interim/uce/<STUDY_NAME>_<species>.h5ad
## Writes  <STUDY_DIR>/interim/<STUDY_NAME>_<species>_annotated.rds
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SingleCellExperiment); library(zellkonverter)
  library(bluster); library(igraph)
})

args <- commandArgs(trailingOnly = TRUE)
STUDY_DIR <- if (length(args) >= 1) args[1] else Sys.getenv("STUDY_DIR", unset = NA)
if (is.na(STUDY_DIR)) stop("Usage: Rscript 02_integrate_and_cluster.R <study_dir>")
STUDY_DIR <- normalizePath(STUDY_DIR, mustWork = TRUE)
source(file.path(STUDY_DIR, "config.R"))
INTERIM <- file.path(STUDY_DIR, "interim")
emb_file <- function(sp) file.path(INTERIM, "uce", sprintf("%s_%s.h5ad", STUDY_NAME, sp))
ann_file <- function(sp) file.path(INTERIM, sprintf("%s_%s_annotated.rds", STUDY_NAME, sp))

get_uce <- function(sce) {
  nm <- reducedDimNames(sce); hit <- nm[grepl("uce", nm, ignore.case = TRUE)]
  if (!length(hit)) stop("No UCE embedding in reducedDims (", paste(nm, collapse = ", "), ").")
  reducedDim(sce, hit[1])
}

species <- c("human", "mouse")[file.exists(c(emb_file("human"), emb_file("mouse")))]
stopifnot(length(species) >= 1)
sces <- setNames(lapply(species, function(sp) { s <- readH5AD(emb_file(sp), use_hdf5 = FALSE); s$species <- sp; s }), species)

## ---- stack embeddings for JOINT (cross-species) clustering ----
embs <- lapply(sces, get_uce)
cellnames <- unlist(lapply(species, function(sp) paste0(sp, "::", colnames(sces[[sp]]))))
E <- do.call(rbind, embs); rownames(E) <- cellnames
message("Joint embedding: ", nrow(E), " cells x ", ncol(E), " dims")

set.seed(1234)  # Louvain and Annoy are randomised; fixed seed (C2)
g <- makeSNNGraph(E, k = K_NN, BNPARAM = BiocNeighbors::AnnoyParam())  # approximate k-NN (S7)
comm <- switch(CLUSTER_FUN,
               louvain  = cluster_louvain(g),
               leiden   = cluster_leiden(g, objective_function = "modularity"),
               walktrap = cluster_walktrap(g))
clusters <- factor(membership(comm))
message("Found ", nlevels(clusters), " joint clusters")

## ---- annotation ----
CELL_TYPE_COL <- get0("CELL_TYPE_COL", ifnotfound = NULL)
if (!is.null(CELL_TYPE_COL)) {
  ## Trusted labels (author annotation or fixture truth). Clusters are kept as a cross-check.
  cell_type <- unlist(lapply(sces, function(s) as.character(colData(s)[[CELL_TYPE_COL]])))
  message("Cell types from obs column '", CELL_TYPE_COL, "'. Joint clusters vs labels:")
  print(table(cluster = clusters, label = cell_type))
} else if (!is.null(REFERENCE)) {
  ref <- readH5AD(if (grepl("^(/|~)", REFERENCE)) REFERENCE else file.path(STUDY_DIR, REFERENCE),
                  use_hdf5 = FALSE)
  refE <- get_uce(ref); refLab <- as.character(colData(ref)[[REF_LABEL_COL]])
  kn <- BiocNeighbors::queryKNN(refE, E, k = K_TRANSFER)
  cell_type <- apply(kn$index, 1L, function(idx) names(sort(table(refLab[idx]), decreasing = TRUE))[1])
  message("Transferred labels; top types:"); print(head(sort(table(cell_type), decreasing = TRUE), 15))
} else {
  cell_type <- paste0("cluster_", as.character(clusters))
  message("No reference: joint clusters used as provisional cell types (relabel with markers).")
}

meta <- data.frame(study = unlist(lapply(sces, function(s) as.character(s$study))),
                   species = unlist(lapply(sces, function(s) as.character(s$species))),
                   cell_type = cell_type, row.names = cellnames)
message("\nStudies per cell type (want > 1):")
print(sort(tapply(meta$study, meta$cell_type, function(x) length(unique(x)))))
message("Species per cell type (2 = co-embedded across species):")
print(sort(tapply(meta$species, meta$cell_type, function(x) length(unique(x)))))
## Expectation: every cell type holds both conditions. A one-condition type cannot be tested in
## Stage 3; a strong disease state (e.g. lupus IFN-high monocytes) does this to unsupervised clusters.
meta$condition <- unlist(lapply(sces, function(s) as.character(s$condition)))
cond_tab <- table(meta$cell_type, meta$condition)
message("Cells per cell type x condition:"); print(cond_tab)
pure <- rownames(cond_tab)[rowSums(cond_tab > 0) < 2]
if (length(pure))
  message(sprintf("WARNING: %d of %d cell types contain one condition only (%s); clusters track condition, not cell type. Set CELL_TYPE_COL or REFERENCE.",
                  length(pure), nrow(cond_tab), paste(pure, collapse = ", ")))

split_lab <- split(meta$cell_type, meta$species)
for (sp in species) {
  lab <- split_lab[[sp]]
  names(lab) <- sub(paste0("^", sp, "::"), "", rownames(meta)[meta$species == sp])
  sce <- sces[[sp]]
  sce$cell_type <- lab[colnames(sce)]
  sce$joint_cluster <- clusters[match(paste0(sp, "::", colnames(sce)), cellnames)]
  saveRDS(sce, ann_file(sp))
  message("Wrote ", ann_file(sp), " (", ncol(sce), " cells)")
}
message("\nNext: Rscript pipeline/03_pseudobulk_de.R ", STUDY_DIR)
