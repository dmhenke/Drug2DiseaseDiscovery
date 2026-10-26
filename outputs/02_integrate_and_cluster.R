#!/usr/bin/env Rscript
## 02_integrate_and_cluster.R
## ---------------------------------------------------------------------------
## Stage 2 of the UCE lung pipeline (R).
##
## Read the UCE-embedded species files, JOINTLY cluster them in the shared UCE
## space (mouse + human co-embed there), assign one consistent cell_type label
## set, and write per-species SCEs (with cell_type) for DE.
##
## All neighbour/cluster work uses the 1,280-dim UCE embedding, not expression.
##
##   Rscript 02_integrate_and_cluster.R
## ---------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(SingleCellExperiment)
  library(zellkonverter)
  library(bluster)
  library(igraph)
})

## ===========================================================================
## CONFIG
## ===========================================================================
EMB_FILES <- c(human = "uce_out/human_lung.h5ad",
               mouse = "uce_out/mouse_lung.h5ad")   # drop 'mouse' if none
OUT_FILES <- c(human = "human_lung_annotated.rds",
               mouse = "mouse_lung_annotated.rds")
K_NN      <- 30            # SNN graph neighbours
CLUSTER_FUN <- "louvain"   # igraph community detection on the SNN graph

## Optional reference atlas for label transfer: a UCE-embedded object (same model)
## with a cell_type column. Set to NULL to use de-novo clusters instead.
REFERENCE     <- NULL       # e.g. "hlca_core_uce.h5ad"
REF_LABEL_COL <- "cell_type"
K_TRANSFER    <- 30
## ===========================================================================

## find the reducedDim holding the UCE embedding (zellkonverter maps obsm here)
get_uce <- function(sce) {
  nm <- reducedDimNames(sce)
  hit <- nm[grepl("uce", nm, ignore.case = TRUE)]
  if (length(hit) == 0)
    stop("No UCE embedding in reducedDims (", paste(nm, collapse = ", "),
         "). Did you run UCE (run_uce.sh)?")
  reducedDim(sce, hit[1])
}

main <- function() {
  EMB_FILES <- EMB_FILES[file.exists(EMB_FILES)]
  stopifnot(length(EMB_FILES) >= 1)
  sces <- lapply(EMB_FILES, readH5AD, use_hdf5 = FALSE)
  for (sp in names(sces)) sces[[sp]]$species <- sp

  ## ---- stack embeddings across species for JOINT clustering ----
  embs <- lapply(sces, get_uce)
  cellnames <- unlist(lapply(names(sces), function(sp)
    paste0(sp, "::", colnames(sces[[sp]]))), use.names = FALSE)
  E <- do.call(rbind, embs)
  rownames(E) <- cellnames
  message("Joint embedding: ", nrow(E), " cells x ", ncol(E), " dims")

  ## ---- SNN graph + community detection on the UCE space ----
  g <- makeSNNGraph(E, k = K_NN)
  comm <- switch(CLUSTER_FUN,
                 louvain = cluster_louvain(g),
                 leiden  = cluster_leiden(g, objective_function = "modularity"),
                 walktrap = cluster_walktrap(g))
  clusters <- factor(membership(comm))
  message("Found ", nlevels(clusters), " joint clusters")

  ## ---- annotation ----
  if (!is.null(REFERENCE)) {
    ref <- readH5AD(REFERENCE, use_hdf5 = FALSE)
    refE <- get_uce(ref); refLab <- as.character(colData(ref)[[REF_LABEL_COL]])
    ## kNN majority-vote label transfer in the shared space
    library(BiocNeighbors)
    kn <- queryKNN(refE, E, k = K_TRANSFER)
    lab <- apply(kn$index, 1L, function(idx) {
      names(sort(table(refLab[idx]), decreasing = TRUE))[1]
    })
    cell_type <- lab
    message("Transferred labels; top types:")
    print(head(sort(table(cell_type), decreasing = TRUE), 15))
  } else {
    cell_type <- paste0("cluster_", as.character(clusters))
    message("No reference: using joint clusters as provisional cell types. ",
            "Relabel with markers before DE.")
  }

  ## ---- mixing diagnostic: do studies/species share cell types? ----
  meta <- data.frame(
    study   = unlist(lapply(sces, function(s) as.character(s$study))),
    species = unlist(lapply(sces, function(s) as.character(s$species))),
    cell_type = cell_type, row.names = cellnames)
  message("\nStudies contributing to each cell type (want > 1):")
  print(sort(tapply(meta$study, meta$cell_type, function(x) length(unique(x)))))
  message("\nSpecies per cell type (2 = co-embedded across species):")
  print(sort(tapply(meta$species, meta$cell_type, function(x) length(unique(x)))))

  ## ---- attach labels back to each per-species SCE and save ----
  split_lab <- split(meta$cell_type, meta$species)
  for (sp in names(sces)) {
    lab_sp <- split_lab[[sp]]
    names(lab_sp) <- sub(paste0("^", sp, "::"), "",
                         rownames(meta)[meta$species == sp])
    sce <- sces[[sp]]
    sce$cell_type <- lab_sp[colnames(sce)]
    sce$joint_cluster <- clusters[match(paste0(sp, "::", colnames(sce)), cellnames)]
    saveRDS(sce, OUT_FILES[[sp]])
    message("Wrote ", OUT_FILES[[sp]], " (", ncol(sce), " cells, counts + cell_type)")
  }
  message("\nNext: 03_pseudobulk_de.R")
}

main()
