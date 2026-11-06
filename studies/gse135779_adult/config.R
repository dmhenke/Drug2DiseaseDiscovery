## config.R -- Nehar-Belaid 2020 GSE135779, adult cohort (7 SLE / 5 healthy on GEO). Replication, analysed alone.
## Cell types are MARKER-BASED (scripts/annotate_markers.py), not author labels: GEO carries none.
STUDY_NAME <- "gse135779_adult"
INPUTS <- list(list(path = "data/gse135779_adult.h5ad", study = "gse135779_adult", species = "human",
  obs_map = c(donor = "sample_id"), symbols_from = NULL, counts_assay = "X"))
CONDITION_MAP <- c(normal = "healthy", healthy = "healthy",
                   "systemic lupus erythematosus" = "diseased", diseased = "diseased")
MIN_GENES_PER_CELL <- 200
MIN_CELLS_PER_GENE <- 10
K_NN <- 30
CLUSTER_FUN <- "louvain"
CELL_TYPE_COL <- "cell_type"
REFERENCE <- NULL
REF_LABEL_COL <- "cell_type"
K_TRANSFER <- 30
MIN_CELLS_PER_SAMPLE <- 10
MIN_SAMPLES_PER_GROUP <- 3
MIN_COUNT <- 10
MIN_TOTAL_SAMPLES <- 2
ALPHA <- 0.05
COMPARE_ORTHOLOGS <- FALSE
