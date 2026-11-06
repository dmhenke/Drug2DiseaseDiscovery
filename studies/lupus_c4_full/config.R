## config.R -- Perez 2022 PBMC, all of Processing_Cohort 4 (26 SLE / 22 healthy female EA donors), author labels.
## Input built by scripts/run_expansion.sh with scripts/subsample_h5ad.py. Run in PCA mode: with CELL_TYPE_COL set
## the DE does not depend on the embedding.
STUDY_NAME <- "lupus_c4_full"
INPUTS <- list(list(path = "data/lupus_c4_full.h5ad", study = "perez2022_c4", species = "human",
  obs_map = c(donor_id = "sample_id", disease = "condition"), symbols_from = "feature_name", counts_assay = "X"))
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
