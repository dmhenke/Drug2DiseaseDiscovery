## config.R -- lupus TEST subsample (MACHINE TEST), Perez 2022 PBMC (CELLxGENE c55dc602).
## Processing_Cohort 1 is healthy-only, so it is excluded; cohorts 2-4 hold both conditions and each
## is its own "study", so ~ study + condition adjusts for processing batch. Female, European American
## only: male and Asian cells are SLE-heavy (obs crosstab, 2026-10-27). Built by scripts/subsample_h5ad.py.
STUDY_NAME <- "lupus_test"
cohort <- function(k) list(path = sprintf("data/lupus_c%d.h5ad", k), study = sprintf("cohort%d", k),
  species = "human", obs_map = c(donor_id = "sample_id", disease = "condition"),
  symbols_from = "feature_name", counts_assay = "X")
INPUTS <- lapply(2:4, cohort)
CONDITION_MAP <- c(normal = "healthy", healthy = "healthy",
                   "systemic lupus erythematosus" = "diseased", diseased = "diseased")
MIN_GENES_PER_CELL <- 200
MIN_CELLS_PER_GENE <- 10
K_NN <- 30
CLUSTER_FUN <- "louvain"
# CELL_TYPE_COL <- "cell_type"  # obs column of trusted labels (CELLxGENE author annotation); NULL = clusters
REFERENCE <- NULL
REF_LABEL_COL <- "cell_type"
K_TRANSFER <- 30
MIN_CELLS_PER_SAMPLE <- 10
MIN_SAMPLES_PER_GROUP <- 3
MIN_COUNT <- 10
MIN_TOTAL_SAMPLES <- 2
ALPHA <- 0.05
COMPARE_ORTHOLOGS <- FALSE
