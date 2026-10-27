## config.R -- lupus FULL cohort (GPU deploy). Input: data/lupus_full.h5ad from
## scripts/subsample_h5ad.py --all-donors --max-cells-per-donor 0.
## CHECK before use: obs names and labels are CELLxGENE conventions, not verified.
STUDY_NAME <- "lupus_full"
INPUTS <- list(
  list(path = "data/lupus_full.h5ad", study = "perez2022", species = "human",
       obs_map = c(donor_id = "sample_id", disease = "condition"),  # CHECK
       symbols_from = "feature_name", counts_assay = "X")            # CHECK
)
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
