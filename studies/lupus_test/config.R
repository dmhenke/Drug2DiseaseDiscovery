## config.R -- lupus TEST subsample (MACHINE TEST). Input: data/lupus_test.h5ad from
## scripts/subsample_h5ad.py. CHECK: obs column names and labels below are CELLxGENE
## conventions, NOT verified against the file. Run the dry-run and edit before use.
STUDY_NAME <- "lupus_test"
INPUTS <- list(
  list(path = "data/lupus_test.h5ad", study = "lupus_test", species = "human",
       obs_map = c(donor_id = "sample_id", disease = "condition"),  # CHECK
       symbols_from = "feature_name", counts_assay = "X")            # CHECK: var layout
)
CONDITION_MAP <- c(normal = "healthy", healthy = "healthy",
                   "systemic lupus erythematosus" = "diseased", diseased = "diseased")
MIN_GENES_PER_CELL <- 200
MIN_CELLS_PER_GENE <- 10
K_NN <- 30
CLUSTER_FUN <- "louvain"
REFERENCE <- NULL
REF_LABEL_COL <- "cell_type"
K_TRANSFER <- 30
MIN_CELLS_PER_SAMPLE <- 10
MIN_SAMPLES_PER_GROUP <- 3
MIN_COUNT <- 10
MIN_TOTAL_SAMPLES <- 2
ALPHA <- 0.05
COMPARE_ORTHOLOGS <- FALSE
