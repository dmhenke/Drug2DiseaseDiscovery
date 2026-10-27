## config.R  --  STUDY TEMPLATE
## ===========================================================================
## Copy this whole folder to studies/<organ_or_illness>/ and edit below.
##   cp -r studies/_template studies/kidney
## Then put raw .h5ad files in that folder's data/ and run the pipeline with
## that study directory (see project README). Nothing you do here affects any
## other study -- separate data/, interim/, and results/ per study.
## ===========================================================================

STUDY_NAME <- "CHANGE_ME"          # e.g. "kidney", "liver_nafld"

## ---- INPUT DATASETS -------------------------------------------------------
INPUTS <- list(
  list(path = "data/example_human.h5ad", study = "studyX", species = "human",
       obs_map = c(donor = "sample_id", disease_state = "condition"),
       symbols_from = NULL, counts_assay = "counts")
  ## add human and/or mouse datasets ...
)

## ---- CONDITION LABELS  (edit diseased terms for THIS illness) -------------
CONDITION_MAP <- c(healthy = "healthy", normal = "healthy", control = "healthy",
                   diseased = "diseased", disease = "diseased")

## ---- QC (Stage 0) ---------------------------------------------------------
MIN_GENES_PER_CELL <- 200
MIN_CELLS_PER_GENE <- 10

## ---- INTEGRATION / CLUSTERING (Stage 2) -----------------------------------
K_NN          <- 30
CLUSTER_FUN   <- "louvain"
# CELL_TYPE_COL <- "cell_type"  # obs column of trusted labels (CELLxGENE author annotation); NULL = clusters
REFERENCE     <- NULL              # path to a UCE-embedded reference atlas, or NULL
REF_LABEL_COL <- "cell_type"
K_TRANSFER    <- 30

## ---- DIFFERENTIAL EXPRESSION (Stage 3) ------------------------------------
MIN_CELLS_PER_SAMPLE  <- 10
MIN_SAMPLES_PER_GROUP <- 2
MIN_COUNT             <- 10
MIN_TOTAL_SAMPLES     <- 2
ALPHA                 <- 0.05
COMPARE_ORTHOLOGS     <- TRUE
