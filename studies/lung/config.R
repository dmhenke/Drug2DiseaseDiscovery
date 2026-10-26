## config.R  --  STUDY: lung
## ===========================================================================
## Per-study configuration. The pipeline scripts in ../../pipeline/ source this
## file; the CODE never changes per organ -- only this config and the data/ do.
## To add a new organ/illness, copy studies/_template/ to studies/<name>/ and
## edit its config.R. Nothing here overwrites another study.
## ===========================================================================

STUDY_NAME <- "lung"

## ---- INPUT DATASETS -------------------------------------------------------
## One entry per input .h5ad. `obs_map` renames colData columns to the canonical
## names (sample_id, condition). `symbols_from` = rowData column with gene
## symbols if rownames are Ensembl IDs (NULL if rownames are already symbols).
## Put the raw files in studies/lung/data/ and reference them relative to it.
INPUTS <- list(
  list(path = "data/studyA_human_lung.h5ad", study = "studyA", species = "human",
       obs_map = c(donor = "sample_id", disease_state = "condition"),
       symbols_from = NULL,           counts_assay = "counts"),
  list(path = "data/studyB_human_lung.h5ad", study = "studyB", species = "human",
       obs_map = c(patient = "sample_id", status = "condition"),
       symbols_from = "feature_name", counts_assay = "counts"),
  list(path = "data/studyC_mouse_lung.h5ad", study = "studyC", species = "mouse",
       obs_map = c(mouse_id = "sample_id", group = "condition"),
       symbols_from = "symbol",       counts_assay = "counts")
  ## ... add every lung dataset (human and mouse) here ...
)

## ---- CONDITION LABELS -----------------------------------------------------
## Map each dataset's raw condition labels to exactly {'healthy','diseased'}.
## Unlisted labels are dropped with a warning. Tailor the diseased terms to the
## illness(es) in THIS study (lung: IPF, COPD, fibrosis, tumor, ...).
CONDITION_MAP <- c(healthy = "healthy", normal = "healthy", control = "healthy",
                   ctrl = "healthy",
                   diseased = "diseased", disease = "diseased", ipf = "diseased",
                   copd = "diseased", fibrosis = "diseased", tumor = "diseased")

## ---- QC (Stage 0) ---------------------------------------------------------
MIN_GENES_PER_CELL <- 200    # paper CxG preprocessing
MIN_CELLS_PER_GENE <- 10

## ---- INTEGRATION / CLUSTERING (Stage 2) -----------------------------------
K_NN          <- 30          # SNN graph neighbours on the UCE embedding
CLUSTER_FUN   <- "louvain"   # louvain | leiden | walktrap
REFERENCE     <- NULL        # e.g. "../../shared/references/hlca_core_uce.h5ad" for label transfer
REF_LABEL_COL <- "cell_type"
K_TRANSFER    <- 30

## ---- DIFFERENTIAL EXPRESSION (Stage 3) ------------------------------------
MIN_CELLS_PER_SAMPLE  <- 10  # cells needed for a (cell_type,sample) pseudobulk
MIN_SAMPLES_PER_GROUP <- 2   # >=2 healthy AND >=2 diseased samples to test
MIN_COUNT             <- 10
MIN_TOTAL_SAMPLES     <- 2
ALPHA                 <- 0.05
COMPARE_ORTHOLOGS     <- TRUE  # human vs mouse ortholog intersection (needs babelgene)
