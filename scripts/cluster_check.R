#!/usr/bin/env Rscript
## cluster_check.R -- do the joint clusters track condition or donor rather than cell type?
## Reads <study>/interim/<STUDY_NAME>_<species>_annotated.rds (Stage 2 output).
## Expectation: every cluster holds both conditions near the overall SLE fraction, and no donor
## dominates a cluster. Within one label, a cluster skewed by condition is a composition shift.
## Usage: Rscript scripts/cluster_check.R studies/<name> [species]
suppressPackageStartupMessages(library(SingleCellExperiment))
args <- commandArgs(trailingOnly = TRUE)
study <- normalizePath(args[1], mustWork = TRUE); sp <- if (length(args) > 1) args[2] else "human"
source(file.path(study, "config.R"))
s <- readRDS(file.path(study, "interim", sprintf("%s_%s_annotated.rds", STUDY_NAME, sp)))
cl <- s$joint_cluster
t <- table(cl, s$condition)
d <- table(cl, s$sample_id)
message(sprintf("overall SLE fraction: %.2f", mean(s$condition == "diseased")))
print(cbind(t, frac_SLE = round(t[, "diseased"] / rowSums(t), 2),
            donors = rowSums(d > 0), top_donor_frac = round(apply(d, 1, max) / rowSums(d), 2)))
for (g in c("CCR7", "SELL", "GZMB", "IFI44L"))  # naive T, naive T, cytotoxic, interferon
  if (g %in% rownames(s)) {
    message("fraction of cells expressing ", g, " per cluster:")
    print(round(tapply(counts(s)[g, ] > 0, cl, mean), 2))
  }
for (lab in names(which(table(s$cell_type) >= 500))) {
  k <- s$cell_type == lab
  message("\n", lab, ": cluster x condition")
  print(table(cluster = droplevels(cl[k]), condition = s$condition[k]))
}
