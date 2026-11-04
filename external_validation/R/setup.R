# =============================================================================
# R/setup.R — Project paths, packages, shared constants
#
# Sourced at the top of every script in scripts/. All paths are relative to the
# external_validation/ folder, so the repository can be cloned anywhere (macOS, Linux, Windows).
# Run scripts from external_validation/ (e.g. `Rscript scripts/03_mitra_2024.R`);
# running from inside scripts/ also works.
# =============================================================================

# -----------------------------------------------------------------------------
# Paths
# -----------------------------------------------------------------------------
repo_root <- normalizePath(if (file.exists("R/setup.R")) "." else "..", winslash = "/")
if (!file.exists(file.path(repo_root, "R", "setup.R"))) {
  stop("Run from the external_validation/ folder (the one containing R/ and scripts/).")
}

paths <- list(
  raw_tcga      = file.path(repo_root, "data", "raw", "sarc_tcga_pan_can_atlas_2018"),
  raw_mitra     = file.path(repo_root, "data", "raw", "mitra_2024"),
  raw_mitra_rna = file.path(repo_root, "data", "raw", "mitra_2024", "E-MTAB-13947"),
  processed     = file.path(repo_root, "data", "processed"),
  figures       = file.path(repo_root, "figures"),
  results       = file.path(repo_root, "results")
)

files <- list(
  mitra_table_s1   = file.path(paths$raw_mitra, "path6347-sup-0002-tables1.xlsx"),
  mitra_clinical   = file.path(paths$processed, "mitra_2024_table_s1.csv"),
  mitra_tpm        = file.path(paths$processed, "mitra_2024_tpm_matrix.rds"),
  tcga             = file.path(paths$processed, "tcga_sarc_data.rds")
)

# -----------------------------------------------------------------------------
# Shared constants
# -----------------------------------------------------------------------------
icluster1_genes <- c("TNC", "TMEM158", "STEAP1", "SERPINE1", "ITGA10", "LRRC15")
icluster2_genes <- c("CLEC3B", "TNXB", "MFAP5", "GPR64", "ABCA8")
hm_genes        <- c(icluster1_genes, icluster2_genes)
gene_set        <- list(iCluster1 = icluster1_genes)
s6_rows         <- c("RPS6|S6_pS235_S236", "RPS6|S6_pS240_S244", "RPS6|S6")

# Sample groups: median split of the iCluster 1 GSVA score within each cohort
group_high   <- "iCluster 1-high"
group_low    <- "iCluster 1-low"
group_colors <- setNames(c("firebrick1", "dodgerblue1"), c(group_high, group_low))

# -----------------------------------------------------------------------------
# Packages and project functions (only when called from an analysis script)
# -----------------------------------------------------------------------------
load_analysis_env <- function() {
  suppressPackageStartupMessages({
    library(data.table)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(cowplot)
    library(survival)
    library(survminer)
    library(GSVA)
    library(ComplexHeatmap)
    library(circlize)
  })
  for (f in c("cohorts.R", "plots.R", "stats.R")) {
    sys.source(file.path(repo_root, "R", f), envir = globalenv())
  }
  # Under Rscript, give stray drawing calls a null device instead of Rplots.pdf
  if (!interactive() && dev.cur() == 1) pdf(NULL)
  invisible(NULL)
}

# Stop with a clear message if an input file is missing
require_files <- function(...) {
  f <- c(...)
  missing <- f[!file.exists(f)]
  if (length(missing)) {
    stop("Missing input file(s):\n  ", paste(missing, collapse = "\n  "),
         "\nSee README.md (Data) and run scripts/01_download_data.R and scripts/02_prepare_data.R.",
         call. = FALSE)
  }
}
