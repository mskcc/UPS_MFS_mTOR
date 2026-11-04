# =============================================================================
# run_all.R — Run the full pipeline in order
#
# From external_validation/:  Rscript run_all.R            (prepare + figures)
#                      Rscript run_all.R download   (also download data first)
# Each step runs in a fresh R process so steps cannot affect each other.
# =============================================================================
args  <- commandArgs(trailingOnly = TRUE)
steps <- c(
  if ("download" %in% args) "scripts/01_download_data.R",
  "scripts/02_prepare_data.R",
  "scripts/03_mitra_2024.R",
  "scripts/04_tcga_sarc.R",
  "scripts/05_survival_summary.R"
)
if (!file.exists("R/setup.R")) stop("Run from the external_validation/ folder.")

rscript <- file.path(R.home("bin"), "Rscript")
for (s in steps) {
  print(paste("==========", s, "=========="))
  status <- system2(rscript, shQuote(s))  # system2 quotes the command itself
  if (status != 0) stop(s, " failed (exit status ", status, ")")
}
print(paste("All steps complete."))
