# =============================================================================
# 05_survival_summary.R — Collect HR / p-values from both cohorts into one file
#
# Inputs:  results/mitra_2024_n136/cox_metfs_hr_pvalue.csv   (03_mitra_2024.R)
#          results/tcga_mfs_ups_n75/cox_os_dss_hr_pvalue.csv (04_tcga_sarc.R)
# Output:  results/survival_HR_pvalues_Mitra_TCGA.xlsx
#
# Run from external_validation/:  Rscript scripts/05_survival_summary.R
# =============================================================================
source(file.path(if (file.exists("R/setup.R")) "." else "..", "R", "setup.R"))

mitra_csv <- file.path(paths$results, "mitra_2024_n136", "cox_metfs_hr_pvalue.csv")
tcga_csv  <- file.path(paths$results, "tcga_mfs_ups_n75", "cox_os_dss_hr_pvalue.csv")
require_files(mitra_csv, tcga_csv)

mitra <- read.csv(mitra_csv, check.names = FALSE)
tcga  <- read.csv(tcga_csv,  check.names = FALSE)

out <- file.path(paths$results, "survival_HR_pvalues_Mitra_TCGA.xlsx")
writexl::write_xlsx(list(All_Cohorts      = rbind(mitra, tcga),
                         Mitra_MetFS_n136 = mitra,
                         TCGA_OS_DSS_n75  = tcga), out)
print(rbind(mitra, tcga))
print(paste("Saved:", out))
