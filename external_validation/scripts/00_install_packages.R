# =============================================================================
# 00_install_packages.R — Install the CRAN and Bioconductor packages used here
#
# Run once:  Rscript scripts/00_install_packages.R
# =============================================================================
cran <- c("data.table", "dplyr", "tidyr", "tibble", "ggplot2", "cowplot",
          "survival", "survminer", "circlize", "readxl", "writexl")
bioc <- c("GSVA", "ComplexHeatmap", "AnnotationDbi", "org.Hs.eg.db")

repos <- "https://cloud.r-project.org"
missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran, repos = repos)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = repos)
missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, update = FALSE, ask = FALSE)

still_missing <- c(cran, bioc)[!vapply(c(cran, bioc), requireNamespace, logical(1), quietly = TRUE)]
if (length(still_missing)) stop("Failed to install: ", paste(still_missing, collapse = ", "))
print(paste("All packages installed."))
