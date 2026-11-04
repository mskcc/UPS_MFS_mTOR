# =============================================================================
# 01_download_data.R — Download the public input data (base R; works on Windows)
#
#   1. TCGA SARC PanCancer Atlas 2018, cBioPortal datahub   4 files, ~36 MB
#   2. Mitra et al. 2024 RNA-seq, ArrayExpress E-MTAB-13947  162 files, ~2.4 GB
#   3. Mitra et al. 2024 Table S1                            MANUAL — see message
#
# Files already present are skipped, so the script can be re-run after an
# interrupted download.
#
# Run from external_validation/:  Rscript scripts/01_download_data.R
# =============================================================================
source(file.path(if (file.exists("R/setup.R")) "." else "..", "R", "setup.R"))

options(timeout = max(3600, getOption("timeout")))

fetch <- function(url, dest) {
  if (file.exists(dest) && file.size(dest) > 0) return(invisible("exists"))
  dir.create(dirname(dest), recursive = TRUE, showWarnings = FALSE)
  tmp <- paste0(dest, ".part")
  status <- tryCatch(download.file(url, tmp, mode = "wb", quiet = TRUE), error = function(e) 1L)
  if (status != 0 || !file.exists(tmp) || file.size(tmp) == 0) {
    unlink(tmp)
    stop("Download failed: ", url, call. = FALSE)
  }
  file.rename(tmp, dest)
  print(paste("downloaded:", basename(dest)))
  invisible("downloaded")
}

# -----------------------------------------------------------------------------
# 1. TCGA SARC PanCancer Atlas 2018
#    https://www.cbioportal.org/study/summary?id=sarc_tcga_pan_can_atlas_2018
# -----------------------------------------------------------------------------
print(paste("TCGA SARC PanCancer Atlas 2018 (cBioPortal datahub)..."))
tcga_url <- "https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/sarc_tcga_pan_can_atlas_2018"
for (f in c("data_mrna_seq_v2_rsem.txt", "data_rppa.txt",
            "data_clinical_patient.txt", "data_clinical_sample.txt")) {
  fetch(paste0(tcga_url, "/", f), file.path(paths$raw_tcga, f))
}

# -----------------------------------------------------------------------------
# 2. Mitra et al. 2024 — E-MTAB-13947, 162 tumor RSEM isoform files
#    https://www.ebi.ac.uk/biostudies/arrayexpress/studies/E-MTAB-13947
#    (the 21 control tumors in the same accession are not used)
# -----------------------------------------------------------------------------
print(paste("Mitra 2024 RNA-seq, E-MTAB-13947 (ArrayExpress)..."))
ebi_url <- "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/947/E-MTAB-13947/Files"
fetch(paste0(ebi_url, "/E-MTAB-13947.sdrf.txt"), file.path(paths$raw_mitra_rna, "E-MTAB-13947.sdrf.txt"))

listing   <- readLines(paste0(ebi_url, "/"), warn = FALSE)
rsem_hits <- regmatches(listing, regexpr("P[0-9]+_fastq_Aligned\\.toTranscriptome\\.out\\.isoforms\\.results", listing))
rsem      <- sort(unique(rsem_hits))
if (length(rsem) != 162) stop("Expected 162 tumor RSEM files in the EBI listing, found ", length(rsem))

todo <- rsem[!file.exists(file.path(paths$raw_mitra_rna, rsem))]
print(paste(length(rsem) - length(todo), "of 162 already present;", length(todo), "to download"))
for (f in todo) fetch(paste0(ebi_url, "/", f), file.path(paths$raw_mitra_rna, f))

# -----------------------------------------------------------------------------
# 3. Mitra et al. 2024 — Table S1 (manual)
# -----------------------------------------------------------------------------
if (!file.exists(files$mitra_table_s1)) {
  message(
    "\nMANUAL STEP — Mitra et al. 2024 Table S1 (the publisher blocks scripted downloads):\n",
    "  1. Open https://doi.org/10.1002/path.6347\n",
    "  2. Go to 'Supporting Information' and download 'path6347-sup-0002-tables1.xlsx'\n",
    "  3. Save it, unchanged, as:\n     ", files$mitra_table_s1, "\n"
  )
} else {
  print(paste("Mitra Table S1 present."))
}

print(paste("01_download_data.R complete. Next: Rscript scripts/02_prepare_data.R"))
