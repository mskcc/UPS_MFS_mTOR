# =============================================================================
# 02_prepare_data.R — Build processed inputs from the raw public data
#
# Inputs  (see README.md, Data):
#   data/raw/mitra_2024/path6347-sup-0002-tables1.xlsx          Mitra Table S1 (manual download)
#   data/raw/mitra_2024/E-MTAB-13947/P*_...isoforms.results     Mitra RSEM files (162)
#   data/raw/sarc_tcga_pan_can_atlas_2018/*.txt                 cBioPortal TCGA SARC files
# Outputs:
#   data/processed/mitra_2024_table_s1.csv     Table S1 as CSV
#   data/processed/mitra_2024_tpm_matrix.rds   gene-level TPM, genes x 162 samples
#   data/processed/tcga_sarc_data.rds          TCGA expression, RPPA, clinical
#
# Run from external_validation/:  Rscript scripts/02_prepare_data.R   (~3 min)
# =============================================================================
source(file.path(if (file.exists("R/setup.R")) "." else "..", "R", "setup.R"))
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(readxl)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})
dir.create(paths$processed, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Mitra Table S1 -> CSV
#    Title row skipped; footnote rows are dropped downstream (R/cohorts.R).
# =============================================================================
print(paste("Mitra 2024 Table S1..."))
require_files(files$mitra_table_s1)
s1 <- read_excel(files$mitra_table_s1, sheet = "Sample Annotations", skip = 1, col_types = "text")
fwrite(s1, files$mitra_clinical)
print(paste("Saved:", files$mitra_clinical))

# =============================================================================
# 2. Mitra RNA-seq (E-MTAB-13947) — aggregate RSEM isoforms to gene-level TPM
# =============================================================================
print(paste("Mitra 2024 gene-level TPM matrix (~3 min)..."))
rsem_files <- list.files(paths$raw_mitra_rna, pattern = "^P[0-9]+_fastq.*\\.results$", full.names = TRUE)
if (length(rsem_files) != 162) {
  stop("Expected 162 RSEM files in ", paths$raw_mitra_rna, ", found ", length(rsem_files),
       ". Run scripts/01_download_data.R.")
}
rsem_ids <- sub("_fastq_Aligned\\.toTranscriptome\\.out\\.isoforms\\.results$", "", basename(rsem_files))

gene_tpm_list <- lapply(rsem_files, function(f) {
  dt <- fread(f, select = c("gene_id", "TPM"))
  dt[, .(TPM = sum(TPM)), by = gene_id]
})
gene_ids_ref <- gene_tpm_list[[1]]$gene_id
tpm_mat <- do.call(cbind, lapply(gene_tpm_list, function(x) {
  stopifnot(identical(x$gene_id, gene_ids_ref))
  x$TPM
}))
rownames(tpm_mat) <- gene_ids_ref
colnames(tpm_mat) <- rsem_ids

# Ensembl ID (version stripped) -> HGNC symbol
ensg_nv  <- sub("\\.[0-9]+$", "", rownames(tpm_mat))
sym_map  <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(ensg_nv),
                                  keytype = "ENSEMBL", columns = "SYMBOL")
sym_map  <- sym_map[!is.na(sym_map$SYMBOL) & !duplicated(sym_map$ENSEMBL), ]
ensg2sym <- setNames(sym_map$SYMBOL, sym_map$ENSEMBL)

# GPR64 is ADGRG2 in current annotation — alias back to the signature name
if (!"GPR64" %in% ensg2sym && "ADGRG2" %in% ensg2sym) {
  ensg2sym[names(ensg2sym)[ensg2sym == "ADGRG2"]] <- "GPR64"
  print(paste("Aliased ADGRG2 -> GPR64"))
}

symbols <- ensg2sym[ensg_nv]
valid   <- !is.na(symbols)
tpm_mat <- tpm_mat[valid, ]
rownames(tpm_mat) <- symbols[valid]

# Duplicate symbols: keep the row with the highest mean expression
mean_e   <- rowMeans(tpm_mat)
keep_idx <- tapply(seq_len(nrow(tpm_mat)), rownames(tpm_mat), function(idx) idx[which.max(mean_e[idx])])
tpm_mat  <- tpm_mat[unlist(keep_idx), ]

print(paste(nrow(tpm_mat), "genes x", ncol(tpm_mat), "samples"))
saveRDS(tpm_mat, files$mitra_tpm)
print(paste("Saved:", files$mitra_tpm))

# =============================================================================
# 3. TCGA SARC PanCancer Atlas 2018 (cBioPortal)
# =============================================================================
print(paste("TCGA SARC PanCancer Atlas 2018..."))
tcga_file <- function(f) file.path(paths$raw_tcga, f)
require_files(tcga_file(c("data_mrna_seq_v2_rsem.txt", "data_rppa.txt",
                          "data_clinical_patient.txt", "data_clinical_sample.txt")))

tcga <- list()
tcga$RPPA            <- fread(tcga_file("data_rppa.txt"))
tcga$ClinicalPatient <- fread(tcga_file("data_clinical_patient.txt"), skip = 4)
tcga$ClinicalSample  <- fread(tcga_file("data_clinical_sample.txt"),  skip = 4)
# Deduplicate genes: keep first occurrence of each Hugo symbol
tcga$mrna_clean <- fread(tcga_file("data_mrna_seq_v2_rsem.txt")) %>%
  filter(!is.na(Hugo_Symbol) & Hugo_Symbol != "") %>%
  distinct(Hugo_Symbol, .keep_all = TRUE)

print(paste("mRNA:", nrow(tcga$mrna_clean), "genes,", ncol(tcga$mrna_clean) - 2, "samples"))
saveRDS(tcga, files$tcga)
print(paste("Saved:", files$tcga))

print(paste("02_prepare_data.R complete."))
