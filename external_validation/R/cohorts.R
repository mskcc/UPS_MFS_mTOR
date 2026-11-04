# =============================================================================
# R/cohorts.R — Cohort loading, filtering, and iCluster 1 GSVA scoring
# =============================================================================

# -----------------------------------------------------------------------------
# Mitra et al. 2024 (J Pathol, doi:10.1002/path.6347)
# -----------------------------------------------------------------------------

# Table S1 rows for the 162 MFS/UPS cases (drops footnote rows)
read_mitra_table_s1 <- function() {
  require_files(files$mitra_clinical)
  raw <- fread(files$mitra_clinical, header = TRUE, fill = TRUE)
  raw[grepl("^P[0-9]+$", raw[["Pat ID"]])]
}

# Base clinical variables (all MFS/UPS cases, primaries and local recurrences)
mitra_clinical_all <- function(raw) {
  raw %>%
    transmute(
      Patient_ID  = `Pat ID`,
      Diagnosis   = case_when(grepl("MFS", Diagnosisa) ~ "MFS",
                              grepl("UPS", Diagnosisa) ~ "UPS",
                              TRUE ~ NA_character_),
      MFS_Time    = suppressWarnings(as.numeric(MFS_Timee)),
      MFS_event   = as.integer(grepl("Met", Fud, ignore.case = FALSE)),
      Sample_type = Samplef
    ) %>%
    filter(!is.na(Diagnosis))
}

# Analytic cohort: primary tumors with MFS follow-up > 0 (n = 136)
mitra_analytic_cohort <- function(raw) {
  mitra_clinical_all(raw) %>%
    filter(grepl("^[Pp]$", Sample_type), MFS_Time > 0)
}

# Counts for the CONSORT diagram
mitra_consort_counts <- function(raw) {
  step1 <- mitra_clinical_all(raw)
  step2 <- step1 %>% filter(grepl("^[Pp]$", Sample_type))
  list(
    n_total = nrow(step1),
    n_lr    = nrow(step1) - nrow(step2),
    n_pt    = nrow(step2),
    n_sync  = sum(step2$MFS_Time <= 0, na.rm = TRUE),
    n_na    = sum(is.na(step2$MFS_Time)),
    n_final = nrow(step2 %>% filter(MFS_Time > 0))
  )
}

# log2(TPM + 1) gene x sample matrix
read_mitra_expression <- function() {
  require_files(files$mitra_tpm)
  log2(readRDS(files$mitra_tpm) + 1)
}

# -----------------------------------------------------------------------------
# TCGA SARC PanCancer Atlas 2018 (cBioPortal) — MFS (MFS) + UPS (MFH)
# -----------------------------------------------------------------------------
read_tcga <- function() {
  require_files(files$tcga)
  readRDS(files$tcga)
}

tcga_mfs_ups_samples <- function(tcga) {
  tcga$ClinicalSample %>% filter(ONCOTREE_CODE %in% c("MFH", "MFS"))
}

# One row per MFS/UPS sample with OS and DSS
tcga_clinical <- function(tcga) {
  mfs_ups <- tcga_mfs_ups_samples(tcga)
  tcga$ClinicalPatient %>%
    filter(PATIENT_ID %in% mfs_ups$PATIENT_ID) %>%
    mutate(
      Patient_ID = PATIENT_ID,
      OS_months  = as.numeric(OS_MONTHS),  OS_event  = as.integer(grepl("^1:", OS_STATUS)),
      DSS_months = as.numeric(DSS_MONTHS), DSS_event = as.integer(grepl("^1:", DSS_STATUS))
    ) %>%
    dplyr::select(Patient_ID, OS_months, OS_event, DSS_months, DSS_event) %>%
    inner_join(tcga$ClinicalSample %>%
                 dplyr::select(PATIENT_ID, SAMPLE_ID, ONCOTREE_CODE) %>%
                 rename(Patient_ID = PATIENT_ID, Sample_ID = SAMPLE_ID),
               by = "Patient_ID")
}

# log2(RSEM + 1) gene x sample matrix for MFS/UPS samples
tcga_expression <- function(tcga) {
  cols <- intersect(tcga_mfs_ups_samples(tcga)$SAMPLE_ID, colnames(tcga$mrna_clean))
  mat  <- as.matrix(tcga$mrna_clean[, ..cols])
  rownames(mat) <- tcga$mrna_clean$Hugo_Symbol
  log2(mat + 1)
}

# -----------------------------------------------------------------------------
# iCluster 1 signature scoring
# -----------------------------------------------------------------------------

# GSVA score per sample (named numeric vector)
score_icluster1 <- function(mat) {
  res <- gsva(gsvaParam(mat, gene_set, kcdf = "Gaussian"))
  if (inherits(res, "SummarizedExperiment")) res <- SummarizedExperiment::assay(res)
  setNames(as.numeric(res["iCluster1", ]), colnames(res))
}

# Median split of the GSVA score: >= median -> iCluster 1-high, else iCluster 1-low
assign_clusters <- function(df) {
  med <- median(df$GSVA_score, na.rm = TRUE)
  df$Group <- ifelse(df$GSVA_score >= med, group_high, group_low)
  attr(df, "median") <- med
  df
}

# Row-wise z-score of signature genes; drops genes absent from the matrix
signature_zscores <- function(mat, samples) {
  genes <- hm_genes[hm_genes %in% rownames(mat)]
  t(scale(t(mat[genes, samples, drop = FALSE])))
}

signature_row_split <- function(genes) {
  factor(ifelse(genes %in% icluster1_genes, "iCluster 1", "iCluster 2"),
         levels = c("iCluster 1", "iCluster 2"))
}
