# =============================================================================
# 03_mitra_2024.R — Mitra et al., J Pathol 2024 cohort (n = 136 primary MFS/UPS)
#
# Figure panels (figures/*.pdf):
#   Fig1D_Mitra2024_MetFS_KM_n136               Metastasis-free survival Kaplan-Meier
#   SuppFig2B_Mitra2024_CONSORT                 Cohort selection (162 -> 145 -> 136)
#   SuppFig2C_Mitra2024_GSVA_ranked_n136        Ranked iCluster 1 GSVA scores, median split
#   SuppFig2D_Mitra2024_signature_heatmap_n136  iCluster 1 / iCluster 2 gene heatmap
# Statistics:
#   results/mitra_2024_n136/cox_metfs_hr_pvalue.csv
#
# Run from external_validation/:  Rscript scripts/03_mitra_2024.R
# =============================================================================
source(file.path(if (file.exists("R/setup.R")) "." else "..", "R", "setup.R"))
load_analysis_env()

out_dir <- file.path(paths$results, "mitra_2024_n136")
fig     <- function(name) file.path(paths$figures, name)

# -----------------------------------------------------------------------------
# Data: clinical (Table S1) + expression (E-MTAB-13947) + GSVA
# -----------------------------------------------------------------------------
raw      <- read_mitra_table_s1()
clinical <- mitra_analytic_cohort(raw)
counts   <- mitra_consort_counts(raw)
mat      <- read_mitra_expression()

gsva_all <- score_icluster1(mat)                         # scored on all 162 samples
scores   <- data.frame(Patient_ID = names(gsva_all), GSVA_score = unname(gsva_all)) %>%
  filter(Patient_ID %in% clinical$Patient_ID) %>%       # median taken within n = 136
  assign_clusters()

surv_df <- inner_join(clinical, scores, by = "Patient_ID") %>%
  mutate(Group = factor(Group, levels = c(group_low, group_high)))
n_final <- nrow(surv_df)
stopifnot(n_final == counts$n_final)
cohort  <- paste0("Mitra et al. 2024, n=", n_final)

# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------
cox_tab <- cox_summary("MFS_Time", "MFS_event", surv_df,
                       cohort = paste0("Mitra et al. 2024 (n=", n_final, ")"),
                       endpoint = "Metastasis-Free Survival")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(cox_tab, file.path(out_dir, "cox_metfs_hr_pvalue.csv"), row.names = FALSE)
print(cox_tab)

# -----------------------------------------------------------------------------
# Fig. 1D — Metastasis-free survival
# -----------------------------------------------------------------------------
metfs_fit <- survfit(Surv(MFS_Time, MFS_event) ~ Group, data = surv_df)
save_plot(plot_km(metfs_fit, surv_df, "Metastasis-Free Survival", paste0("MetFS (", cohort, ")")),
          fig(paste0("Fig1D_Mitra2024_MetFS_KM_n", n_final)), 2.6, 3)

# -----------------------------------------------------------------------------
# Supp. Fig. 2B — CONSORT diagram
# -----------------------------------------------------------------------------
save_plot(plot_consort(counts), fig("SuppFig2B_Mitra2024_CONSORT"), 4.5, 3.6)

# -----------------------------------------------------------------------------
# Supp. Fig. 2C — Ranked GSVA scores
# -----------------------------------------------------------------------------
save_plot(plot_lollipop_ranked(scores, paste0("iCluster 1 Signature (", cohort, ")")),
          fig(paste0("SuppFig2C_Mitra2024_GSVA_ranked_n", n_final)), 2.6, 3.4)

# -----------------------------------------------------------------------------
# Supp. Fig. 2D — Signature heatmap (z-scores over the 136 samples, ordered by GSVA)
# -----------------------------------------------------------------------------
hm_anno <- scores %>%
  left_join(clinical %>% dplyr::select(Patient_ID, Diagnosis), by = "Patient_ID") %>%
  arrange(GSVA_score)
hm_z    <- signature_zscores(mat, surv_df$Patient_ID)[, hm_anno$Patient_ID]
dx_cols <- c("MFS" = "#8DD69E", "UPS" = "#FAA44B")[na.omit(unique(hm_anno$Diagnosis))]

save_plot(signature_heatmap(hm_z, hm_anno, "Diagnosis", dx_cols,
                            paste0("iCluster Signatures (", cohort, ")"),
                            row_font = 5, title_font = 5),
          fig(paste0("SuppFig2D_Mitra2024_signature_heatmap_n", n_final)), 8.5, 2.5)

print(paste("03_mitra_2024.R complete."))
