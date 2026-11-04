# =============================================================================
# 04_tcga_sarc.R — TCGA SARC PanCancer Atlas, MFS + UPS (n = 75)
#
# Figure panels (figures/*.pdf):
#   Fig1E_TCGA_MFS_UPS_OS_KM_n75                      Overall survival Kaplan-Meier
#   SuppFig2E_TCGA_MFS_UPS_GSVA_ranked_n75            Ranked iCluster 1 GSVA scores, median split
#   SuppFig2F_TCGA_MFS_UPS_signature_heatmap_n75      Gene heatmap + RPPA phospho-S6
#   SuppFig2G_TCGA_MFS_UPS_DSS_KM_n75                 Disease-specific survival Kaplan-Meier
# Statistics:
#   results/tcga_mfs_ups_n75/cox_os_dss_hr_pvalue.csv
#
# Run from external_validation/:  Rscript scripts/04_tcga_sarc.R
# =============================================================================
source(file.path(if (file.exists("R/setup.R")) "." else "..", "R", "setup.R"))
load_analysis_env()

out_dir <- file.path(paths$results, "tcga_mfs_ups_n75")
fig     <- function(name) file.path(paths$figures, name)

# -----------------------------------------------------------------------------
# Data: clinical + log2(RSEM + 1) + GSVA
# -----------------------------------------------------------------------------
tcga     <- read_tcga()
clinical <- tcga_clinical(tcga)
mat      <- tcga_expression(tcga)

gsva_all <- score_icluster1(mat)
scores   <- data.frame(Sample_ID = names(gsva_all), GSVA_score = unname(gsva_all)) %>%
  assign_clusters()
scores$Patient_ID <- gsub("-01$", "", scores$Sample_ID)

surv_df <- inner_join(clinical, scores, by = c("Patient_ID", "Sample_ID")) %>%
  mutate(Group = factor(Group, levels = c(group_low, group_high)))
surv_dss <- surv_df %>% filter(!is.na(DSS_months) & !is.na(DSS_event))
n_final  <- nrow(surv_df)
cohort   <- paste0("TCGA MFS/UPS, n=", n_final)

# -----------------------------------------------------------------------------
# Statistics
# -----------------------------------------------------------------------------
cox_tab <- rbind(
  cox_summary("OS_months",  "OS_event",  surv_df,  paste0("TCGA (n=", n_final, ")"), "Overall Survival"),
  cox_summary("DSS_months", "DSS_event", surv_dss, paste0("TCGA (n=", n_final, ")"), "Disease-Specific Survival")
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write.csv(cox_tab, file.path(out_dir, "cox_os_dss_hr_pvalue.csv"), row.names = FALSE)
print(cox_tab)

# -----------------------------------------------------------------------------
# Fig. 1E — Overall survival
# -----------------------------------------------------------------------------
os_fit <- survfit(Surv(OS_months, OS_event) ~ Group, data = surv_df)
save_plot(plot_km(os_fit, surv_df, "Overall Survival", paste0("OS (", cohort, ")")),
          fig(paste0("Fig1E_TCGA_MFS_UPS_OS_KM_n", n_final)), 2.6, 3)

# -----------------------------------------------------------------------------
# Supp. Fig. 2E — Ranked GSVA scores
# -----------------------------------------------------------------------------
save_plot(plot_lollipop_ranked(scores, paste0("iCluster 1 Signature (", cohort, ")")),
          fig(paste0("SuppFig2E_TCGA_MFS_UPS_GSVA_ranked_n", n_final)), 2.6, 3.4)

# -----------------------------------------------------------------------------
# Supp. Fig. 2F — Signature heatmap + RPPA phospho-S6 (z-scores over all 75 samples)
# -----------------------------------------------------------------------------
hm_anno <- scores %>%
  left_join(surv_df %>% dplyr::select(Sample_ID, ONCOTREE = ONCOTREE_CODE), by = "Sample_ID") %>%
  arrange(GSVA_score)
hm_z <- signature_zscores(mat, colnames(mat))[, hm_anno$Sample_ID]

hm <- signature_heatmap(hm_z, hm_anno, "ONCOTREE", c("MFH" = "#FAA44B", "MFS" = "#8DD69E"),
                        paste0("iCluster Signatures (", cohort, ")"), row_font = 6, title_font = 7) %v%
  rppa_heatmap(tcga$RPPA, hm_anno$Sample_ID, row_font = 6)
save_plot(hm, fig(paste0("SuppFig2F_TCGA_MFS_UPS_signature_heatmap_n", n_final)), 7.5, 2.7)

# -----------------------------------------------------------------------------
# Supp. Fig. 2G — Disease-specific survival
# -----------------------------------------------------------------------------
dss_fit <- survfit(Surv(DSS_months, DSS_event) ~ Group, data = surv_dss)
save_plot(plot_km(dss_fit, surv_dss, "Disease-Specific Survival", paste0("DSS (", cohort, ")")),
          fig(paste0("SuppFig2G_TCGA_MFS_UPS_DSS_KM_n", n_final)), 2.6, 3)

print(paste("04_tcga_sarc.R complete."))
