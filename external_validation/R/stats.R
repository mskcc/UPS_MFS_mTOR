# =============================================================================
# R/stats.R — Survival statistics
# =============================================================================

# Cox model + log-rank test for iCluster 1-high vs iCluster 1-low (reference).
# `data` must have Group as a factor with iCluster 1-low as the first level.
cox_summary <- function(time, event, data, cohort, endpoint) {
  stopifnot(levels(data$Group)[1] == group_low)
  form <- as.formula(paste0("Surv(", time, ", ", event, ") ~ Group"))
  fit  <- coxph(form, data = data)
  s    <- summary(fit)
  lr   <- survdiff(form, data = data)
  data.frame(
    Cohort      = cohort,
    Endpoint    = endpoint,
    N           = s$n,
    Events      = s$nevent,
    Comparison  = paste(group_high, "vs", group_low),
    Reference   = group_low,
    HR          = round(s$conf.int[1, "exp(coef)"], 3),
    CI_lower_95 = round(s$conf.int[1, "lower .95"], 3),
    CI_upper_95 = round(s$conf.int[1, "upper .95"], 3),
    Logrank_p   = signif(1 - pchisq(lr$chisq, df = 1), 3),
    Cox_Wald_p  = signif(s$coefficients[1, "Pr(>|z|)"], 3),
    Cox_LRT_p   = signif(s$logtest[["pvalue"]], 3)
  )
}
