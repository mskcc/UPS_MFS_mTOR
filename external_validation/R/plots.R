# =============================================================================
# R/plots.R — Themes, plot builders, and save helpers
# =============================================================================

# -----------------------------------------------------------------------------
# Themes
# NOTE: ggsurvplot font.x / font.y / font.tickslab / font.legend must be
#       integers (e.g. 6) — decimals cause ggplot2 >= 3.5 validation errors.
# -----------------------------------------------------------------------------
base_theme <- theme_linedraw(base_size = 5.5) +
  theme(
    panel.grid.major  = element_line(color = "grey85", linewidth = 0.2),
    panel.grid.minor  = element_blank(),
    panel.background  = element_blank(),
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", colour = "black", size = 5.5),
    plot.title        = element_text(size = 5.5, face = "bold"),
    axis.text         = element_text(size = 5.5),
    axis.title        = element_text(size = 5.5),
    legend.text       = element_text(size = 5.5),
    legend.title      = element_text(size = 5.5, face = "bold")
  )

surv_theme <- base_theme + theme(plot.margin = margin(1, 1, 0, 1, "mm"))

risk_table_theme <- theme_linedraw(base_size = 5.5) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank(),
    plot.title       = element_text(size = 5.5, face = "bold"),
    axis.text.x      = element_blank(),
    axis.text.y      = element_text(size = 5.5),
    axis.title       = element_blank(),
    axis.ticks       = element_blank(),
    plot.margin      = margin(0, 1, 1, 1, "mm")
  )

# Text drawn in figures. The standard pdf() device renders "-" as a minus sign,
# so hyphens in labels are drawn with a soft hyphen (U+00AD), which both the PDF
# and PNG devices render as a normal hyphen. Data values keep the plain "-".
plot_label <- function(x) gsub("-", "\u00ad", x, fixed = TRUE)

# -----------------------------------------------------------------------------
# CONSORT diagram (Mitra cohort)
# -----------------------------------------------------------------------------
plot_consort <- function(n) {
  boxes <- data.frame(
    x = c(0, 0, 0), y = c(3, 2, 1),
    label = c(
      paste0("Mitra et al. 2024\nAll MFS/UPS patients\nn = ", n$n_total),
      paste0("Primary tumors only\nn = ", n$n_pt),
      paste0("Final analytic cohort\nn = ", n$n_final)
    )
  )
  exclusions <- data.frame(
    x = c(2.2, 2.2), y = c(2.5, 1.5),
    label = plot_label(c(
      paste0("Excluded (n = ", n$n_lr, ")\nLocal recurrences\n(LR1/LR2/LR3)"),
      paste0("Excluded (n = ", n$n_sync + n$n_na, ")\nSynchronous metastases (n = ", n$n_sync,
             ")\nNo MFS follow-up (n = ", n$n_na, ")")
    ))
  )
  arr <- arrow(length = unit(1.5, "mm"), type = "closed")
  ggplot() +
    geom_rect(data = boxes, aes(xmin = x - 1.1, xmax = x + 1.1, ymin = y - 0.3, ymax = y + 0.3),
              fill = "white", color = "black", linewidth = 0.4) +
    geom_text(data = boxes, aes(x = x, y = y, label = label), size = 2.5, lineheight = 1.1) +
    geom_rect(data = exclusions, aes(xmin = x - 1.2, xmax = x + 1.2, ymin = y - 0.25, ymax = y + 0.25),
              fill = "grey95", color = "grey40", linewidth = 0.3, linetype = "dashed") +
    geom_text(data = exclusions, aes(x = x, y = y, label = label), size = 2, color = "grey30", lineheight = 1.1) +
    geom_segment(aes(x = 0, xend = 0, y = 2.7, yend = 2.3), arrow = arr, linewidth = 0.3) +
    geom_segment(aes(x = 0, xend = 0, y = 1.7, yend = 1.3), arrow = arr, linewidth = 0.3) +
    geom_segment(aes(x = 0.4, xend = 1.0, y = 2.5, yend = 2.5), arrow = arr, linewidth = 0.3, color = "grey40") +
    geom_segment(aes(x = 0.4, xend = 1.0, y = 1.5, yend = 1.5), arrow = arr, linewidth = 0.3, color = "grey40") +
    geom_point(aes(x = 0, y = 2.5), shape = 16, size = 1) +
    geom_point(aes(x = 0, y = 1.5), shape = 16, size = 1) +
    scale_x_continuous(limits = c(-1.8, 4)) +
    scale_y_continuous(limits = c(0.5, 3.5)) +
    theme_void()
}

# -----------------------------------------------------------------------------
# Ranked GSVA lollipops
# -----------------------------------------------------------------------------
lollipop_df <- function(scores) {
  scores %>%
    arrange(desc(GSVA_score)) %>%
    mutate(Group = factor(Group, levels = c(group_high, group_low)))
}

# Horizontal ranked lollipop, colored by group, dashed line at the median
plot_lollipop_ranked <- function(scores, title) {
  lp_df <- lollipop_df(scores)
  med   <- attr(scores, "median")
  ggplot(lp_df, aes(x = reorder(Patient_ID, GSVA_score), y = GSVA_score)) +
    geom_segment(aes(xend = reorder(Patient_ID, GSVA_score), y = 0, yend = GSVA_score, color = Group),
                 linewidth = 0.3) +
    geom_point(aes(color = Group), size = 0.8, shape = 16, stroke = 0) +
    geom_hline(yintercept = med, linetype = "dashed", linewidth = 0.3, color = "grey30") +
    annotate("text", y = med, x = nrow(lp_df) * 0.95,
             label = paste0("Median = ", round(med, 3)),
             hjust = -0.05, vjust = 0.5, size = 1.5, color = "grey30") +
    scale_color_manual(values = group_colors, labels = plot_label) +
    coord_flip() +
    labs(x = NULL, y = "GSVA Enrichment Score", title = title) +
    base_theme + theme(axis.text.y = element_blank(),
                       axis.ticks.y = element_blank(),
                       legend.position = "bottom", legend.direction = "horizontal",
                       legend.key.size = unit(0.4, "lines"),
                       panel.grid.major = element_blank(), panel.grid.minor = element_blank())
}

# -----------------------------------------------------------------------------
# Kaplan-Meier (iCluster 1-low = reference level, plotted blue)
# -----------------------------------------------------------------------------
plot_km <- function(fit, data, ylab, title, font = 6, pval_size = 2,
                    table_fontsize = 2, table_height = 0.20, legend = "top") {
  km <- ggsurvplot(
    fit, data = data, pval = TRUE, pval.size = pval_size,
    risk.table = TRUE, risk.table.height = table_height, risk.table.fontsize = table_fontsize,
    risk.table.y.text = TRUE,
    palette = c("dodgerblue1", "firebrick1"), size = 0.5, censor.size = 1.5,
    surv.median.line = "hv",
    legend = legend,
    legend.title = "Group", legend.labs = plot_label(c(group_low, group_high)),
    xlab = "Time (months)", ylab = plot_label(ylab),
    title = plot_label(title),
    font.x = font, font.y = font, font.tickslab = font, font.legend = font,
    ggtheme = surv_theme, tables.theme = risk_table_theme
  )
  thin_median_lines(km)
}

# Thin median lines; drop survminer's colored strata rectangles in the risk table
thin_median_lines <- function(km_obj, lwd = 0.2) {
  km_obj$plot$layers <- lapply(km_obj$plot$layers, function(lyr) {
    if (inherits(lyr$geom, "GeomSegment")) lyr$aes_params$linewidth <- lwd
    lyr
  })
  km_obj$table$layers <- Filter(
    function(lyr) !inherits(lyr$geom, c("GeomTile", "GeomRect")),
    km_obj$table$layers
  )
  km_obj$table <- km_obj$table + theme(
    axis.ticks        = element_blank(),
    axis.ticks.length = unit(0, "pt")
  )
  km_obj
}

# KM curve + risk table stacked into one grob
km_grob <- function(km) {
  aligned <- cowplot::align_plots(km$plot, km$table, align = "v", axis = "lr")
  cowplot::plot_grid(aligned[[1]], aligned[[2]], ncol = 1, rel_heights = c(0.80, 0.20))
}

# -----------------------------------------------------------------------------
# Heatmaps
# -----------------------------------------------------------------------------
hm_style <- list(anno_name = 5, anno_size = 2, leg_title = 5, leg_label = 4, leg_height = 15,
                 leg_width = 2.5, disc_grid = 2, dend_width = 4)

# z: genes x samples z-scores (columns already ordered by GSVA score)
# anno: data.frame with GSVA_score and a categorical column `group_col`
signature_heatmap <- function(z, anno, group_col, group_colors, title, row_font, title_font) {
  s <- hm_style
  leg <- function() list(title_gp = gpar(fontsize = s$leg_title), labels_gp = gpar(fontsize = s$leg_label),
                         legend_height = unit(s$leg_height, "mm"), grid_width = unit(s$leg_width, "mm"))
  anno_cols <- list(
    `GSVA Score` = colorRamp2(c(min(anno$GSVA_score), median(anno$GSVA_score), max(anno$GSVA_score)),
                              c("#4B0082", "white", "#1a7a1a"))
  )
  anno_cols[[group_col]] <- group_colors
  anno_args <- list(`GSVA Score` = anno$GSVA_score)
  anno_args[[group_col]] <- anno[[group_col]]
  legend_params <- list(`GSVA Score` = leg())
  legend_params[[group_col]] <- list(title_gp = gpar(fontsize = s$leg_title),
                                     labels_gp = gpar(fontsize = s$leg_label),
                                     grid_height = unit(s$disc_grid, "mm"),
                                     grid_width = unit(s$disc_grid, "mm"))
  ha <- do.call(HeatmapAnnotation, c(anno_args, list(
    col = anno_cols,
    annotation_name_gp = gpar(fontsize = s$anno_name),
    simple_anno_size   = unit(s$anno_size, "mm"),
    annotation_legend_param = legend_params
  )))
  Heatmap(
    z, name = plot_label("Z-score"),
    col = colorRamp2(c(-2, 0, 2), c("dodgerblue1", "white", "firebrick1")),
    top_annotation = ha,
    cluster_columns = FALSE, cluster_rows = TRUE,
    row_split = signature_row_split(rownames(z)), row_title_gp = gpar(fontsize = row_font, fontface = "bold"),
    row_gap = unit(1, "mm"), show_column_names = FALSE,
    row_names_gp = gpar(fontsize = row_font, fontface = "italic"), border = TRUE,
    row_dend_width = unit(s$dend_width, "mm"), row_dend_gp = gpar(lwd = 0.5),
    column_title = title,
    column_title_gp = gpar(fontsize = title_font, fontface = "bold"),
    heatmap_legend_param = leg()
  )
}

# RPPA phospho-S6 rows aligned to the expression heatmap columns (NA if no RPPA)
rppa_heatmap <- function(rppa, sample_order, row_font) {
  s <- hm_style
  wide <- rppa %>%
    filter(Composite.Element.REF %in% s6_rows) %>%
    dplyr::select(Composite.Element.REF, all_of(intersect(sample_order, colnames(rppa)))) %>%
    tibble::column_to_rownames("Composite.Element.REF") %>%
    as.matrix()
  m <- matrix(NA_real_, nrow = nrow(wide), ncol = length(sample_order),
              dimnames = list(rownames(wide), sample_order))
  shared <- intersect(colnames(m), colnames(wide))
  m[, shared] <- wide[, shared]
  m <- t(scale(t(m)))
  # Space around "|" so it is not read as "l" in Helvetica (display only)
  rownames(m) <- gsub("|", " | ", rownames(m), fixed = TRUE)
  Heatmap(
    m, name = plot_label("RPPA\n(z-score)"),
    col = colorRamp2(c(-2, 0, 2), c("#7B2D8B", "white", "#DAA520")),
    cluster_columns = FALSE, cluster_rows = FALSE,
    show_column_names = FALSE, border = TRUE,
    row_names_gp = gpar(fontsize = row_font, fontface = "italic"),
    row_title = plot_label("Phospho-S6"), row_title_gp = gpar(fontsize = row_font, fontface = "bold"),
    na_col = "grey90",
    heatmap_legend_param = list(title_gp = gpar(fontsize = s$leg_title), labels_gp = gpar(fontsize = s$leg_label),
                                legend_height = unit(s$leg_height, "mm"), grid_width = unit(s$leg_width, "mm"))
  )
}

# -----------------------------------------------------------------------------
# Saving (PDF + 300 dpi PNG)
# -----------------------------------------------------------------------------
# draw_fn: a function with no arguments that draws the figure on the open device
save_figure <- function(draw_fn, path_no_ext, w, h) {
  dir.create(dirname(path_no_ext), recursive = TRUE, showWarnings = FALSE)
  # Standard pdf() device, as used by ggplot2::ggsave(): Helvetica, editable text
  pdf(paste0(path_no_ext, ".pdf"), width = w, height = h, useDingbats = FALSE)
  draw_fn()
  dev.off()
  png(paste0(path_no_ext, ".png"), width = w, height = h, units = "in", res = 300)
  draw_fn()
  dev.off()
  invisible(path_no_ext)
}

# ggplot / ggsurvplot / ComplexHeatmap object
save_plot <- function(p, path_no_ext, w, h) {
  draw_fn <- if (inherits(p, "ggsurvplot")) {
    function() print(km_grob(p))
  } else if (inherits(p, c("Heatmap", "HeatmapList"))) {
    function() draw(p, heatmap_legend_side = "right", merge_legend = TRUE)
  } else {
    function() print(p)
  }
  save_figure(draw_fn, path_no_ext, w, h)
}
