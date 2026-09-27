#!/usr/bin/env Rscript

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
outdir <- file.path(root, "manuscript", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cochlea_baseline <- read.csv(
  file.path(root, "limma_severity_trend.csv"),
  check.names = FALSE
)
central_baseline <- read.csv(
  file.path(root, "gse49522_severity_trend.csv"),
  check.names = FALSE
)
cochlea_sex <- read.csv(
  file.path(root, "sensitivity_analysis", "cochlea_severity_sex_adjusted.csv"),
  check.names = FALSE
)
cochlea_age <- read.csv(
  file.path(root, "sensitivity_analysis", "cochlea_severity_age_adjusted.csv"),
  check.names = FALSE
)
central_sex <- read.csv(
  file.path(root, "sensitivity_analysis", "central_severity_sex_adjusted.csv"),
  check.names = FALSE
)
central_age <- read.csv(
  file.path(root, "sensitivity_analysis", "central_severity_age_adjusted.csv"),
  check.names = FALSE
)
bootstrap <- read.csv(
  file.path(root, "sensitivity_analysis", "candidate_severity_bootstrap.csv"),
  check.names = FALSE
)

draw <- function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE))
  par(family = "sans", mar = c(4.2, 4.4, 2.0, 1.0), mgp = c(2.5, 0.75, 0), las = 1)

  counts <- matrix(
    c(
      sum(cochlea_baseline$adj.P.Val < 0.05),
      sum(cochlea_sex$adj.P.Val < 0.05),
      sum(cochlea_age$adj.P.Val < 0.05),
      sum(central_baseline$central_fdr_bh < 0.05),
      sum(central_sex$adj.P.Val < 0.05),
      sum(central_age$adj.P.Val < 0.05)
    ),
    nrow = 3,
    byrow = FALSE
  )
  colnames(counts) <- c("Cochlea", "Inferior colliculus")
  rownames(counts) <- c("Severity only", "Sex adjusted", "Age adjusted")
  barplot(
    counts,
    beside = TRUE,
    col = c("#4C78A8", "#E6B84A", "#C43C39"),
    border = NA,
    ylab = "Genes at FDR < 0.05",
    cex.names = 0.72,
    cex.axis = 0.7
  )
  legend("topright", legend = rownames(counts), fill = c("#4C78A8", "#E6B84A", "#C43C39"), border = NA, bty = "n", cex = 0.58)
  panel_label <- function(label) mtext(label, side = 3, adj = 0.01, line = -1.4, font = 2)
  panel_label("A")

  cochlea_base_sex <- cochlea_baseline$logFC[match(cochlea_sex$gene, cochlea_baseline$gene)]
  plot(
    cochlea_base_sex,
    cochlea_sex$logFC,
    pch = 21,
    bg = adjustcolor("#4C78A8", alpha.f = 0.35),
    col = NA,
    cex = 0.35,
    xlab = "Severity-only effect",
    ylab = "Sex-adjusted effect",
    cex.axis = 0.68,
    cex.lab = 0.75
  )
  abline(0, 1, lty = 3, col = "#777777")
  legend("topleft", legend = paste0("r = ", sprintf("%.3f", cor(cochlea_base_sex, cochlea_sex$logFC, use = "complete.obs"))), bty = "n", cex = 0.68)
  panel_label("B")

  central_base_sex <- central_baseline$central_severity_correlation[match(central_sex$gene, central_baseline$gene)]
  plot(
    central_base_sex,
    central_sex$logFC,
    pch = 21,
    bg = adjustcolor("#B24C63", alpha.f = 0.35),
    col = NA,
    cex = 0.35,
    xlab = "Severity-only correlation",
    ylab = "Sex-adjusted effect",
    cex.axis = 0.68,
    cex.lab = 0.75
  )
  abline(h = 0, v = 0, lty = 3, col = "#777777")
  legend("topleft", legend = paste0("r = ", sprintf("%.3f", cor(central_base_sex, central_sex$logFC, use = "complete.obs"))), bty = "n", cex = 0.68)
  panel_label("C")

  selected_genes <- c("H2-Aa", "Cd74", "Ctss", "Fcgr3", "Cd68", "Tyrobp", "Mpeg1", "C1qc", "C3ar1", "Clec7a")
  plot_data <- bootstrap[bootstrap$gene %in% selected_genes, ]
  plot_data <- plot_data[order(plot_data$tissue, match(plot_data$gene, selected_genes)), ]
  y <- seq_len(nrow(plot_data))
  colors <- ifelse(plot_data$tissue == "Cochlea", "#4C78A8", "#B24C63")
  plot(
    plot_data$median,
    y,
    xlim = range(c(plot_data$lower, plot_data$upper, 0), na.rm = TRUE),
    ylim = c(0.5, nrow(plot_data) + 0.5),
    pch = 21,
    bg = colors,
    col = "white",
    yaxt = "n",
    xlab = "Bootstrap severity coefficient",
    ylab = "",
    cex.axis = 0.62,
    cex.lab = 0.72
  )
  segments(plot_data$lower, y, plot_data$upper, y, col = colors, lwd = 1.4)
  abline(v = 0, lty = 3, col = "#777777")
  axis(2, at = y, labels = paste(plot_data$gene, plot_data$tissue, sep = " | "), las = 2, cex.axis = 0.52)
  panel_label("D")
}

png(
  file.path(outdir, "figure8_sensitivity_analysis.png"),
  width = 7.2,
  height = 7.2,
  units = "in",
  res = 300,
  bg = "white"
)
draw()
dev.off()

pdf(
  file.path(outdir, "figure8_sensitivity_analysis.pdf"),
  width = 7.2,
  height = 7.2,
  bg = "white"
)
draw()
dev.off()

cat("Wrote Figure 8 to", outdir, "\n")
