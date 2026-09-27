#!/usr/bin/env Rscript

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
outdir <- file.path(root, "manuscript", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

umap <- read.csv(
  file.path(root, "myeloid_subclustering", "myeloid_umap_coordinates.csv"),
  row.names = 1,
  check.names = FALSE
)
proportions <- read.csv(
  file.path(root, "myeloid_subclustering", "myeloid_cluster_proportions_by_age.csv"),
  check.names = FALSE
)
module_summary <- read.csv(
  file.path(root, "myeloid_subclustering", "myeloid_module_scores_by_state_age.csv"),
  check.names = FALSE
)
bootstrap <- read.csv(
  file.path(root, "myeloid_subclustering", "myeloid_module_bootstrap_samples.csv"),
  check.names = FALSE
)

draw <- function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE))
  par(family = "sans", mar = c(4.0, 4.5, 2.0, 0.8), mgp = c(2.5, 0.75, 0), las = 1)
  palette <- c("Myeloid_C0" = "#4C78A8", "Myeloid_C1" = "#C43C39")

  plot(
    umap$umap_1,
    umap$umap_2,
    type = "n",
    xlab = "UMAP1",
    ylab = "UMAP2",
    cex.axis = 0.68,
    cex.lab = 0.75
  )
  for (state in names(palette)) {
    subset <- umap$myeloid_state == state
    points(
      umap$umap_1[subset],
      umap$umap_2[subset],
      pch = 21,
      bg = adjustcolor(palette[state], alpha.f = 0.75),
      col = "white",
      cex = 0.9
    )
  }
  legend("topright", legend = names(palette), pt.bg = palette, pch = 21, col = "white", bty = "n", cex = 0.62)
  mtext("A", side = 3, adj = 0.01, line = -1.4, font = 2)

  ages <- sort(unique(proportions$age_months))
  states <- sort(unique(proportions$myeloid_state))
  matrix_prop <- matrix(0, nrow = length(states), ncol = length(ages))
  rownames(matrix_prop) <- states
  colnames(matrix_prop) <- ages
  for (i in seq_len(nrow(proportions))) {
    matrix_prop[proportions$myeloid_state[i], as.character(proportions$age_months[i])] <-
      proportions$proportion[i]
  }
  barplot(
    matrix_prop,
    col = palette[states],
    border = NA,
    ylim = c(0, 1),
    ylab = "Cell proportion",
    xlab = "Age (months)",
    cex.axis = 0.68,
    cex.lab = 0.75
  )
  mtext("B", side = 3, adj = 0.01, line = -1.4, font = 2)

  score_columns <- grep("_score1$", colnames(module_summary), value = TRUE)
  score_matrix <- as.matrix(module_summary[, score_columns])
  rownames(score_matrix) <- paste0(module_summary$myeloid_state, " ", module_summary$age_months, "M")
  z <- scale(score_matrix)
  z[!is.finite(z)] <- 0
  image(
    x = seq_len(ncol(z)),
    y = seq_len(nrow(z)),
    z = t(z),
    col = colorRampPalette(c("#2E6F9E", "white", "#C43C39"))(101),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = ""
  )
  axis(1, at = seq_len(ncol(z)), labels = sub("_score1$", "", colnames(z)), las = 2, cex.axis = 0.55)
  axis(2, at = seq_len(nrow(z)), labels = rownames(z), las = 2, cex.axis = 0.55)
  box()
  mtext("C", side = 3, adj = 0.01, line = -1.4, font = 2)

  bootstrap_list <- split(bootstrap, bootstrap$module)
  summary_bootstrap <- do.call(rbind, lapply(names(bootstrap_list), function(module) {
    values <- bootstrap_list[[module]]$delta_24m_minus_3m
    data.frame(
      module = module,
      lower = quantile(values, 0.025, na.rm = TRUE),
      median = median(values, na.rm = TRUE),
      upper = quantile(values, 0.975, na.rm = TRUE)
    )
  }))
  summary_bootstrap <- summary_bootstrap[order(summary_bootstrap$median), ]
  y <- seq_len(nrow(summary_bootstrap))
  colors <- ifelse(summary_bootstrap$median > 0, "#C43C39", "#4C78A8")
  plot(
    summary_bootstrap$median,
    y,
    xlim = range(c(summary_bootstrap$lower, summary_bootstrap$upper, 0), na.rm = TRUE),
    ylim = c(0.5, nrow(summary_bootstrap) + 0.5),
    pch = 21,
    bg = colors,
    col = "white",
    yaxt = "n",
    xlab = "Delta module score (24M - 3M)",
    ylab = "",
    cex.axis = 0.65,
    cex.lab = 0.72
  )
  segments(summary_bootstrap$lower, y, summary_bootstrap$upper, y, col = colors, lwd = 1.5)
  abline(v = 0, lty = 3, col = "#777777")
  axis(2, at = y, labels = summary_bootstrap$module, las = 2, cex.axis = 0.56)
  mtext("D", side = 3, adj = 0.01, line = -1.4, font = 2)
}

png(
  file.path(outdir, "figure9_myeloid_subclustering.png"),
  width = 7.2,
  height = 7.0,
  units = "in",
  res = 300,
  bg = "white"
)
draw()
dev.off()

pdf(
  file.path(outdir, "figure9_myeloid_subclustering.pdf"),
  width = 7.2,
  height = 7.0,
  bg = "white"
)
draw()
dev.off()

cat("Wrote Figure 9 to", outdir, "\n")
