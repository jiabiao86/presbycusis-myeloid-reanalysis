#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run with Rscript so repository paths can be resolved.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
source(file.path(dirname(script_file), "common_v14.R"))
configure_v14_library()
repo_root <- v14_repo_root()
prepared_dir <- v14_path(
  "V14_SINGLE_CELL_OUTPUT_DIR",
  file.path(
    repo_root, "..", "S11_myeloid_cell_audit", "gse274279_four_matrix"
  )
)
outdir <- v14_figure_dir(repo_root)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

umap <- read.csv(
  file.path(prepared_dir, "myeloid_four_matrix_umap_coordinates.csv"),
  row.names = 1,
  check.names = FALSE
)
proportions <- read.csv(
  file.path(
    prepared_dir, "myeloid_four_matrix_cluster_proportions_by_age.csv"
  ),
  check.names = FALSE
)
module_summary <- read.csv(
  file.path(
    prepared_dir, "myeloid_four_matrix_module_scores_by_state_age.csv"
  ),
  check.names = FALSE
)
bootstrap <- read.csv(
  file.path(prepared_dir, "myeloid_four_matrix_cell_bootstrap_samples.csv"),
  check.names = FALSE
)

panel_label <- function(label) {
  mtext(label, side = 3, adj = 0.01, line = 0.55, font = 2, xpd = NA)
}

draw <- function() {
  layout(
    matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE),
    heights = c(1, 1, 0.30)
  )
  par(family = "sans", mar = c(4.0, 4.5, 2.8, 0.8), mgp = c(2.5, 0.75, 0), las = 1)
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
  panel_label("A")

  ages <- sort(unique(proportions$age_months))
  states <- sort(unique(proportions$myeloid_state))
  matrix_prop <- matrix(0, nrow = length(states), ncol = length(ages))
  rownames(matrix_prop) <- states
  colnames(matrix_prop) <- ages
  for (i in seq_len(nrow(proportions))) {
    matrix_prop[proportions$myeloid_state[i], as.character(proportions$age_months[i])] <-
      proportions$proportion[i]
  }
  par(mar = c(4.0, 4.5, 4.2, 0.8), xpd = NA)
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
  legend_usr <- par("usr")
  legend(
    x = mean(legend_usr[1:2]),
    y = legend_usr[4] + diff(legend_usr[3:4]) * 0.04,
    xjust = 0.5,
    yjust = 0,
    legend = names(palette),
    fill = unname(palette),
    border = NA,
    horiz = TRUE,
    bty = "n",
    cex = 0.58,
    xpd = NA
  )
  panel_label("B")

  par(mar = c(4.0, 4.5, 2.8, 0.8), xpd = FALSE)
  score_columns <- grep("_score1$", colnames(module_summary), value = TRUE)
  score_matrix <- as.matrix(module_summary[, score_columns])
  rownames(score_matrix) <- paste0(module_summary$myeloid_state, " ", module_summary$age_months, "M")
  z <- t(scale(t(score_matrix)))
  z[!is.finite(z)] <- 0
  z <- pmax(pmin(z, 0.8), -0.8)
  heat_colors <- colorRampPalette(c("#2E6F9E", "white", "#C43C39"))(128)
  image(
    x = seq_len(ncol(z)),
    y = seq_len(nrow(z)),
    z = t(z),
    col = heat_colors,
    zlim = c(-0.8, 0.8),
    axes = FALSE,
    xlab = "",
    ylab = "",
    main = ""
  )
  axis(1, at = seq_len(ncol(z)), labels = sub("_score1$", "", colnames(z)), las = 2, cex.axis = 0.55)
  axis(2, at = seq_len(nrow(z)), labels = rownames(z), las = 2, cex.axis = 0.55)
  box()
  mtext("Row-scaled module score; display clipped to +/-0.8", side = 3, adj = 0.99,
        line = 0.45, cex = 0.46)
  panel_label("C")

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
  bootstrap_limits <- extendrange(
    range(c(summary_bootstrap$lower, summary_bootstrap$upper, 0), na.rm = TRUE),
    f = 0.06
  )
  plot(
    summary_bootstrap$median,
    y,
    xlim = bootstrap_limits,
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
  mtext(
    "500 nucleus-level bootstrap iterations; descriptive interval",
    side = 3, line = 0.55, cex = 0.48
  )
  panel_label("D")

  par(mar = c(3.0, 3.5, 1.3, 3.5), mgp = c(2.0, 0.55, 0), las = 1)
  breaks <- seq(-0.8, 0.8, length.out = length(heat_colors) + 1)
  plot.new()
  plot.window(xlim = c(-0.8, 0.8), ylim = c(0, 1))
  rect(head(breaks, -1), 0.20, tail(breaks, -1), 0.78, col = heat_colors, border = NA)
  box(lwd = 0.7)
  axis(
    3, at = seq(-0.8, 0.8, by = 0.2),
    labels = formatC(seq(-0.8, 0.8, by = 0.2), format = "f", digits = 1),
    cex.axis = 0.58, tick = FALSE, line = -0.55
  )
  mtext("Row-scaled module score", side = 1, line = 0.45, cex = 0.52)
  mtext(
    "Descriptive for sampled nuclei; not an animal-level confidence interval",
    side = 1,
    line = 2.05,
    cex = 0.44
  )
}

png(
  file.path(outdir, "figureS5_myeloid_subclustering.png"),
  width = 7.2,
  height = 7.0,
  units = "in",
  res = 300,
  bg = "white"
)
draw()
dev.off()

pdf(
  file.path(outdir, "figureS5_myeloid_subclustering.pdf"),
  width = 7.2,
  height = 7.0,
  bg = "white"
)
draw()
dev.off()

cat("Wrote Supplementary Figure S5 to", outdir, "\n")
