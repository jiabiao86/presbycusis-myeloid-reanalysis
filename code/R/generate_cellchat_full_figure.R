#!/usr/bin/env Rscript
# Supplementary Figure S3: full CellChatDB.mouse communication across ages.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run with Rscript so repository paths can be resolved.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
source(file.path(dirname(script_file), "common_v14.R"))
configure_v14_library()
repo_root <- v14_repo_root()
analysis_dir <- file.path(
  v14_prepared_dir(repo_root), "cellchat_full", "analysis"
)
outdir <- v14_figure_dir(repo_root)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

age_summary <- read.csv(
  file.path(analysis_dir, "CellChat_full_age_summary.csv"),
  check.names = FALSE
)
pathway_strength <- read.csv(
  file.path(analysis_dir, "CellChat_full_pathway_strength.csv"),
  check.names = FALSE
)
celltype_strength <- read.csv(
  file.path(analysis_dir, "CellChat_full_celltype_strength.csv"),
  check.names = FALSE
)

panel_label <- function(label) {
  mtext(label, side = 3, adj = 0.01, line = 0.55, font = 2, xpd = NA)
}

ages <- sort(unique(age_summary$age_months))
age_labels <- paste0(ages, "M")
age_colors <- c("#4C78A8", "#F2A93B", "#C43C39")

top_pathways <- unique(unlist(lapply(ages, function(age) {
  frame <- pathway_strength[pathway_strength$age_months == age, , drop = FALSE]
  head(frame$pathway[order(-frame$total_probability)], 8)
})))
top_pathways <- top_pathways[top_pathways %in% pathway_strength$pathway]
pathway_matrix <- matrix(
  0,
  nrow = length(top_pathways),
  ncol = length(ages),
  dimnames = list(top_pathways, age_labels)
)
for (index in seq_len(nrow(pathway_strength))) {
  row <- pathway_strength[index, ]
  if (row$pathway %in% top_pathways) {
    pathway_matrix[row$pathway, paste0(row$age_months, "M")] <- row$total_probability
  }
}

draw <- function() {
  layout(
    matrix(c(1, 2, 3, 4, 5, 5), nrow = 3, byrow = TRUE),
    heights = c(1, 1, 0.18)
  )
  par(
    family = "sans",
    mar = c(4.2, 4.6, 2.8, 0.9),
    mgp = c(2.6, 0.75, 0),
    las = 1
  )

  edges <- age_summary$significant_edges
  interactions <- age_summary$significant_interactions
  pathways <- age_summary$significant_pathways
  midpoints <- barplot(
    edges,
    names.arg = age_labels,
    col = age_colors,
    border = NA,
    ylim = c(0, max(edges) * 1.25),
    ylab = "Significant edges",
    xlab = "Age (months)",
    cex.axis = 0.68,
    cex.lab = 0.75,
    main = ""
  )
  text(
    midpoints,
    edges,
    labels = paste0(interactions, " L-R\n", pathways, " pathways"),
    pos = 3,
    cex = 0.56,
    xpd = NA
  )
  panel_label("A")

  totals <- age_summary$total_significant_probability
  midpoints <- barplot(
    totals,
    names.arg = age_labels,
    col = age_colors,
    border = NA,
    ylim = c(0, max(totals) * 1.2),
    ylab = "Total communication probability",
    xlab = "Age (months)",
    cex.axis = 0.68,
    cex.lab = 0.75
  )
  text(
    midpoints,
    totals,
    labels = formatC(totals, format = "f", digits = 1),
    pos = 3,
    cex = 0.62,
    xpd = NA
  )
  panel_label("B")

  z <- pathway_matrix
  z[z == 0] <- NA
  z_max <- max(z, na.rm = TRUE)
  heat_colors <- colorRampPalette(c("#F7FBFF", "#9ECAE1", "#3182BD", "#08306B"))(128)
  image(
    x = seq_len(ncol(z)),
    y = seq_len(nrow(z)),
    z = t(z[nrow(z):1, , drop = FALSE]),
    col = heat_colors,
    zlim = c(0, z_max),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  axis(1, at = seq_len(ncol(z)), labels = colnames(z), cex.axis = 0.68)
  axis(2, at = seq_len(nrow(z)), labels = rev(rownames(z)), las = 2, cex.axis = 0.58)
  box()
  mtext("Pathway total probability", side = 3, adj = 0.99, line = 0.45, cex = 0.46)
  panel_label("C")

  wide <- reshape(
    pathway_strength[, c("pathway", "age_months", "total_probability")],
    idvar = "pathway",
    timevar = "age_months",
    direction = "wide"
  )
  names(wide) <- sub("total_probability\\.", "p", names(wide))
  for (column in c("p3", "p24")) {
    if (!column %in% names(wide)) wide[[column]] <- 0
  }
  wide$p3[is.na(wide$p3)] <- 0
  wide$p24[is.na(wide$p24)] <- 0
  plot(
    wide$p3,
    wide$p24,
    pch = 21,
    bg = adjustcolor("#4C78A8", alpha.f = 0.7),
    col = "white",
    xlab = "3M probability",
    ylab = "24M probability",
    cex.axis = 0.68,
    cex.lab = 0.75,
    xlim = c(0, max(c(wide$p3, wide$p24)) * 1.05),
    ylim = c(0, max(c(wide$p3, wide$p24)) * 1.05)
  )
  abline(0, 1, lty = 2, col = "grey45")
  label_genes <- c("PTN", "EPHA", "NRG", "PTPRM / NRXN")
  label_set <- wide[match(c("PTN", "EPHA", "NRG", "PTPRM"), wide$pathway), , drop = FALSE]
  label_x <- c(3.05, 2.18, 1.78, 0.72)
  label_y <- c(2.72, 1.58, 1.28, 0.95)
  segments(
    label_set$p3,
    label_set$p24,
    label_x,
    label_y,
    col = "#4A4A4A",
    lwd = 0.55
  )
  nrxn <- wide[match("NRXN", wide$pathway), , drop = FALSE]
  segments(
    nrxn$p3,
    nrxn$p24,
    label_x[4],
    label_y[4],
    col = "#4A4A4A",
    lwd = 0.55
  )
  text(
    label_x,
    label_y,
    labels = label_genes,
    adj = c(0, 0, 0, 0),
    cex = 0.56,
    col = "#1F3B57"
  )
  panel_label("D")

  par(mar = c(0.7, 6.0, 1.3, 6.0), mgp = c(2.0, 0.55, 0), las = 1)
  breaks <- seq(0, z_max, length.out = length(heat_colors) + 1)
  plot.new()
  plot.window(xlim = c(0, z_max), ylim = c(0, 1))
  rect(head(breaks, -1), 0.20, tail(breaks, -1), 0.78, col = heat_colors, border = NA)
  box(lwd = 0.7)
  ticks <- pretty(c(0, z_max), n = 5)
  ticks <- ticks[ticks >= 0 & ticks <= z_max]
  axis(
    3, at = ticks, labels = formatC(ticks, format = "f", digits = 2),
    cex.axis = 0.62, tick = FALSE, line = -0.55
  )
  mtext("Pathway total probability", side = 2, line = 1.0, cex = 0.58)
}

png(
  file.path(outdir, "figureS3_cellchat_age_overview.png"),
  width = 7.2,
  height = 7.0,
  units = "in",
  res = 300,
  bg = "white"
)
draw()
dev.off()

pdf(
  file.path(outdir, "figureS3_cellchat_age_overview.pdf"),
  width = 7.2,
  height = 7.0,
  bg = "white"
)
draw()
dev.off()

cat("Wrote Supplementary Figure S3 to", outdir, "\n")
