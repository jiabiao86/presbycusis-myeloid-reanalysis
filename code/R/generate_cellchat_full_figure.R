#!/usr/bin/env Rscript
# Figure 10: full CellChatDB.mouse communication across presbycusis stages.

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
outdir <- file.path(root, "manuscript", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

age_summary <- read.csv(
  file.path(root, "cellchat_full", "analysis", "CellChat_full_age_summary.csv"),
  check.names = FALSE
)
pathway_strength <- read.csv(
  file.path(root, "cellchat_full", "analysis", "CellChat_full_pathway_strength.csv"),
  check.names = FALSE
)
celltype_strength <- read.csv(
  file.path(root, "cellchat_full", "analysis", "CellChat_full_celltype_strength.csv"),
  check.names = FALSE
)

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
  layout(matrix(1:4, nrow = 2, byrow = TRUE))
  par(
    family = "sans",
    mar = c(4.2, 4.6, 2.2, 0.9),
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
  mtext("A", side = 3, adj = 0.01, line = -1.4, font = 2)

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
  mtext("B", side = 3, adj = 0.01, line = -1.4, font = 2)

  z <- pathway_matrix
  z[z == 0] <- NA
  image(
    x = seq_len(ncol(z)),
    y = seq_len(nrow(z)),
    z = t(z[nrow(z):1, , drop = FALSE]),
    col = colorRampPalette(c("#F7FBFF", "#9ECAE1", "#3182BD", "#08306B"))(64),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  axis(1, at = seq_len(ncol(z)), labels = colnames(z), cex.axis = 0.68)
  axis(2, at = seq_len(nrow(z)), labels = rev(rownames(z)), las = 2, cex.axis = 0.58)
  box()
  mtext("C", side = 3, adj = 0.01, line = -1.4, font = 2)

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
  label_set <- wide[order(-(wide$p3 + wide$p24)), , drop = FALSE]
  label_set <- head(label_set, 8)
  text(
    label_set$p3,
    label_set$p24,
    labels = label_set$pathway,
    pos = 4,
    cex = 0.5,
    col = "#1F3B57"
  )
  mtext("D", side = 3, adj = 0.01, line = -1.4, font = 2)
}

png(
  file.path(outdir, "figure10_cellchat_full_database.png"),
  width = 7.2,
  height = 7.0,
  units = "in",
  res = 300,
  bg = "white"
)
draw()
dev.off()

pdf(
  file.path(outdir, "figure10_cellchat_full_database.pdf"),
  width = 7.2,
  height = 7.0,
  bg = "white"
)
draw()
dev.off()

cat("Wrote Figure 10 to", outdir, "\n")
