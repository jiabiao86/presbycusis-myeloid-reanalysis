#!/usr/bin/env Rscript

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
outdir <- file.path(root, "manuscript", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

cellchat <- read.csv(
  file.path(root, "cellchat_formal", "CellChat_all_communications.csv"),
  check.names = FALSE
)
tf <- read.csv(
  file.path(root, "enhancement_analysis", "Severity_positive_FDR005_transcription_factors.csv"),
  check.names = FALSE
)
mirna <- read.csv(
  file.path(root, "enhancement_analysis", "Severity_positive_FDR005_micrornas.csv"),
  check.names = FALSE
)
go <- read.csv(
  file.path(root, "enhancement_analysis", "Severity_positive_FDR005_go_biological_process.csv"),
  check.names = FALSE
)
kegg <- read.csv(
  file.path(root, "enhancement_analysis", "Severity_positive_FDR005_kegg_mouse.csv"),
  check.names = FALSE
)
drug <- read.csv(
  file.path(root, "enhancement_analysis", "Severity_positive_FDR005_drug_repositioning.csv"),
  check.names = FALSE
)

draw_figure <- function() {
  layout(matrix(1:6, nrow = 2, byrow = TRUE))
  par(family = "sans", mar = c(4.0, 9.5, 2.0, 1.0), mgp = c(2.4, 0.7, 0), las = 1)

  cellchat$age_months <- factor(cellchat$age_months, levels = c(3, 12, 24))
  probabilities <- cellchat$prob[match(c(3, 12, 24), cellchat$age_months)]
  plot(
    c(3, 12, 24),
    probabilities,
    type = "b",
    pch = 21,
    bg = "#B24C63",
    col = "#B24C63",
    lwd = 2,
    ylim = c(0, max(probabilities) * 1.25),
    xlab = "Age (months)",
    ylab = "CellChat communication probability",
    xaxt = "n",
    cex.axis = 0.7,
    cex.lab = 0.78
  )
  axis(1, at = c(3, 12, 24), labels = c("3M", "12M", "24M"), cex.axis = 0.72)
  text(
    12,
    probabilities[2] * 1.12,
    paste0("PTPRC-MRC1\nmacrophage/microglia"),
    cex = 0.62
  )
  mtext("A", side = 3, adj = 0.01, line = -1.4, font = 2)

  tf_top <- head(tf[order(tf$adjusted_p_value), ], 10)
  tf_top <- tf_top[nrow(tf_top):1, ]
  barplot(
    -log10(tf_top$adjusted_p_value),
    horiz = TRUE,
    names.arg = substr(tf_top$term, 1, 38),
    col = "#4C78A8",
    border = NA,
    xlab = "-log10 FDR",
    cex.names = 0.50,
    cex.axis = 0.65
  )
  mtext("B", side = 3, adj = 0.01, line = -1.4, font = 2)

  mirna_top <- head(mirna[order(mirna$p_value), ], 10)
  mirna_top <- mirna_top[nrow(mirna_top):1, ]
  barplot(
    -log10(mirna_top$p_value),
    horiz = TRUE,
    names.arg = mirna_top$term,
    col = ifelse(mirna_top$adjusted_p_value < 0.05, "#D89A3C", "#BDBDBD"),
    border = NA,
    xlab = "-log10 P; none FDR<0.05",
    cex.names = 0.52,
    cex.axis = 0.65
  )
  mtext("C", side = 3, adj = 0.01, line = -1.4, font = 2)

  go_top <- head(go[order(go$adjusted_p_value), ], 10)
  go_top <- go_top[nrow(go_top):1, ]
  barplot(
    -log10(go_top$adjusted_p_value),
    horiz = TRUE,
    names.arg = substr(go_top$term, 1, 38),
    col = "#6BAED6",
    border = NA,
    xlab = "-log10 FDR",
    cex.names = 0.50,
    cex.axis = 0.65
  )
  mtext("D", side = 3, adj = 0.01, line = -1.4, font = 2)

  kegg_top <- head(kegg[order(kegg$adjusted_p_value), ], 10)
  kegg_top <- kegg_top[nrow(kegg_top):1, ]
  barplot(
    -log10(kegg_top$adjusted_p_value),
    horiz = TRUE,
    names.arg = substr(kegg_top$term, 1, 38),
    col = "#4C8C4A",
    border = NA,
    xlab = "-log10 FDR",
    cex.names = 0.50,
    cex.axis = 0.65
  )
  mtext("E", side = 3, adj = 0.01, line = -1.4, font = 2)

  drug_top <- head(drug[order(drug$adjusted_p_value), ], 10)
  drug_top <- drug_top[nrow(drug_top):1, ]
  barplot(
    -log10(drug_top$adjusted_p_value),
    horiz = TRUE,
    names.arg = substr(drug_top$term, 1, 38),
    col = "#8C6BB1",
    border = NA,
    xlab = "-log10 FDR",
    cex.names = 0.50,
    cex.axis = 0.65
  )
  mtext("F", side = 3, adj = 0.01, line = -1.4, font = 2)
}

png(
  file.path(outdir, "figure7_enhancement_analysis.png"),
  width = 7.2,
  height = 7.8,
  units = "in",
  res = 300,
  bg = "white"
)
draw_figure()
dev.off()

pdf(
  file.path(outdir, "figure7_enhancement_analysis.pdf"),
  width = 7.2,
  height = 7.8,
  bg = "white"
)
draw_figure()
dev.off()

cat("Wrote Figure 7 to", outdir, "\n")
