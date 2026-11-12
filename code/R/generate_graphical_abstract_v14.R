#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run with Rscript so repository paths can be resolved.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
source(file.path(dirname(script_file), "common_v14.R"))
REPO_ROOT <- v14_repo_root()
analysis_root <- v14_prepared_dir(REPO_ROOT)
results_root <- v14_primary_results_dir(REPO_ROOT)
sensitivity_results_root <- v14_sensitivity_results_dir(REPO_ROOT)
outdir <- v14_figure_dir(REPO_ROOT)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

counts <- read.csv(
  file.path(results_root, "significant_counts.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
shared <- read.csv(
  file.path(results_root, "shared_genes.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
shared_sensitivity <- read.csv(
  file.path(sensitivity_results_root, "shared_genes.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

get_count <- function(tissue, contrast) {
  value <- counts$n_significant[
    counts$tissue == tissue & counts$contrast == contrast
  ]
  stopifnot(length(value) == 1)
  value
}

counts_matrix <- rbind(
  `MA vs YC` = c(
    Cochlea = get_count("cochlea", "MA_vs_YC"),
    `Inferior colliculus` = get_count("inferior_colliculus", "MA_vs_YC")
  ),
  `SP vs MP` = c(
    Cochlea = get_count("cochlea", "SP_vs_MP"),
    `Inferior colliculus` = get_count("inferior_colliculus", "SP_vs_MP")
  ),
  `SP vs YC` = c(
    Cochlea = get_count("cochlea", "SP_vs_YC"),
    `Inferior colliculus` = get_count("inferior_colliculus", "SP_vs_YC")
  )
)

shared_sp_yc <- shared$gene[shared$contrast == "SP_vs_YC"]
stopifnot(length(shared_sp_yc) == 7)
shared_sp_yc_sensitivity <- shared_sensitivity$gene[
  shared_sensitivity$contrast == "SP_vs_YC"
]
stopifnot(
  length(shared_sp_yc_sensitivity) == 8,
  "Tyrobp" %in% shared_sp_yc_sensitivity
)

draw_panel <- function(x0, x1, y0, y1, title) {
  rect(
    x0, y0, x1, y1,
    col = "#F7F8FA", border = "#D8DDE3", lwd = 1.2
  )
  rect(
    x0, y1 - 5.2, x1, y1,
    col = "#E9EEF3", border = NA
  )
  text(
    x0 + 1.3, y1 - 2.8, title,
    adj = c(0, 0.5), cex = 1.28, font = 2, col = "#203040"
  )
}

draw_graphical_abstract <- function() {
  par(mar = c(0, 0, 0, 0), family = "sans")
  plot.new()
  plot.window(xlim = c(0, 100), ylim = c(0, 52), xaxs = "i", yaxs = "i")

  text(
    2, 49.5,
    "Audited reanalysis of myeloid-associated genes in two auditory tissues",
    adj = c(0, 0.5), cex = 1.72, font = 2, col = "#183B56"
  )
  text(
    2, 45.6,
    "Unordered cohort labels; limited sex adjustment; tissues modeled separately",
    adj = c(0, 0.5), cex = 1.12, col = "#4A5563"
  )

  draw_panel(2, 30, 14.0, 42.5, "Cohort audit")
  text(
    3.5, 36.0, "Observed samples",
    adj = c(0, 0.5), cex = 1.05, font = 2, col = "#34495E"
  )

  bar_y <- c(32.2, 26.9)
  bar_values <- c(41, 39)
  bar_labels <- c("Cochlea", "Inferior colliculus")
  bar_colors <- c("#2E6F9E", "#A34D63")
  for (i in seq_along(bar_values)) {
    rect(
      3.5, bar_y[i] - 1.35, 3.5 + bar_values[i] * 0.48,
      bar_y[i] + 1.35, col = bar_colors[i], border = NA
    )
    text(
      4.2, bar_y[i],
      paste(bar_values[i], "samples"),
      adj = c(0, 0.5), cex = 0.90, font = 2, col = "white"
    )
    text(
      3.5, bar_y[i] + 1.85, bar_labels[i],
      adj = c(0, 0.5), cex = 0.92, col = "#4A5563"
    )
  }

  text(
    3.5, 22.8,
    "Animal level age, ABR, and DPOAE values: NA",
    adj = c(0, 0.5), cex = 0.86, col = "#7A3E3E"
  )
  text(
    3.5, 19.5,
    "Unmatched cochlea only: ID 13 (SP), ID 16 (YC)",
    adj = c(0, 0.5), cex = 0.76, col = "#4A5563"
  )
  text(
    3.5, 16.4,
    "Sample scenarios: 5 to 8 genes\nLeave one animal out: 5 to 10 genes\nNo-sex sensitivity: 7 intersection genes",
    adj = c(0, 0.5), cex = 0.70, col = "#4A5563"
  )

  draw_panel(33, 65.5, 14.0, 42.5, "Frozen contrast results")
  x_positions <- c(43.2, 56.2)
  colnames_short <- c("Cochlea", "Inferior colliculus")
  text(
    x_positions, 37.8, colnames_short,
    cex = 1.02, font = 2, col = "#34495E"
  )
  row_y <- c(32.7, 27.2, 21.7)
  row_labels <- c("MA vs YC", "SP vs MP", "SP vs YC")
  for (i in seq_along(row_y)) {
    text(
      34.4, row_y[i], row_labels[i],
      adj = c(0, 0.5), cex = 1.02, font = 2, col = "#4A5563"
    )
    for (j in seq_along(x_positions)) {
      value <- counts_matrix[i, j]
      fill <- if (value == 0) "#E8ECEF" else "#F2D8B8"
      text_col <- if (value == 0) "#59636D" else "#6B3B16"
      symbols(
        x_positions[j], row_y[i],
        circles = 2.25, inches = FALSE,
        bg = fill, fg = "#C7CDD3", lwd = 1.0, add = TRUE
      )
      text(
        x_positions[j], row_y[i], value,
        cex = 1.28, font = 2, col = text_col
      )
    }
  }
  text(
    34.4, 17.4,
    "Numbers are genes meeting BH FDR < 0.05 and |log2FC| >= 1",
    adj = c(0, 0.5), cex = 0.84, col = "#4A5563"
  )

  draw_panel(68.5, 98, 14.0, 42.5, "7 primary intersection genes")
  text(
    70.2, 35.7,
    "Intersection of separately\nsignificant genes (SP vs YC)",
    adj = c(0, 0.5), cex = 0.82, font = 2, col = "#34495E"
  )
  gene_rows <- 4
  gene_x <- c(70.2, 84.2)
  for (i in seq_along(shared_sp_yc)) {
    column <- if (i <= gene_rows) 1 else 2
    row <- if (i <= gene_rows) i else i - gene_rows
    x <- gene_x[column]
    y <- 32.4 - (row - 1) * 3.1
    points(x, y, pch = 24, bg = "#C46A2D", col = "#8D4C22", cex = 1.0)
    text(
      x + 1.8, y, shared_sp_yc[i],
      adj = c(0, 0.5), cex = 1.04, font = 2, col = "#243746"
    )
  }
  text(
    70.2, 20.6,
    "Leave one animal out intersection: 5 to 10 genes\n(median 7; Tyrobp in 6 of 41 models)",
    adj = c(0, 0.5), cex = 0.70, col = "#7A3E3E", font = 2
  )
  text(
    70.2, 17.4,
    "0 intersect genes for MA vs YC and SP vs MP\n0 biological WGCNA modules for those contrasts",
    adj = c(0, 0.5), cex = 0.60, col = "#4A5563"
  )

  rect(
    2, 0.8, 98, 11.4,
    col = "#EEF2F5", border = NA
  )
  text(
    3.5, 8.4,
    "Interpretation boundary",
    adj = c(0, 0.5), cex = 0.96, font = 2, col = "#7A3E3E"
  )
  text(
    3.5, 4.4,
    paste(
      "SP vs YC combines age and hearing status differences; MP and SP age ranges overlap.\n",
      "Single-nucleus detection is within marker-defined myeloid cells, not independent localization.\n",
      "Age and 10x batch are inseparable; CellChat is descriptive, not functional validation."
    ),
    adj = c(0, 0.5), cex = 0.88, col = "#4A5563"
  )
}

render_png <- function(path) {
  png(path, width = 12, height = 6.4, units = "in", res = 300, bg = "white")
  draw_graphical_abstract()
  dev.off()
}

render_pdf <- function(path) {
  pdf(path, width = 12, height = 6.4, bg = "white", useDingbats = FALSE)
  draw_graphical_abstract()
  dev.off()
}

render_tiff <- function(path) {
  tiff(
    path, width = 12, height = 6.4, units = "in", res = 300,
    compression = "lzw", bg = "white"
  )
  draw_graphical_abstract()
  dev.off()
}

render_png(file.path(outdir, "graphical_abstract_v14.png"))
render_pdf(file.path(outdir, "graphical_abstract_v14.pdf"))
render_tiff(file.path(outdir, "graphical_abstract_v14.tiff"))

cat("Wrote graphical_abstract_v14.png, .pdf, and .tiff\n")

