#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
set.seed(20261028)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run with Rscript so repository paths can be resolved.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
source(file.path(dirname(script_file), "common_v14.R"))
REPO_ROOT <- v14_repo_root()
analysis_root <- v14_prepared_dir(REPO_ROOT)
four_matrix_dir <- file.path(
  REPO_ROOT,
  "..",
  "S11_myeloid_cell_audit",
  "gse274279_four_matrix"
)
external_results_dir <- file.path(REPO_ROOT, "results", "external")
results_root <- v14_primary_results_dir(REPO_ROOT)
sensitivity_results_root <- v14_path(
  "V14_SENSITIVITY_RESULTS_DIR",
  file.path(
    REPO_ROOT, "..", "S9_significant_and_intersection",
    "sensitivity_exclude_GSM1201568"
  )
)
exclude_1602_results_root <- v14_path(
  "V14_EXCLUDE_1602_RESULTS_DIR",
  file.path(
    REPO_ROOT, "..", "S9_significant_and_intersection",
    "sensitivity_exclude_GSM1201602"
  )
)
exclude_both_results_root <- v14_path(
  "V14_EXCLUDE_BOTH_RESULTS_DIR",
  file.path(
    REPO_ROOT, "..", "S9_significant_and_intersection",
    "sensitivity_exclude_both_unmatched"
  )
)
loo_results_root <- v14_path(
  "V14_LEAVE_ONE_ANIMAL_OUT_RESULTS_DIR",
  file.path(
    REPO_ROOT, "..", "S9_significant_and_intersection", "leave_one_animal_out"
  )
)
no_sex_results_root <- file.path(
  REPO_ROOT, "..", "S9_significant_and_intersection", "sensitivity_no_sex"
)
source_figure_dir <- v14_figure_dir(REPO_ROOT)
outdir <- v14_figure_dir(REPO_ROOT)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
archived_figure_dir <- file.path(REPO_ROOT, "figures", "figures_v14")

groups <- c("YC", "MA", "MP", "SP")
group_labels <- c(
  YC = "Young control",
  MA = "Middle-aged",
  MP = "MP (source label)",
  SP = "SP (source label)"
)
group_short <- c(YC = "YC", MA = "MA", MP = "MP", SP = "SP")
group_colors <- c(
  YC = "#4C78A8",
  MA = "#D9A441",
  MP = "#D66B3D",
  SP = "#A93F55"
)
tissue_colors <- c(
  Cochlea = "#2E6F9E",
  `Inferior colliculus` = "#A34D63"
)

read_csv <- function(path) {
  read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
}

panel <- function(label, cex = 1.18, line = 0.4) {
  mtext(label, side = 3, adj = 0.01, line = line, font = 2, cex = cex, xpd = NA)
}

wrap <- function(x, width = 28) {
  vapply(x, function(v) paste(strwrap(v, width = width), collapse = "\n"),
         character(1))
}

smart_labels <- function(x, y, labels, cex = 0.55, col = "#333333",
                         seg_col = "#8A8A8A") {
  usr <- par("usr")
  x_range <- diff(usr[1:2])
  y_range <- diff(usr[3:4])
  widths <- strwidth(labels, cex = cex, units = "user")
  heights <- strheight(labels, cex = cex, units = "user")
  x_offsets <- c(0, 0.035, 0.070, -0.035, -0.070,
                 0.035, 0.070, -0.035, -0.070, 0) * x_range
  y_offsets <- c(0.055, 0.055, 0.055, 0.055, 0.055,
                 -0.055, -0.055, -0.055, -0.055, 0.100) * y_range
  placed <- matrix(numeric(0), ncol = 4)

  for (index in order(y, decreasing = TRUE)) {
    candidates <- data.frame(
      cx = x[index] + x_offsets,
      cy = y[index] + y_offsets
    )
    candidates$left <- candidates$cx - widths[index] / 2
    candidates$right <- candidates$cx + widths[index] / 2
    candidates$bottom <- candidates$cy - heights[index] / 2
    candidates$top <- candidates$cy + heights[index] / 2
    candidates$inside <-
      candidates$left > usr[1] + 0.01 * x_range &
      candidates$right < usr[2] - 0.01 * x_range &
      candidates$bottom > usr[3] + 0.01 * y_range &
      candidates$top < usr[4] - 0.01 * y_range

    overlaps <- rep(FALSE, nrow(candidates))
    if (nrow(placed) > 0) {
      for (candidate_index in seq_len(nrow(candidates))) {
        candidate <- candidates[candidate_index, ]
        overlaps[candidate_index] <- any(
          candidate$left < placed[, 3] &
            candidate$right > placed[, 1] &
            candidate$bottom < placed[, 4] &
            candidate$top > placed[, 2]
        )
      }
    }

    valid <- which(candidates$inside & !overlaps)
    chosen <- if (length(valid) > 0) valid[1] else {
      which(candidates$inside)[1]
    }
    if (is.na(chosen)) {
      chosen <- 1
    }

    selected <- candidates[chosen, ]
    placed <- rbind(
      placed,
      c(selected$left, selected$bottom, selected$right, selected$top)
    )
    anchor_x <- min(max(x[index], selected$left), selected$right)
    anchor_y <- min(max(y[index], selected$bottom), selected$top)
    segments(
      x[index], y[index], anchor_x, anchor_y,
      col = seg_col, lwd = 0.65
    )
    text(selected$cx, selected$cy, labels[index], cex = cex, col = col)
  }
}

render <- function(stem, draw, width = 7.2, height = 6.4) {
  message("Rendering ", stem)
  png(
    file.path(outdir, paste0(stem, ".png")),
    width = width, height = height, units = "in", res = 300, bg = "white"
  )
  draw()
  dev.off()
  pdf(
    file.path(outdir, paste0(stem, ".pdf")),
    width = width, height = height, bg = "white"
  )
  draw()
  dev.off()
}

sample_master <- read_csv(file.path(analysis_root, "sample_master.csv"))
counts <- read_csv(file.path(results_root, "significant_counts.csv"))
shared <- read_csv(file.path(results_root, "shared_genes.csv"))
shared_sensitivity <- read_csv(
  file.path(sensitivity_results_root, "shared_genes.csv")
)
wgcna_summary <- read_csv(
  file.path(results_root, "wgcna", "wgcna_summary.csv")
)
module_assoc <- read_csv(
  file.path(results_root, "wgcna", "module_contrast_associations.csv")
)
cochlea_result <- read_csv(
  file.path(results_root, "limma", "cochlea_SP_vs_YC.csv")
)
ic_result <- read_csv(
  file.path(results_root, "limma", "inferior_colliculus_SP_vs_YC.csv")
)
age_ranges <- read_csv(file.path(results_root, "age_comparability.csv"))
cell_counts <- read_csv(
  file.path(four_matrix_dir, "GSE274279_four_matrix_cell_type_counts.csv")
)
umap <- read_csv(
  file.path(four_matrix_dir, "GSE274279_four_matrix_umap_coordinates.csv")
)
candidate_cell <- read_csv(
  file.path(
    four_matrix_dir,
    "GSE274279_four_matrix_candidate_expression_by_cell_type.csv"
  )
)
myeloid_age_summary <- read_csv(
  file.path(four_matrix_dir, "GSE274279_myeloid_age_descriptive_summary.csv")
)
myeloid_library_summary <- read_csv(
  file.path(
    four_matrix_dir,
    "GSE274279_myeloid_library_descriptive_summary.csv"
  )
)
single_cell_qc <- read_csv(
  file.path(four_matrix_dir, "GSE274279_four_matrix_QC_by_library.csv")
)
myeloid_umap <- read_csv(
  file.path(four_matrix_dir, "myeloid_four_matrix_umap_coordinates.csv")
)
myeloid_cluster_scores <- read_csv(
  file.path(
    four_matrix_dir,
    "myeloid_four_matrix_cluster_module_scores.csv"
  )
)
myeloid_age_proportions <- read_csv(
  file.path(
    four_matrix_dir,
    "myeloid_four_matrix_cluster_proportions_by_age.csv"
  )
)
myeloid_bootstrap <- read_csv(
  file.path(
    four_matrix_dir,
    "myeloid_four_matrix_cell_bootstrap_samples.csv"
  )
)
cellchat_summary <- read_csv(
  file.path(analysis_root, "cellchat_full", "analysis",
            "CellChat_full_age_summary.csv")
)
cellchat_focus <- read_csv(
  file.path(analysis_root, "cellchat_full", "analysis",
            "CellChat_full_focus_interactions.csv")
)

shared_genes <- c("Cd68", "Clec7a", "Fcgr3", "Lyz1", "Lyz2", "Mpeg1", "Ptprc")
shared_main_sp_yc <- shared[
  shared$contrast == "SP_vs_YC" & shared$gene %in% shared_genes, ,
  drop = FALSE
]
shared_main_sp_yc <- shared_main_sp_yc[
  match(shared_genes, shared_main_sp_yc$gene), ,
  drop = FALSE
]
shared_sensitivity_sp_yc <- shared_sensitivity[
  shared_sensitivity$contrast == "SP_vs_YC", ,
  drop = FALSE
]
shared_sensitivity_1602 <- read_csv(
  file.path(exclude_1602_results_root, "shared_genes.csv")
)
shared_sensitivity_1602 <- shared_sensitivity_1602[
  shared_sensitivity_1602$contrast == "SP_vs_YC", ,
  drop = FALSE
]
shared_sensitivity_both <- read_csv(
  file.path(exclude_both_results_root, "shared_genes.csv")
)
shared_sensitivity_both <- shared_sensitivity_both[
  shared_sensitivity_both$contrast == "SP_vs_YC", ,
  drop = FALSE
]
shared_union_genes <- unique(c(
  shared_main_sp_yc$gene,
  shared_sensitivity_sp_yc$gene
))
lookup_main_cochlea <- function(gene, column) {
  value <- cochlea_result[[column]][cochlea_result$gene == gene]
  if (length(value) == 0) NA_real_ else value[1]
}
shared_comparison <- data.frame(
  gene = shared_union_genes,
  logFC_cochlea_main = shared_main_sp_yc$logFC_cochlea[
    match(shared_union_genes, shared_main_sp_yc$gene)
  ],
  logFC_cochlea_sensitivity = shared_sensitivity_sp_yc$logFC_cochlea[
    match(shared_union_genes, shared_sensitivity_sp_yc$gene)
  ],
  adj.P.Val_cochlea_main = shared_main_sp_yc$adj.P.Val_cochlea[
    match(shared_union_genes, shared_main_sp_yc$gene)
  ],
  adj.P.Val_cochlea_sensitivity =
    shared_sensitivity_sp_yc$adj.P.Val_cochlea[
      match(shared_union_genes, shared_sensitivity_sp_yc$gene)
    ],
  logFC_inferior_colliculus =
    shared_sensitivity_sp_yc$logFC_inferior_colliculus[
      match(shared_union_genes, shared_sensitivity_sp_yc$gene)
    ],
  adj.P.Val_inferior_colliculus =
    shared_sensitivity_sp_yc$adj.P.Val_inferior_colliculus[
      match(shared_union_genes, shared_sensitivity_sp_yc$gene)
    ],
  shared_main = shared_union_genes %in% shared_genes,
  shared_sensitivity = shared_union_genes %in% shared_sensitivity_sp_yc$gene,
  stringsAsFactors = FALSE
)
missing_main_index <- is.na(shared_comparison$logFC_cochlea_main)
shared_comparison$logFC_cochlea_main[missing_main_index] <- vapply(
  shared_comparison$gene[missing_main_index],
  lookup_main_cochlea,
  numeric(1),
  column = "logFC"
)
missing_main_p_index <- is.na(shared_comparison$adj.P.Val_cochlea_main)
shared_comparison$adj.P.Val_cochlea_main[missing_main_p_index] <- vapply(
  shared_comparison$gene[missing_main_p_index],
  lookup_main_cochlea,
  numeric(1),
  column = "adj.P.Val"
)

figure1 <- function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE))
  par(family = "sans", mar = c(4.4, 4.5, 2.2, 1.0),
      mgp = c(2.5, 0.75, 0), las = 1)

  plot.new()
  plot.window(c(0, 1), c(0, 1))
  rect(0.02, 0.63, 0.31, 0.93, border = "#2E6F9E", lwd = 1.6, col = "#EDF4F8")
  rect(0.35, 0.63, 0.64, 0.93, border = "#A34D63", lwd = 1.6, col = "#F8EEF1")
  rect(0.68, 0.63, 0.98, 0.93, border = "#5C7D4E", lwd = 1.6, col = "#F0F5EC")
  text(0.165, 0.81, "GSE49543\nCochlea\n41 samples", font = 2, cex = 0.72)
  text(0.495, 0.81, "GSE49522\nInferior\ncolliculus\n39 samples", font = 2, cex = 0.66)
  text(
    0.83, 0.81,
    "GSE274279\nC57BL/6J\nCochlea snRNA-seq\n14,152 nuclei\n4 libraries",
    font = 2, cex = 0.56
  )
  rect(0.06, 0.16, 0.28, 0.40, border = "#777777", col = "#F7F7F7")
  rect(0.39, 0.16, 0.61, 0.40, border = "#777777", col = "#F7F7F7")
  rect(0.72, 0.16, 0.94, 0.40, border = "#777777", col = "#F7F7F7")
  text(0.17, 0.29, "Frozen limma\nand WGCNA", cex = 0.66)
  text(0.50, 0.29, "Intersection of\nseparately\nsignificant genes", cex = 0.58)
  text(
    0.83, 0.29,
    "Myeloid detection\n(four libraries)\nCellChat\n(three-library subset)",
    cex = 0.52
  )
  arrows(0.17, 0.62, 0.17, 0.41, length = 0.06, lwd = 1.2)
  arrows(0.83, 0.62, 0.83, 0.41, length = 0.06, lwd = 1.2)
  segments(c(0.10, 0.495), 0.62, c(0.10, 0.495), 0.55,
           lwd = 1.2, col = "#555555")
  segments(0.10, 0.55, 0.495, 0.55, lwd = 1.2, col = "#555555")
  arrows(0.50, 0.55, 0.50, 0.41, length = 0.06, lwd = 1.2)
  panel("A")

  count_matrix <- table(
    factor(
      ifelse(sample_master$tissue == "cochlea", "Cochlea", "Inferior colliculus"),
      levels = names(tissue_colors)
    ),
    factor(sample_master$cohort_group, levels = groups)
  )
  barplot(
    count_matrix, beside = TRUE, col = tissue_colors[rownames(count_matrix)],
    border = NA, ylab = "Observed samples", ylim = c(0, 20),
    names.arg = group_short[groups], cex.names = 0.78, cex.axis = 0.72
  )
  legend("topright", legend = rownames(count_matrix), fill = tissue_colors,
         border = NA, bty = "n", cex = 0.58)
  panel("B")

  sex_table <- aggregate(
    sample_id ~ tissue + sex,
    data = sample_master,
    FUN = length
  )
  names(sex_table)[names(sex_table) == "sample_id"] <- "n"
  sex_matrix <- xtabs(n ~ tissue + sex, data = sex_table)
  barplot(
    t(sex_matrix), beside = TRUE,
    col = c(Female = "#B55A78", Male = "#436B8E"),
    border = NA, ylab = "Observed samples",
    names.arg = c("Cochlea", "Inferior\ncolliculus"),
    cex.names = 0.68, cex.axis = 0.72
  )
  legend("topright", legend = c("Female", "Male"),
         fill = c("#B55A78", "#436B8E"), border = NA, bty = "n", cex = 0.62)
  panel("C")

  age_ranges$cohort_group <- factor(age_ranges$cohort_group, levels = groups)
  x <- seq_len(nrow(age_ranges))
  plot(
    x, age_ranges$age_mean_months, type = "n", xaxt = "n",
    xlim = c(0.6, 4.4), ylim = c(0, 35),
    ylab = "Published cohort age (months)", xlab = "",
    cex.axis = 0.72, cex.lab = 0.78
  )
  segments(
    x, age_ranges$mean_minus_one_sd_months,
    x, age_ranges$mean_plus_one_sd_months,
    col = group_colors[as.character(age_ranges$cohort_group)], lwd = 2
  )
  points(
    x, age_ranges$age_mean_months, pch = 21,
    bg = group_colors[as.character(age_ranges$cohort_group)],
    col = "white", cex = 1.3
  )
  axis(1, at = x, labels = group_short[as.character(age_ranges$cohort_group)],
       cex.axis = 0.78)
  rect(2.55, 28.7, 4.45, 31.1, border = NA, col = adjustcolor("#A93F55", 0.10))
  text(3.5, 33.0, "MP/SP intersection: 28.7-31.1 mo", cex = 0.54)
  panel("D")
}

volcano_panel <- function(data, color, label_genes) {
  significant <- data$significant
  plot(
    data$logFC, -log10(pmax(data$adj.P.Val, 1e-15)),
    pch = 21, cex = 0.45,
    bg = ifelse(significant, adjustcolor(color, 0.72), "#C9C9C9"),
    col = ifelse(significant, "white", NA),
    xlab = "log2 fold change (SP vs YC)",
    ylab = "-log10 adjusted P",
    cex.axis = 0.68, cex.lab = 0.76
  )
  abline(v = c(-1, 1), lty = 3, col = "#777777")
  abline(h = -log10(0.05), lty = 3, col = "#777777")
  selected <- data[data$gene %in% label_genes & significant, ]
  if (nrow(selected) > 0) {
    smart_labels(
      selected$logFC, -log10(pmax(selected$adj.P.Val, 1e-15)),
      selected$gene, cex = 0.50, col = color
    )
  }
}

figure2 <- function() {
  layout(matrix(1:3, nrow = 1), widths = c(1.25, 1, 1))
  par(family = "sans", mar = c(5.0, 4.2, 2.3, 0.8),
      mgp = c(2.5, 0.75, 0), las = 1)

  count_subset <- counts
  count_subset$label <- paste(count_subset$tissue, count_subset$contrast, sep = "\n")
  count_subset$label <- factor(count_subset$label, levels = count_subset$label)
  centers <- barplot(
    count_subset$n_significant,
    col = ifelse(count_subset$tissue == "cochlea", "#2E6F9E", "#A34D63"),
    border = NA, ylim = c(0, 55), ylab = "Genes meeting both criteria",
    names.arg = rep("", nrow(count_subset)), cex.axis = 0.68
  )
  text(
    centers, par("usr")[3] - 1.2,
    labels = count_subset$label, srt = 38, adj = 1,
    cex = 0.50, xpd = NA
  )
  text(centers, count_subset$n_significant, count_subset$n_significant,
       pos = 3, cex = 0.68)
  panel("A")

  par(mar = c(4.2, 4.2, 2.3, 0.8))
  volcano_panel(cochlea_result, "#2E6F9E", shared_genes)
  panel("B")
  volcano_panel(ic_result, "#A34D63", shared_genes)
  panel("C")
}

figure3 <- function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE))
  par(family = "sans", mar = c(4.3, 4.8, 2.2, 0.8),
      mgp = c(2.6, 0.75, 0), las = 1)

  y <- rev(seq_len(nrow(wgcna_summary)))
  plot(
    wgcna_summary$n_modules_excluding_grey, y,
    pch = 21, bg = tissue_colors[c("Cochlea", "Inferior colliculus")],
    col = "white", cex = 1.5, xlim = c(0, 36), yaxt = "n",
    xlab = "Biological modules", ylab = "",
    cex.axis = 0.70, cex.lab = 0.76
  )
  axis(2, at = y, labels = c("Cochlea", "Inferior\ncolliculus"),
       cex.axis = 0.68, las = 1)
  text(
    wgcna_summary$n_modules_excluding_grey, y,
    labels = paste0(
      wgcna_summary$n_modules_excluding_grey,
      " modules; power ", wgcna_summary$selected_power
    ),
    pos = 4, cex = 0.58
  )
  panel("A")

  module_panel <- function(tissue, color, main) {
    subset <- module_assoc[module_assoc$tissue == tissue, ]
    contrasts <- c("MA_vs_YC", "SP_vs_MP", "SP_vs_YC")
    contrast_colors <- c(
      MA_vs_YC = "#D9A441",
      SP_vs_MP = "#D66B3D",
      SP_vs_YC = color
    )
    plot(
      subset$logFC, -log10(pmax(subset$adj.P.Val, 1e-15)),
      type = "n", xlab = "Module eigengene log2 fold change",
      ylab = "-log10 adjusted P", main = main,
      ylim = c(0, max(-log10(pmax(subset$adj.P.Val, 1e-15))) * 1.10),
      cex.main = 0.80, cex.axis = 0.66, cex.lab = 0.72
    )
    for (contrast in contrasts) {
      current <- subset[subset$contrast == contrast, ]
      points(
        current$logFC, -log10(pmax(current$adj.P.Val, 1e-15)),
        pch = 21, col = "white",
        bg = adjustcolor(contrast_colors[contrast], 0.70), cex = 0.80
      )
    }
    abline(h = -log10(0.05), lty = 3, col = "#777777")
    legend(
      "bottomright", legend = c("MA vs YC", "SP vs MP", "SP vs YC"),
      pt.bg = contrast_colors, pch = 21, col = "white",
      bty = "n", cex = 0.48
    )
  }
  module_panel("cochlea", "#2E6F9E", "Cochlea")
  panel("B")
  module_panel("inferior_colliculus", "#A34D63", "Inferior colliculus")
  panel("C")

  sig_modules <- module_assoc[
    module_assoc$adj.P.Val < 0.05 & !module_assoc$is_grey_unassigned,
  ]
  sig_modules <- sig_modules[order(sig_modules$logFC), ]
  if (nrow(sig_modules) == 0) {
    plot.new()
    text(0.5, 0.5, "No module met FDR < 0.05")
  } else {
    y <- seq_len(nrow(sig_modules))
    plot(
      sig_modules$logFC, y, pch = 21,
      bg = ifelse(sig_modules$tissue == "cochlea", "#2E6F9E", "#A34D63"),
      col = "white", yaxt = "n", xlim = c(-3, 3),
      xlab = "Module eigengene log2 fold change",
      ylab = "", cex.axis = 0.66, cex.lab = 0.72
    )
    abline(v = 0, lty = 3, col = "#777777")
    axis(
      2, at = y,
      labels = paste(sig_modules$module, sig_modules$contrast, sep = "\n"),
      las = 2, cex.axis = 0.48
    )
  }
  panel("D")
}

figure3_intersection_sensitivity <- function() {
  layout(matrix(1:4, nrow = 2, byrow = TRUE))
  par(family = "sans", mar = c(4.3, 5.0, 2.4, 0.8),
      mgp = c(2.5, 0.75, 0), las = 1)

  shared_files <- c(
    Primary = file.path(results_root, "shared_genes.csv"),
    `Exclude\nGSM1201568` = file.path(sensitivity_results_root, "shared_genes.csv"),
    `Exclude\nGSM1201602` = file.path(exclude_1602_results_root, "shared_genes.csv"),
    `Exclude\nboth` = file.path(exclude_both_results_root, "shared_genes.csv")
  )
  shared_sets <- lapply(shared_files, function(path) {
    data <- read_csv(path)
    data$gene[data$contrast == "SP_vs_YC"]
  })
  primary_genes <- shared_sets[[1]]
  scenario_counts <- vapply(shared_sets, length, integer(1))
  message("  Figure 3A")
  centers <- barplot(
    scenario_counts,
    col = c("#3B6D8C", "#9B6A9F", "#C46A2D", "#5C7D4E"),
    border = NA, ylim = c(0, 10),
    ylab = "Intersection size",
    names.arg = names(scenario_counts),
    cex.names = 0.57, cex.axis = 0.68
  )
  text(centers, scenario_counts, scenario_counts, pos = 3, cex = 0.72)
  mtext("SP versus YC", side = 3, line = 0.8, cex = 0.72, font = 2)
  panel("A")

  membership_genes <- sort(unique(c(
    unlist(shared_sets),
    c("Cd68", "Clec7a", "Fcgr3", "Lyz1", "Lyz2", "Mpeg1", "Ptprc", "Tyrobp")
  )))
  membership <- sapply(shared_sets, function(genes) {
    membership_genes %in% genes
  })
  membership_order <- c(
    "Cd68", "Clec7a", "Fcgr3", "Lyz1",
    "Lyz2", "Mpeg1", "Ptprc", "Tyrobp"
  )
  membership <- membership[match(
    membership_order,
    membership_genes
  ), , drop = FALSE]
  rownames(membership) <- membership_order
  image(
    x = seq_len(ncol(membership)), y = seq_len(nrow(membership)),
    z = t(membership[nrow(membership):1, , drop = FALSE] * 1),
    col = c("#F1F1F1", "#7A5195"),
    breaks = c(-0.5, 0.5, 1.5),
    axes = FALSE, xlab = "", ylab = ""
  )
  axis(1, at = seq_len(ncol(membership)), labels = names(shared_sets),
       cex.axis = 0.50, las = 2)
  axis(2, at = rev(seq_len(nrow(membership))),
       labels = rownames(membership), cex.axis = 0.64, las = 2)
  for (i in seq_len(nrow(membership))) {
    for (j in seq_len(ncol(membership))) {
      text(j, rev(seq_len(nrow(membership)))[i],
           ifelse(membership[i, j], "Y", ""),
           cex = 0.48, col = "white", font = 2)
    }
  }
  box()
  message("  Figure 3B")
  mtext("Membership under each sample set", side = 3, line = 0.3, cex = 0.62)
  panel("B")

  message("  Figure 3C")
  frequency <- read_csv(
    file.path(loo_results_root, "leave_one_animal_out_gene_frequency.csv")
  )
  frequency <- frequency[match(rev(membership_order), frequency$gene), ]
  par(mar = c(4.3, 5.0, 2.4, 0.8))
  centers <- barplot(
    frequency$frequency * 100,
    horiz = TRUE,
    col = ifelse(frequency$frequency == 1, "#7A5195", "#D9A441"),
    border = NA,
    xlim = c(0, 108),
    names.arg = frequency$gene,
    xlab = "Models with gene in intersection (%)",
    cex.names = 0.62, cex.axis = 0.62, cex.lab = 0.68, las = 1
  )
  text(
    pmax(4.2, frequency$frequency * 100 - 4.2),
    centers,
    sprintf("%.1f", frequency$frequency * 100),
    cex = 0.53, adj = 0.5, col = "white", font = 2
  )
  mtext("Leave-one-animal-out frequency", side = 3, line = 0.8, cex = 0.68, font = 2)
  panel("C")

  message("  Figure 3D")
  loo_summary <- jsonlite::fromJSON(
    file.path(loo_results_root, "leave_one_animal_out_summary.json")
  )
  distribution <- loo_summary$shared_gene_count$distribution
  centers <- barplot(
    distribution$n_models,
    names.arg = distribution$shared_gene_count,
    col = "#3B6D8C",
    border = NA,
    ylim = c(0, max(distribution$n_models) * 1.2),
    xlab = "Intersection size in each model",
    ylab = "Number of models",
    cex.names = 0.68, cex.axis = 0.66, cex.lab = 0.68
  )
  text(centers, distribution$n_models, distribution$n_models, pos = 3, cex = 0.64)
  mtext(
    paste0(
      "41 models; range ", loo_summary$shared_gene_count$min,
      "-", loo_summary$shared_gene_count$max,
      "; median ", loo_summary$shared_gene_count$median
    ),
    side = 3, line = 0.8, cex = 0.60, font = 2
  )
  panel("D")
}

figure4 <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  layout(matrix(1:2, nrow = 1), widths = c(1.65, 1))
  par(family = "sans", mar = c(4.7, 5.1, 2.6, 0.8),
      mgp = c(2.6, 0.75, 0), las = 1)

  library_order <- c("3M", "12M", "24M", "24M_Cochlea_BU")
  library_labels <- c(
    "3M\nGSM8446833",
    "12M\nGSM8446834",
    "24M\nGSM8446835",
    "24M_BU\nGSM8446836"
  )
  measured_genes <- c("Cd68", "Fcgr3", "Lyz1", "Lyz2", "Mpeg1", "Ptprc")
  myeloid_library_summary$library <- factor(
    myeloid_library_summary$library,
    levels = library_order
  )
  myeloid_library_summary$gene <- factor(
    myeloid_library_summary$gene,
    levels = rev(measured_genes)
  )
  myeloid_library_summary$library_index <- match(
    as.character(myeloid_library_summary$library),
    library_order
  )

  plot.new()
  plot.window(
    xlim = c(0.45, 4.55),
    ylim = c(0.45, length(measured_genes) + 0.55),
    xaxs = "i", yaxs = "i"
  )
  abline(v = seq_along(library_order), col = "#ECECEC", lwd = 10)
  point_cex <- 0.65 + 2.7 * myeloid_library_summary$percent_expressing / 100
  expression_scale <- pmax(
    1,
    pmin(100, 1 + 99 * myeloid_library_summary$mean_log_expression / 1.7)
  )
  point_fill <- colorRampPalette(
    c("#F3E7C9", "#D7A33D", "#A93F55")
  )(100)[expression_scale]
  points(
    myeloid_library_summary$library_index,
    as.numeric(myeloid_library_summary$gene),
    pch = 21,
    cex = point_cex,
    bg = point_fill,
    col = "#563441",
    lwd = 0.7
  )
  text(
    myeloid_library_summary$library_index + 0.13,
    as.numeric(myeloid_library_summary$gene),
    sprintf("%.0f%%", myeloid_library_summary$percent_expressing),
    cex = 0.43,
    col = "#333333",
    adj = c(0, 0.5)
  )
  axis(1, at = seq_along(library_order), labels = library_labels,
       cex.axis = 0.55, line = 0.2)
  axis(
    2,
    at = seq_along(measured_genes),
    labels = rev(measured_genes),
    cex.axis = 0.68,
    las = 1
  )
  box(col = "#B9B9B9")
  mtext(
    "Size: nuclei detecting; color: mean normalized log expression",
    side = 1,
    line = 3.75,
    cex = 0.43
  )
  mtext(
    "C57BL/6J; four libraries; Clec7a not detected",
    side = 3,
    line = 0.75,
    cex = 0.53,
    font = 2
  )
  panel("A", line = 1.7)

  library_counts <- unique(
    myeloid_library_summary[
      myeloid_library_summary$gene == "Cd68",
      c("library", "cells_total")
    ]
  )
  library_counts$library <- factor(
    library_counts$library,
    levels = library_order
  )
  library_counts <- library_counts[order(library_counts$library), ]
  centers <- barplot(
    library_counts$cells_total,
    col = c("#4C78A8", "#D9A441", "#A34D63", "#7A5195"),
    border = NA,
    ylim = c(0, 105),
    ylab = "Marker-defined myeloid nuclei",
    names.arg = library_labels,
    cex.names = 0.46,
    cex.axis = 0.66,
    cex.lab = 0.72,
    las = 2
  )
  text(
    centers,
    library_counts$cells_total,
    library_counts$cells_total,
    pos = 3,
    cex = 0.68
  )
  mtext(
    paste0(
      "All four matrices; total myeloid = ",
      sum(library_counts$cells_total),
      " of ",
      sum(single_cell_qc$retained_nuclei),
      " nuclei"
    ),
    side = 3,
    line = 0.75,
    cex = 0.53,
    font = 2
  )
  panel("B", line = 1.7)
}

figureS2 <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)

  clip_limit <- 9
  layout(
    matrix(c(1, 2, 3, 4), nrow = 1),
    widths = c(4.25, 0.85, 4.2, 1.45)
  )
  par(family = "sans", mar = c(7.0, 5.2, 3.0, 0.4),
      mgp = c(2.5, 0.75, 0), las = 1)

  external_233798 <- read_csv(
    file.path(external_results_dir, "GSE233798_all_genes.csv")
  )
  external_ihc <- read_csv(
    file.path(external_results_dir, "GSE153882_IHC_all_genes.csv")
  )
  external_ohc <- read_csv(
    file.path(external_results_dir, "GSE153882_OHC_all_genes.csv")
  )
  external_stria <- read_csv(
    file.path(external_results_dir, "GSE154833_all_genes.csv")
  )

  get_effect <- function(data, gene, column) {
    row <- data[data$gene == gene, , drop = FALSE]
    if (nrow(row) == 0) NA_real_ else row[[column]][1]
  }
  get_optional <- function(data, gene, column) {
    if (!column %in% names(data)) return(NA_real_)
    get_effect(data, gene, column)
  }
  contexts <- c("GSE233798\ncochlea", "GSE153882\nIHC", "GSE153882\nOHC",
                "GSE154833\nstria")
  effect_raw <- sapply(shared_genes, function(gene) c(
    get_effect(external_233798, gene, "log2FC"),
    get_effect(external_ihc, gene, "log2FC"),
    get_effect(external_ohc, gene, "log2FC"),
    get_effect(external_stria, gene, "log2FC")
  ))
  p_values <- sapply(shared_genes, function(gene) c(
    get_optional(external_233798, gene, "p_value"),
    get_optional(external_ihc, gene, "p_value"),
    get_optional(external_ohc, gene, "p_value"),
    get_optional(external_stria, gene, "p_value")
  ))
  fdr_values <- sapply(shared_genes, function(gene) c(
    get_optional(external_233798, gene, "fdr_bh"),
    get_optional(external_ihc, gene, "fdr_bh"),
    get_optional(external_ohc, gene, "fdr_bh"),
    get_optional(external_stria, gene, "fdr_bh")
  ))
  effect <- effect_raw
  effect[is.finite(effect) & effect > clip_limit] <- clip_limit
  effect[is.finite(effect) & effect < -clip_limit] <- -clip_limit
  adjusted_p <- is.finite(fdr_values) & fdr_values < 0.05
  nominal_only <- is.finite(p_values) & p_values < 0.05 & !adjusted_p
  missing_matrix <- !is.finite(effect_raw)
  zero_matrix <- is.finite(effect_raw) & effect_raw == 0
  heat_colors <- colorRampPalette(c("#2E6F9E", "white", "#A93F55"))(257)
  heat_breaks <- seq(-clip_limit, clip_limit, length.out = 258)
  image(
    x = seq_len(ncol(effect)), y = seq_len(nrow(effect)),
    z = t(effect),
    col = heat_colors,
    breaks = heat_breaks,
    zlim = c(-clip_limit, clip_limit),
    axes = FALSE, xlab = "", ylab = ""
  )
  axis(1, at = seq_len(ncol(effect)), labels = colnames(effect),
       srt = 45, adj = 1, cex.axis = 0.80)
  axis(2, at = seq_len(nrow(effect)), labels = contexts,
       cex.axis = 0.76, las = 2)
  box()

  missing_cells <- which(!is.finite(effect), arr.ind = TRUE)
  if (nrow(missing_cells) > 0) {
    for (index in seq_len(nrow(missing_cells))) {
      row <- missing_cells[index, "row"]
      column <- missing_cells[index, "col"]
      rect(
        column - 0.5, row - 0.5, column + 0.5, row + 0.5,
        col = "#D9D9D9", border = "white", lwd = 1.2
      )
      segments(
        column - 0.32, row - 0.32, column + 0.32, row + 0.32,
        col = "#777777", lwd = 0.8
      )
      segments(
        column - 0.32, row + 0.32, column + 0.32, row - 0.32,
        col = "#777777", lwd = 0.8
      )
    }
  }

  zero_cells <- which(zero_matrix, arr.ind = TRUE)
  if (nrow(zero_cells) > 0) {
    points(
      zero_cells[, "col"], zero_cells[, "row"],
      pch = 1, cex = 0.85, lwd = 1.2, col = "black"
    )
  }
  adjusted_cells <- which(adjusted_p, arr.ind = TRUE)
  if (nrow(adjusted_cells) > 0) {
    text(
      adjusted_cells[, "col"] + 0.31,
      adjusted_cells[, "row"] + 0.31,
      "*", cex = 0.75, font = 2, col = "black"
    )
  }
  nominal_cells <- which(nominal_only, arr.ind = TRUE)
  if (nrow(nominal_cells) > 0) {
    points(
      nominal_cells[, "col"] + 0.31,
      nominal_cells[, "row"] + 0.31,
      pch = 3, cex = 0.58, lwd = 1.0, col = "#222222"
    )
  }

  mtext("External log2FC (fixed +/-9)", side = 3, line = 0.2, cex = 0.78)
  panel("A", cex = 0.98, line = 1.2)

  par(mar = c(7.0, 0.5, 3.0, 0.4), mgp = c(2.5, 0.75, 0), las = 1)
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(-clip_limit, clip_limit))
  rect(
    0, head(heat_breaks, -1), 1, tail(heat_breaks, -1),
    col = heat_colors, border = NA
  )
  rect(0, -clip_limit, 1, clip_limit, border = "#444444", lwd = 0.7)
  axis(
    4, at = c(-9, -6, -3, 0, 3, 6, 9),
    labels = c("-9", "-6", "-3", "0", "3", "6", "9"),
    cex.axis = 0.66, tick = FALSE, line = -0.6
  )
  mtext("log2FC", side = 3, line = 0.2, cex = 0.72)
  mtext(
    "Values beyond +/-9 are clipped",
    side = 1, line = 2.0, cex = 0.54
  )

  available <- colSums(is.finite(effect_raw))
  positive <- colSums(effect_raw > 0, na.rm = TRUE)
  negative <- colSums(effect_raw < 0, na.rm = TRUE)
  exact_zero <- colSums(effect_raw == 0, na.rm = TRUE)
  direction_counts <- data.frame(
    gene = shared_genes,
    available = available,
    positive = positive,
    negative = negative,
    exact_zero = exact_zero,
    unavailable = nrow(effect_raw) - available
  )
  write.csv(
    direction_counts,
    file.path(outdir, "figureS2_external_direction_counts.csv"),
    row.names = FALSE
  )

  external_long <- do.call(rbind, lapply(seq_along(contexts), function(index) {
    data.frame(
      gene = shared_genes,
      context = gsub("\n", " ", contexts[index]),
      log2FC_raw = effect_raw[index, ],
      log2FC_plotted = effect[index, ],
      p_value = p_values[index, ],
      fdr_bh = fdr_values[index, ],
      significance = ifelse(
        adjusted_p[index, ],
        "BH FDR < 0.05",
        ifelse(nominal_only[index, ], "Nominal P < 0.05 only", "Not significant")
      ),
      exact_zero = zero_matrix[index, ],
      missing = missing_matrix[index, ],
      stringsAsFactors = FALSE
    )
  }))
  write.csv(
    external_long,
    file.path(outdir, "figureS2_external_values.csv"),
    row.names = FALSE
  )

  y <- rev(seq_along(shared_genes))
  par(
    family = "sans", mar = c(7.0, 5.2, 3.0, 0.6),
    mgp = c(2.5, 0.75, 0), las = 1
  )
  plot(
    positive, y, pch = 21, bg = "#A93F55", col = "white",
    xlim = c(0, 4.25), ylim = c(0.5, length(shared_genes) + 0.5),
    yaxt = "n", type = "n",
    xlab = "Number of external contexts",
    ylab = "", cex.axis = 0.82, cex.lab = 0.88
  )
  positive_index <- positive > 0
  negative_index <- negative > 0
  zero_index <- exact_zero > 0
  points(
    positive[positive_index], y[positive_index],
    pch = 21, bg = "#A93F55", col = "white", cex = 1.08
  )
  points(
    negative[negative_index], y[negative_index],
    pch = 21, bg = "#2E6F9E", col = "white", cex = 1.08
  )
  points(
    exact_zero[zero_index], y[zero_index],
    pch = 1, cex = 1.1, lwd = 1.3, col = "#555555"
  )
  axis(2, at = y, labels = rev(shared_genes), cex.axis = 0.80, las = 2)
  legend("bottomright", legend = c("Positive", "Negative", "Exact zero"),
         pt.bg = c("#A93F55", "#2E6F9E", NA),
         pch = c(21, 21, 1), col = c("white", "white", "#555555"),
         bty = "n", cex = 0.66)
  mtext("Direction across available contexts", side = 3, line = 0.2, cex = 0.78)
  panel("B", cex = 0.98, line = 1.2)

  par(mar = c(7.0, 0.2, 3.0, 0.2))
  plot.new()
  plot.window(c(0, 1), c(0, 1))
  legend(
    "center",
    legend = c(
      "Asterisk: BH FDR < 0.05",
      "Plus: nominal P < 0.05 only",
      "Cross: missing; not imputed",
      "Circle: exact zero"
    ),
    title = "Significance and missingness",
    pch = c(42, 43, 4, 1),
    pt.cex = c(0.9, 0.8, 0.9, 0.85),
    col = c("#111111", "#222222", "#666666", "#111111"),
    bty = "n",
    cex = 0.62,
    xpd = NA
  )
}

figureS3 <- function() {
  for (extension in c("png", "pdf")) {
    candidates <- c(
      file.path(
        source_figure_dir,
        paste0("figureS3_cellchat_age_overview.", extension)
      ),
      file.path(
        source_figure_dir,
        paste0("figureS3_cellchat_full_database.", extension)
      ),
      file.path(
        archived_figure_dir,
        paste0("figure7_cellchat_age_overview.", extension)
      )
    )
    candidate <- candidates[file.exists(candidates)][1]
    if (is.na(candidate)) {
      stop("Supplementary Figure S3 source not found for extension: ", extension)
    }
    destination <- file.path(
      outdir,
      paste0("figureS3_cellchat_age_overview.", extension)
    )
    if (normalizePath(candidate, mustWork = TRUE) !=
        normalizePath(destination, mustWork = FALSE)) {
      file.copy(candidate, destination, overwrite = TRUE)
    }
  }
}

figureS4 <- function() {
  par(family = "sans", mfrow = c(1, 2), mar = c(4.5, 5.2, 2.5, 0.8),
      mgp = c(2.7, 0.8, 0), las = 1)

  ptprc <- cellchat_focus[cellchat_focus$interaction == "PTPRC_MRC1", ]
  ptprc <- ptprc[order(ptprc$age_months), ]
  plot(
    ptprc$age_months, ptprc$prob, type = "b",
    pch = 21, bg = "#A93F55", col = "#A93F55", lwd = 1.8,
    xaxt = "n", ylim = c(0, 0.11),
    xlab = "Age (months)", ylab = "CellChat communication probability",
    cex.axis = 0.92, cex.lab = 0.92
  )
  axis(1, at = c(3, 12, 24), labels = c("3M", "12M", "24M"), cex.axis = 0.95)
  text(ptprc$age_months, ptprc$prob, sprintf("%.3f", ptprc$prob),
       pos = c(4, 3, 2), cex = 0.78)
  mtext("PTPRC-MRC1\nMacrophages/microglia -> macrophages/microglia",
        side = 3, line = 0.2, cex = 0.76)
  panel("A", line = 1.2)

  app <- cellchat_focus[
    cellchat_focus$interaction == "APP_CD74" &
      cellchat_focus$target == "Macrophages/Microglia",
  ]
  app <- app[order(app$age_months, app$prob), ]
  par(mar = c(4.5, 9.0, 2.5, 0.8))
  age_colors <- c(`3` = "#4C78A8", `12` = "#D9A441", `24` = "#A93F55")
  y_positions <- seq_len(nrow(app))
  plot(
    app$prob, y_positions, pch = 21,
    bg = age_colors[as.character(app$age_months)], col = "white",
    yaxt = "n", xlim = c(0, 0.17),
    xlab = "CellChat communication probability",
    ylab = "", cex.axis = 0.88, cex.lab = 0.88
  )
  axis(
    2, at = y_positions,
    labels = paste(
      c(
        `Inner hair cells` = "IHC",
        `Outer hair cells` = "OHC",
        `Supporting cells` = "Supporting",
        Fibrocytes = "Fibrocytes",
        `Macrophages/Microglia` = "Macrophage/microglia",
        `Glia/Schwann` = "Glia/Schwann",
        `Spiral ganglion neurons` = "SGN"
      )[app$source],
      paste0(app$age_months, "M"),
      sep = " | "
    ),
    las = 2, cex.axis = 0.68
  )
  legend("bottomright", legend = c("3M", "12M", "24M"),
         pt.bg = unname(age_colors), pch = 21, col = "white",
         bty = "n", cex = 0.70)
  mtext("Inferred APP-CD74 interactions to myeloid targets",
        side = 3, line = 0.2, cex = 0.76)
  panel("B", line = 1.2)
}

figure5 <- function() {
  layout(matrix(1:2, nrow = 1), widths = c(1.05, 1))
  par(family = "sans", mar = c(4.5, 4.8, 2.3, 0.8),
      mgp = c(2.5, 0.75, 0), las = 1)

  sex_counts <- as.data.frame(
    table(
      tissue = factor(
        ifelse(
          sample_master$tissue == "cochlea",
          "Cochlea", "Inferior colliculus"
        ),
        levels = c("Cochlea", "Inferior colliculus")
      ),
      cohort_group = factor(sample_master$cohort_group, levels = groups),
      sex = factor(sample_master$sex, levels = c("Female", "Male"))
    )
  )
  sex_matrix <- xtabs(Freq ~ cohort_group + sex + tissue, data = sex_counts)
  sex_by_group <- matrix(
    c(
      sex_matrix[, "Female", "Cochlea"],
      sex_matrix[, "Male", "Cochlea"],
      sex_matrix[, "Female", "Inferior colliculus"],
      sex_matrix[, "Male", "Inferior colliculus"]
    ),
    nrow = 4, byrow = TRUE
  )
  combined_plot <- rbind(
    c(sex_by_group[1, ], sex_by_group[3, ]),
    c(sex_by_group[2, ], sex_by_group[4, ])
  )
  colnames(combined_plot) <- paste(rep(c("C", "IC"), each = 4), groups, sep = "-")
  rownames(combined_plot) <- c("Female", "Male")
  centers <- barplot(
    combined_plot,
    col = c("#B55A78", "#436B8E"),
    border = NA,
    ylim = c(0, max(colSums(combined_plot)) * 1.28),
    ylab = "Observed samples",
    names.arg = rep(groups, 2),
    cex.names = 0.66, cex.axis = 0.66, cex.lab = 0.74
  )
  text(
    centers, colSums(combined_plot) + 0.45,
    paste0(
      combined_plot[1, ], "F/", combined_plot[2, ], "M"
    ),
    cex = 0.46, xpd = NA,
    col = ifelse(seq_along(centers) == 8, "#A93F55", "#333333"),
    font = ifelse(seq_along(centers) == 8, 2, 1)
  )
  mtext(
    "Cochlea                                         Inferior colliculus",
    side = 3, line = 1.65, cex = 0.58, font = 2
  )
  legend(
    "topright", legend = c("Female", "Male"),
    fill = c("#B55A78", "#436B8E"), border = NA, bty = "n", cex = 0.58
  )
  usr <- par("usr")
  text(
    usr[1] + 0.18, usr[4] - 0.10, "A",
    adj = c(0, 1), cex = 1.18, font = 2, xpd = NA
  )

  comparison <- read_csv(
    file.path(
      no_sex_results_root, "no_sex_vs_sex_adjusted_all_genes.csv"
    )
  )
  comparison <- comparison[comparison$contrast == "SP_vs_YC", ]
  tissue_plot_colors <- c(
    cochlea = "#2E6F9E", inferior_colliculus = "#A34D63"
  )
  plot(
    comparison$logFC_sex_adjusted,
    comparison$logFC_no_sex,
    type = "n",
    xlab = "Sex-adjusted SP-versus-YC log2FC",
    ylab = "No-sex sensitivity log2FC",
    cex.axis = 0.68, cex.lab = 0.76
  )
  abline(a = 0, b = 1, col = "#777777", lty = 2, lwd = 1)
  abline(h = 0, v = 0, col = "#DDDDDD", lwd = 0.8)
  for (tissue in names(tissue_plot_colors)) {
    current <- comparison$tissue == tissue
    points(
      comparison$logFC_sex_adjusted[current],
      comparison$logFC_no_sex[current],
      pch = ifelse(tissue == "cochlea", 16, 17),
      col = adjustcolor(tissue_plot_colors[tissue], 0.38),
      cex = 0.28
    )
  }
  primary <- comparison[
    comparison$gene %in% shared_genes & comparison$significant_sex_adjusted,
  ]
  points(
    primary$logFC_sex_adjusted, primary$logFC_no_sex,
    pch = 21, bg = "#F2D8B8", col = "#6B3B16", cex = 0.95, lwd = 0.8
  )
  concordance <- read_csv(
    file.path(no_sex_results_root, "no_sex_sp_vs_yc_effect_concordance.csv")
  )
  intersection_summary <- read_csv(
    file.path(no_sex_results_root, "no_sex_intersection_summary.csv")
  )
  sp_summary <- intersection_summary[
    intersection_summary$contrast == "SP_vs_YC", ,
    drop = FALSE
  ]
  legend(
    "topleft",
    legend = c(
      "Cochlea", "Inferior colliculus", "7 SP-vs-YC intersection genes"
    ),
    col = c(
      tissue_plot_colors[c("cochlea", "inferior_colliculus")],
      "#6B3B16"
    ),
    pch = c(16, 17, 21), pt.bg = c(NA, NA, "#F2D8B8"),
    bty = "n", cex = 0.54
  )
  usr <- par("usr")
  text(
    usr[1] + 0.35, usr[3] + 0.35,
    paste0(
      "SP vs YC intersection: ",
      sp_summary$sex_adjusted_intersection_n, " with sex; ",
      sp_summary$no_sex_intersection_n, " without sex\n",
      "Spearman rho: cochlea ",
      sprintf("%.3f", concordance$spearman_logFC[
        concordance$tissue == "cochlea"
      ]),
      "; inferior colliculus ",
      sprintf("%.3f", concordance$spearman_logFC[
        concordance$tissue == "inferior_colliculus"
      ])
    ),
    adj = c(0, 0), cex = 0.46, font = 2, col = "#3B3B3B"
  )
  mtext("Sex-adjustment sensitivity", side = 3, line = 0.75,
        cex = 0.62, font = 2)
  panel("B")
}

figureS5 <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  layout(matrix(1:4, nrow = 2, byrow = TRUE))
  par(family = "sans", mar = c(4.2, 4.6, 2.6, 0.9),
      mgp = c(2.5, 0.75, 0), las = 1)

  state_colors <- c(
    Myeloid_C0 = "#4C78A8",
    Myeloid_C1 = "#D66B3D"
  )
  state_levels <- names(state_colors)

  plot(
    myeloid_umap$umap_1,
    myeloid_umap$umap_2,
    type = "n",
    xlab = "UMAP1",
    ylab = "UMAP2",
    cex.axis = 0.64,
    cex.lab = 0.72
  )
  for (state in state_levels) {
    selected <- myeloid_umap$myeloid_state == state
    points(
      myeloid_umap$umap_1[selected],
      myeloid_umap$umap_2[selected],
      pch = 16,
      col = adjustcolor(state_colors[[state]], 0.78),
      cex = 0.55
    )
  }
  legend(
    "topright",
    legend = state_levels,
    col = state_colors,
    pch = 16,
    bty = "n",
    cex = 0.56
  )
  mtext(
    "245 marker-defined myeloid nuclei",
    side = 3,
    line = 0.7,
    cex = 0.62,
    font = 2
  )
  panel("A")

  age_columns <- c("3", "12", "24")
  proportion_matrix <- xtabs(
    proportion ~ myeloid_state + age_months,
    data = myeloid_age_proportions
  )
  proportion_matrix <- proportion_matrix[state_levels, age_columns, drop = FALSE]
  barplot(
    proportion_matrix,
    col = state_colors[rownames(proportion_matrix)],
    border = NA,
    ylim = c(0, 1.04),
    ylab = "Proportion of myeloid nuclei",
    xlab = "Age (months)",
    names.arg = c("3M", "12M", "24M"),
    cex.names = 0.68,
    cex.axis = 0.63,
    cex.lab = 0.68
  )
  legend(
    "bottomright",
    legend = state_levels,
    fill = state_colors,
    border = NA,
    bty = "n",
    cex = 0.52
  )
  mtext(
    "Exploratory states; age and library are confounded",
    side = 3,
    line = 0.7,
    cex = 0.56
  )
  panel("B")

  module_columns <- setdiff(
    names(myeloid_cluster_scores),
    c("myeloid_state", "nuclei")
  )
  score_matrix <- as.matrix(
    myeloid_cluster_scores[
      match(state_levels, myeloid_cluster_scores$myeloid_state),
      module_columns,
      drop = FALSE
    ]
  )
  rownames(score_matrix) <- state_levels
  displayed_scores <- pmax(pmin(score_matrix, 0.8), -0.8)
  module_labels <- sub("_score1$", "", module_columns)
  module_labels <- gsub("_", "-", module_labels)
  palette_values <- colorRampPalette(
    c("#2F5D8A", "#F7F4EA", "#A93F55")
  )(101)
  image(
    x = seq_along(module_columns),
    y = seq_len(nrow(displayed_scores)),
    z = t(displayed_scores[nrow(displayed_scores):1, , drop = FALSE]),
    col = palette_values,
    zlim = c(-0.8, 0.8),
    axes = FALSE,
    xlab = "",
    ylab = ""
  )
  axis(
    1,
    at = seq_along(module_columns),
    labels = module_labels,
    cex.axis = 0.53,
    las = 2
  )
  axis(
    2,
    at = seq_len(nrow(displayed_scores)),
    labels = rev(state_levels),
    cex.axis = 0.62,
    las = 1
  )
  for (row_index in seq_len(nrow(displayed_scores))) {
    for (column_index in seq_along(module_columns)) {
      text(
        column_index,
        nrow(displayed_scores) - row_index + 1,
        sprintf(
          "%.2f",
          score_matrix[row_index, column_index]
        ),
        cex = 0.42,
        col = "#222222"
      )
    }
  }
  box(col = "#B9B9B9")
  mtext(
    "Mean module score (display clipped to +/-0.8)",
    side = 3,
    line = 0.5,
    cex = 0.54,
    font = 2
  )
  panel("C")

  module_labels_raw <- sub("_score1$", "", module_columns)
  myeloid_bootstrap$module <- factor(
    myeloid_bootstrap$module,
    levels = module_labels_raw
  )
  boxplot(
    delta_24m_minus_3m ~ module,
    data = myeloid_bootstrap,
    outline = FALSE,
    col = "#E4E8EB",
    border = "#4C555C",
    las = 2,
    ylab = "24M minus 3M module score",
    cex.axis = 0.48,
    cex.lab = 0.62
  )
  abline(h = 0, lty = 2, col = "#A93F55", lwd = 1.2)
  mtext(
    "500 nucleus-level bootstrap iterations",
    side = 3,
    line = 0.7,
    cex = 0.58,
    font = 2
  )
  mtext(
    "Descriptive for sampled nuclei; not an animal-level confidence interval",
    side = 1,
    line = 2.70,
    cex = 0.43
  )
  panel("D")
}

figureS6 <- function() {
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  layout(matrix(1:2, nrow = 1), widths = c(1.55, 1))
  par(family = "sans", mar = c(4.4, 4.8, 2.5, 0.8),
      mgp = c(2.5, 0.75, 0), las = 1)

  cell_type_colors <- c(
    `Supporting cells` = "#6F4E7C",
    `Inner hair cells` = "#C7A76B",
    `Outer hair cells` = "#D66B3D",
    `Spiral ganglion neurons` = "#4C78A8",
    `Stria vascularis` = "#3A8D8F",
    Endothelial = "#7E9D45",
    `Macrophages/Microglia` = "#A93F55"
  )
  plot(
    umap$umap_1,
    umap$umap_2,
    type = "n",
    xlab = "UMAP1",
    ylab = "UMAP2",
    cex.axis = 0.62,
    cex.lab = 0.70
  )
  for (cell_type in names(cell_type_colors)) {
    selected <- umap$cell_type == cell_type
    points(
      umap$umap_1[selected],
      umap$umap_2[selected],
      pch = 16,
      col = adjustcolor(cell_type_colors[[cell_type]], 0.68),
      cex = ifelse(cell_type == "Macrophages/Microglia", 0.62, 0.25)
    )
  }
  legend(
    "topright",
    legend = names(cell_type_colors),
    col = cell_type_colors,
    pch = 16,
    bty = "n",
    cex = 0.45
  )
  mtext(
    "C57BL/6J; all four GSE274279 libraries",
    side = 3,
    line = 0.7,
    cex = 0.60,
    font = 2
  )
  panel("A")

  library_order <- c("3M", "12M", "24M", "24M_Cochlea_BU")
  library_labels <- c(
    "3M\nGSM8446833",
    "12M\nGSM8446834",
    "24M\nGSM8446835",
    "24M_BU\nGSM8446836"
  )
  cell_counts$library <- factor(
    cell_counts$library,
    levels = library_order
  )
  cell_counts$cell_type <- factor(
    cell_counts$cell_type,
    levels = names(cell_type_colors)
  )
  count_matrix <- xtabs(
    nuclei ~ cell_type + library,
    data = cell_counts
  )
  if (!identical(colnames(count_matrix), library_order)) {
    stop("Supplementary Figure S6B library order does not match labels.")
  }
  observed_totals <- colSums(count_matrix)
  expected_totals <- c("3M" = 2916, "12M" = 5532, "24M" = 2954,
                       "24M_Cochlea_BU" = 2750)
  if (!isTRUE(all.equal(as.numeric(observed_totals), as.numeric(expected_totals)))) {
    stop("Supplementary Figure S6B library totals changed: ",
         paste(names(observed_totals), observed_totals, sep = "=",
               collapse = ", "))
  }
  names(library_labels) <- library_order
  bar_midpoints <- barplot(
    count_matrix,
    col = cell_type_colors[rownames(count_matrix)],
    border = NA,
    ylim = c(0, max(colSums(count_matrix)) * 1.12),
    ylab = "Nuclei",
    xlab = "10x library",
    names.arg = unname(library_labels[colnames(count_matrix)]),
    cex.names = 0.49,
    cex.axis = 0.62,
    cex.lab = 0.68
  )
  text(
    bar_midpoints,
    observed_totals + max(observed_totals) * 0.025,
    labels = format(observed_totals, big.mark = ",", trim = TRUE),
    cex = 0.56,
    font = 2
  )
  mtext(
    "Descriptive cell-type counts; one library per age except 24M",
    side = 3,
    line = 0.7,
    cex = 0.53
  )
  panel("B")
}

render_requested <- function(stem) {
  requested <- trimws(strsplit(Sys.getenv("V14_RENDER_TARGETS"), ",", fixed = TRUE)[[1]])
  !nzchar(Sys.getenv("V14_RENDER_TARGETS")) || stem %in% requested
}

if (identical(Sys.getenv("V14_RENDER_FIGURE1_ONLY"), "true")) {
  render("figure1_study_audit", figure1)
} else {
  if (render_requested("figure1_study_audit")) {
    render("figure1_study_audit", figure1)
  }
  if (render_requested("figure2_frozen_limma")) {
    render("figure2_frozen_limma", figure2, width = 10.2, height = 4.7)
  }
  if (render_requested("figure3_intersection_sensitivity")) {
    render(
      "figure3_intersection_sensitivity",
      figure3_intersection_sensitivity,
      width = 8.0,
      height = 7.2
    )
  }
  if (render_requested("figure4_cellular_context")) {
    render(
      "figure4_cellular_context",
      figure4,
      width = 10.0,
      height = 5.0
    )
  }
  if (render_requested("figure5_model_identifiability")) {
    render(
      "figure5_model_identifiability",
      figure5,
      width = 10.4,
      height = 4.9
    )
  }
  if (render_requested("figureS1_wgcna")) {
    render("figureS1_wgcna", figure3)
  }
  if (render_requested("figureS2_external_context")) {
    render(
      "figureS2_external_context",
      figureS2,
    width = 10.4,
    height = 4.8
    )
  }
  if (render_requested("figureS3_cellchat_age_overview")) {
    figureS3()
  }
  if (render_requested("figureS4_cellchat_focus")) {
    render(
      "figureS4_cellchat_focus",
      figureS4,
      width = 9.2,
      height = 4.8
    )
  }
  if (render_requested("figureS5_myeloid_subclustering")) {
    render(
      "figureS5_myeloid_subclustering",
      figureS5,
      width = 8.2,
      height = 7.0
    )
  }
  if (render_requested("figureS6_full_atlas")) {
    render(
      "figureS6_full_atlas",
      figureS6,
      width = 10.0,
      height = 5.0
    )
  }
}

cat("V14 figures written to", outdir, "\n")
