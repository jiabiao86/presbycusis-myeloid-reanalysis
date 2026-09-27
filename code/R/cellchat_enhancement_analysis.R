#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(Matrix))

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
out <- file.path(root, "cellchat_enhancement")
data_dir <- "/tmp/ear_study_learning/GSE274279"

cell_order <- c(
  "Supporting cells",
  "Fibrocytes",
  "Glia/Schwann",
  "Outer hair cells",
  "Inner hair cells",
  "Macrophages/Microglia",
  "Spiral ganglion neurons"
)
age_files <- list(
  "3" = c("3M_matrix.mtx.gz", "3M_barcodes.tsv.gz", "-0"),
  "12" = c("12M_matrix.mtx.gz", "12M_barcodes.tsv.gz", "-1"),
  "24" = c("24M_matrix.mtx.gz", "24M_barcodes.tsv.gz", "-2")
)

interactions <- read.csv(
  file.path(out, "cellchatdb_mouse_interactions.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
complexes <- read.csv(
  file.path(out, "cellchatdb_mouse_complexes.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
complex_map <- setNames(
  lapply(seq_len(nrow(complexes)), function(index) {
    values <- unlist(complexes[index, c("subunit_1", "subunit_2", "subunit_3", "subunit_4")])
    unique(values[!is.na(values) & nzchar(values)])
  }),
  complexes$complex_name
)

entity_subunits <- function(entity) {
  entity <- as.character(entity)
  if (entity %in% names(complex_map)) {
    return(complex_map[[entity]])
  }
  entity
}

features <- read.delim(
  gzfile(file.path(data_dir, "features.tsv.gz")),
  header = FALSE,
  stringsAsFactors = FALSE
)
gene_symbols <- features[[2]]

umap <- read.csv(
  file.path(root, "GSE274279_umap_coordinates.csv"),
  row.names = 1,
  check.names = FALSE
)

required_genes <- unique(unlist(lapply(seq_len(nrow(interactions)), function(index) {
  c(
    entity_subunits(interactions$ligand[index]),
    entity_subunits(interactions$receptor[index])
  )
})))

normalize_matrix <- function(matrix) {
  totals <- Matrix::colSums(matrix)
  totals[totals <= 0] <- 1
  matrix <- t(t(matrix) * (1e4 / totals))
  matrix@x <- log1p(matrix@x)
  matrix
}

aggregate_expression <- function(matrix, labels) {
  labels <- factor(labels, levels = cell_order)
  present <- levels(droplevels(labels))
  means <- matrix(0, nrow = nrow(matrix), ncol = length(present))
  percents <- means
  rownames(means) <- rownames(matrix)
  colnames(means) <- present
  rownames(percents) <- rownames(matrix)
  colnames(percents) <- present
  for (cell_type in present) {
    mask <- labels == cell_type
    subset <- matrix[, mask, drop = FALSE]
    means[, cell_type] <- Matrix::rowMeans(subset)
    percents[, cell_type] <- Matrix::rowSums(subset > 0) / sum(mask) * 100
  }
  list(means = means, percents = percents)
}

score_rows <- list()
for (age in names(age_files)) {
  matrix_name <- age_files[[age]][1]
  barcode_name <- age_files[[age]][2]
  suffix <- age_files[[age]][3]
  barcodes <- readLines(gzfile(file.path(data_dir, barcode_name)))
  index <- paste0(barcodes, suffix)

  age_umap <- umap[umap$age_months == as.numeric(age), , drop = FALSE]
  age_umap <- age_umap[index[index %in% rownames(age_umap)], , drop = FALSE]
  keep_cells <- match(rownames(age_umap), index)

  matrix <- Matrix::readMM(gzfile(file.path(data_dir, matrix_name)))
  matrix <- as(matrix, "dgCMatrix")
  rownames(matrix) <- gene_symbols
  matrix <- matrix[rownames(matrix) %in% required_genes, keep_cells, drop = FALSE]
  matrix <- normalize_matrix(matrix)
  stats <- aggregate_expression(matrix, age_umap$cell_type)

  age_rows <- vector("list", nrow(interactions))
  counter <- 0L
  for (interaction_index in seq_len(nrow(interactions))) {
    interaction <- interactions[interaction_index, ]
    ligands <- entity_subunits(interaction$ligand)
    receptors <- entity_subunits(interaction$receptor)
    if (!all(c(ligands, receptors) %in% rownames(stats$means))) {
      next
    }
    ligand_mean <- apply(stats$means[ligands, , drop = FALSE], 2, min)
    receptor_mean <- apply(stats$means[receptors, , drop = FALSE], 2, min)
    ligand_percent <- apply(stats$percents[ligands, , drop = FALSE], 2, min)
    receptor_percent <- apply(stats$percents[receptors, , drop = FALSE], 2, min)

    for (sender in colnames(stats$means)) {
      for (receiver in colnames(stats$means)) {
        score <- sqrt(max(ligand_mean[sender], 0) * max(receptor_mean[receiver], 0)) *
          sqrt(max(ligand_percent[sender], 0) / 100 * max(receptor_percent[receiver], 0) / 100)
        if (score <= 0) {
          next
        }
        counter <- counter + 1L
        age_rows[[counter]] <- data.frame(
          age_months = as.numeric(age),
          sender = sender,
          receiver = receiver,
          interaction_name = interaction$interaction_name,
          pathway_name = interaction$pathway_name,
          ligand = interaction$ligand,
          receptor = interaction$receptor,
          annotation = interaction$annotation,
          ligand_mean = ligand_mean[sender],
          receptor_mean = receptor_mean[receiver],
          ligand_percent = ligand_percent[sender],
          receptor_percent = receptor_percent[receiver],
          communication_score = score,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  score_rows[[age]] <- do.call(rbind, age_rows[seq_len(counter)])
  cat("Age", age, "M:", nrow(age_umap), "cells,", counter, "interaction scores\n")
  rm(matrix)
  gc()
}

all_scores <- do.call(rbind, score_rows)
write.csv(
  all_scores,
  file.path(out, "cellchatdb_communication_scores.csv"),
  row.names = FALSE
)

pathway_scores <- aggregate(
  communication_score ~ age_months + sender + receiver + pathway_name,
  data = all_scores,
  FUN = sum
)
write.csv(
  pathway_scores,
  file.path(out, "cellchatdb_pathway_scores.csv"),
  row.names = FALSE
)

key_columns <- c("sender", "receiver", "interaction_name", "pathway_name", "ligand", "receptor", "annotation")
comparison <- aggregate(
  communication_score ~ age_months + sender + receiver + interaction_name + pathway_name + ligand + receptor + annotation,
  data = all_scores,
  FUN = mean
)
age3 <- comparison[comparison$age_months == 3, c(key_columns, "communication_score")]
age24 <- comparison[comparison$age_months == 24, c(key_columns, "communication_score")]
changed <- merge(age3, age24, by = key_columns, all = TRUE, suffixes = c("_3m", "_24m"))
changed$communication_score_3m[is.na(changed$communication_score_3m)] <- 0
changed$communication_score_24m[is.na(changed$communication_score_24m)] <- 0
changed$delta_24m_minus_3m <- changed$communication_score_24m - changed$communication_score_3m
changed$log2FC_24m_vs_3m <- log2(
  (changed$communication_score_24m + 1e-6) /
    (changed$communication_score_3m + 1e-6)
)
changed$involving_myeloid <- changed$sender == "Macrophages/Microglia" |
  changed$receiver == "Macrophages/Microglia"
changed <- changed[order(changed$delta_24m_minus_3m, decreasing = TRUE), ]
write.csv(
  changed,
  file.path(out, "cellchatdb_interaction_changes.csv"),
  row.names = FALSE
)
write.csv(
  head(changed[changed$involving_myeloid, ], 100),
  file.path(out, "cellchatdb_top_myeloid_changes.csv"),
  row.names = FALSE
)

cat("Wrote CellChatDB communication tables to", out, "\n")
