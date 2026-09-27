#!/usr/bin/env Rscript

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
lib <- "/Users/jijiabiao/Desktop/耳聋课题/Rlibs"
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
})

data_dir <- "/tmp/ear_study_learning/GSE274279"
out <- file.path(root, "myeloid_subclustering")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

age_files <- list(
  "3" = c("3M_matrix.mtx.gz", "3M_barcodes.tsv.gz", "-0"),
  "12" = c("12M_matrix.mtx.gz", "12M_barcodes.tsv.gz", "-1"),
  "24" = c("24M_matrix.mtx.gz", "24M_barcodes.tsv.gz", "-2")
)

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

matrices <- list()
metadata <- list()
for (age in names(age_files)) {
  matrix_name <- age_files[[age]][1]
  barcode_name <- age_files[[age]][2]
  suffix <- age_files[[age]][3]
  barcodes <- readLines(gzfile(file.path(data_dir, barcode_name)))
  barcode_index <- paste0(barcodes, suffix)
  age_meta <- umap[umap$age_months == as.numeric(age), , drop = FALSE]
  keep <- barcode_index %in% rownames(age_meta)
  barcode_index <- barcode_index[keep]
  barcodes <- barcodes[keep]
  myeloid_keep <- age_meta[barcode_index, "cell_type"] == "Macrophages/Microglia"
  barcodes <- barcodes[myeloid_keep]
  barcode_index <- barcode_index[myeloid_keep]

  counts <- Matrix::readMM(gzfile(file.path(data_dir, matrix_name)))
  counts <- as(counts, "CsparseMatrix")
  counts <- counts[, keep, drop = FALSE]
  counts <- counts[, myeloid_keep, drop = FALSE]
  rownames(counts) <- gene_symbols
  colnames(counts) <- barcode_index
  counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]

  matrices[[age]] <- counts
  metadata[[age]] <- data.frame(
    cell_id = barcode_index,
    age_months = as.numeric(age),
    broad_cell_type = "Macrophages/Microglia",
    row.names = barcode_index,
    stringsAsFactors = FALSE
  )
}

common_genes <- Reduce(intersect, lapply(matrices, rownames))
combined_counts <- do.call(cbind, lapply(matrices, function(matrix) matrix[common_genes, , drop = FALSE]))
combined_metadata <- do.call(rbind, metadata)
rownames(combined_metadata) <- colnames(combined_counts)

seurat <- CreateSeuratObject(
  counts = combined_counts,
  meta.data = combined_metadata,
  min.cells = 1,
  min.features = 50
)
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^mt-")
seurat <- NormalizeData(seurat, verbose = FALSE)
seurat <- FindVariableFeatures(seurat, nfeatures = 1000, verbose = FALSE)
seurat <- ScaleData(
  seurat,
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose = FALSE
)
seurat <- RunPCA(seurat, npcs = 20, verbose = FALSE)
seurat <- FindNeighbors(seurat, dims = 1:15, verbose = FALSE)
seurat <- FindClusters(seurat, resolution = 0.4, verbose = FALSE)
seurat <- RunUMAP(seurat, dims = 1:15, verbose = FALSE)

module_sets <- list(
  Homeostatic = c("P2ry12", "Tmem119", "Cx3cr1", "Hexb"),
  MHCII = c("Cd74", "H2-Aa", "H2-Eb1", "H2-Ab1"),
  Complement = c("C1qa", "C1qb", "C1qc"),
  Phagocytic = c("Mertk", "Axl", "Gas6", "Cd68", "Mpeg1"),
  Interferon = c("Ifit1", "Ifit3", "Irf7", "Stat1"),
  APP_CD74 = c("App", "Cd74"),
  PTPRC_MRC1 = c("Ptprc", "Mrc1")
)

for (module in names(module_sets)) {
  genes <- intersect(module_sets[[module]], rownames(seurat))
  if (length(genes) >= 2) {
    seurat <- AddModuleScore(
      seurat,
      features = list(genes),
      name = paste0(module, "_score"),
      assay = "RNA"
    )
  }
}

module_columns <- paste0(names(module_sets), "_score1")
cluster_scores <- seurat@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(
    across(all_of(module_columns), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
write.csv(
  cluster_scores,
  file.path(out, "myeloid_cluster_module_scores.csv"),
  row.names = FALSE
)

cluster_labels <- setNames(
  paste0("Myeloid_C", cluster_scores$seurat_clusters),
  cluster_scores$seurat_clusters
)
seurat$myeloid_state <- unname(cluster_labels[as.character(seurat$seurat_clusters)])

cluster_proportions <- seurat@meta.data %>%
  count(age_months, myeloid_state) %>%
  group_by(age_months) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()
write.csv(
  cluster_proportions,
  file.path(out, "myeloid_cluster_proportions_by_age.csv"),
  row.names = FALSE
)

markers <- FindAllMarkers(
  seurat,
  only.pos = TRUE,
  min.pct = 0.1,
  logfc.threshold = 0.25,
  verbose = FALSE
)
write.csv(
  markers,
  file.path(out, "myeloid_cluster_markers.csv"),
  row.names = FALSE
)

module_summary <- seurat@meta.data %>%
  group_by(myeloid_state, age_months) %>%
  summarise(
    cells = n(),
    across(all_of(module_columns), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
write.csv(
  module_summary,
  file.path(out, "myeloid_module_scores_by_state_age.csv"),
  row.names = FALSE
)

write.csv(
  seurat@meta.data,
  file.path(out, "myeloid_cell_metadata.csv"),
  row.names = TRUE
)

embeddings <- as.data.frame(Embeddings(seurat, "umap"))
embeddings$myeloid_state <- seurat$myeloid_state
embeddings$age_months <- seurat$age_months
write.csv(
  embeddings,
  file.path(out, "myeloid_umap_coordinates.csv"),
  row.names = TRUE
)

set.seed(2026)
n_bootstrap <- 500
bootstrap_rows <- list()
for (bootstrap_index in seq_len(n_bootstrap)) {
  sampled_cells <- unlist(lapply(c(3, 12, 24), function(age) {
    cells <- rownames(seurat@meta.data)[seurat$age_months == age]
    sample(cells, min(25, length(cells)), replace = TRUE)
  }))
  sampled <- seurat@meta.data[sampled_cells, , drop = FALSE]
  for (module in module_columns) {
    means <- tapply(sampled[[module]], sampled$age_months, mean, na.rm = TRUE)
    bootstrap_rows[[length(bootstrap_rows) + 1]] <- data.frame(
      bootstrap = bootstrap_index,
      module = sub("_score1$", "", module),
      mean_3m = means[["3"]],
      mean_12m = means[["12"]],
      mean_24m = means[["24"]],
      delta_24m_minus_3m = means[["24"]] - means[["3"]]
    )
  }
}
bootstrap_frame <- do.call(rbind, bootstrap_rows)
write.csv(
  bootstrap_frame,
  file.path(out, "myeloid_module_bootstrap_samples.csv"),
  row.names = FALSE
)

saveRDS(seurat, file.path(out, "GSE274279_myeloid_subclustered.rds"))

cat("Myeloid subclustering complete:", out, "\n")
cat("Cells:", ncol(seurat), "Clusters:", length(unique(seurat$seurat_clusters)), "\n")
