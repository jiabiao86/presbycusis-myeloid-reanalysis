#!/usr/bin/env Rscript

# Descriptive cross-strain single-nucleus context for all four GSE274279
# cochlear matrices. Age and 10x library remain explicitly confounded.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run with Rscript so repository paths can be resolved.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
source(file.path(dirname(script_file), "common_v14.R"))
configure_v14_library()
repo_root <- v14_repo_root()

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(dplyr)
  library(jsonlite)
})

input_dir <- v14_path(
  "V14_SINGLE_CELL_DIR",
  file.path(repo_root, "inputs", "single_cell", "GSE274279_10x")
)
output_dir <- v14_path(
  "V14_SINGLE_CELL_OUTPUT_DIR",
  file.path(v14_derived_dir(repo_root), "GSE274279_four_matrix")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)

samples <- data.frame(
  gsm = c("GSM8446833", "GSM8446834", "GSM8446835", "GSM8446836"),
  library = c("3M", "12M", "24M", "24M_Cochlea_BU"),
  age_months = c(3, 12, 24, 24),
  stringsAsFactors = FALSE
)

candidates <- c(
  "Cd74", "H2-Aa", "H2-Eb1", "Setd1a", "Nfatc2", "Ctss", "Fcgr3",
  "Cd68", "Tyrobp", "Clec4d", "Mpeg1", "Lgals3", "C1qa", "C1qb",
  "C1qc", "C3ar1", "Ms4a7", "Apoe", "Mif", "Ddx1", "Lamp3", "Tbk1",
  "Fgf16", "Csf1r", "P2ry12", "Trem2", "Clec7a", "Lyz1", "Lyz2",
  "Ptprc"
)

markers <- list(
  "Hair cells" = c("Pou4f3", "Myo7a", "Atoh1", "Pcp4", "Gfi1"),
  "Inner hair cells" = c("Slc17a8", "Ocm", "Calb2"),
  "Outer hair cells" = c("Slc26a5", "Kcnq4", "Cdh23"),
  "Supporting cells" = c("Gjb2", "Sox2", "Fbxo2", "Slc1a3", "S100a6"),
  "Spiral ganglion neurons" = c("Nefh", "Snap25", "Tubb3", "Nefl"),
  "Glia/Schwann" = c("Mog", "Mpz", "Plp1", "Sox10"),
  "Stria vascularis" = c("Kcnq1", "Slc12a2", "Atp1a1", "Dct", "Ednrb"),
  "Fibrocytes" = c("Slc26a4", "Col1a1", "Col1a2", "Dcn"),
  "Endothelial" = c("Pecam1", "Cldn5", "Vwf"),
  "Macrophages/Microglia" = c(
    "Ptprc", "Cd68", "Cd74", "C1qa", "C1qb", "Tyrobp", "Laptm5",
    "Fcgr3", "Ctss", "Ms4a7", "Mpeg1", "Csf1r", "P2ry12", "Trem2"
  ),
  "T cells" = c("Cd3d", "Cd3e", "Cd3g", "Lck"),
  "B cells" = c("Cd79a", "Cd79b", "Ms4a1")
)

find_source_file <- function(gsm, suffix) {
  matches <- list.files(
    input_dir,
    pattern = paste0("^", gsm, ".*", suffix, "$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(matches) != 1L) {
    stop(
      "Expected one ", suffix, " file for ", gsm,
      "; found ", length(matches), " under ", input_dir
    )
  }
  normalizePath(matches, mustWork = TRUE)
}

collapse_duplicate_symbols <- function(counts, symbols) {
  symbols <- trimws(as.character(symbols))
  valid <- !is.na(symbols) & nzchar(symbols) & symbols != "---"
  counts <- counts[valid, , drop = FALSE]
  symbols <- symbols[valid]
  keep <- !duplicated(symbols)
  counts <- counts[keep, , drop = FALSE]
  rownames(counts) <- symbols[keep]
  counts
}

read_sample <- function(gsm, library, age_months) {
  matrix_path <- find_source_file(gsm, "matrix\\.mtx\\.gz$")
  barcode_path <- find_source_file(gsm, "barcodes\\.tsv\\.gz$")
  feature_path <- find_source_file(gsm, "features\\.tsv\\.gz$")
  features <- read.delim(
    gzfile(feature_path),
    header = FALSE,
    stringsAsFactors = FALSE
  )
  barcodes <- readLines(gzfile(barcode_path))
  counts <- Matrix::readMM(gzfile(matrix_path))
  if (nrow(counts) != nrow(features) || ncol(counts) != length(barcodes)) {
    stop("Matrix dimensions do not match features and barcodes for ", gsm)
  }
  counts <- as(counts, "CsparseMatrix")
  counts <- collapse_duplicate_symbols(counts, features[[2]])
  colnames(counts) <- paste0(library, "_", barcodes)

  total_counts <- Matrix::colSums(counts)
  detected_genes <- Matrix::colSums(counts > 0)
  mito_rows <- grep("^mt-", rownames(counts), ignore.case = TRUE)
  mito_counts <- if (length(mito_rows)) {
    Matrix::colSums(counts[mito_rows, , drop = FALSE])
  } else {
    rep(0, ncol(counts))
  }
  percent_mt <- ifelse(total_counts > 0, mito_counts / total_counts * 100, NA)
  keep <- (
    detected_genes >= 200 &
      detected_genes < 8000 &
      is.finite(percent_mt) &
      percent_mt < 20 &
      total_counts > 0
  )

  qc <- data.frame(
    gsm = gsm,
    library = library,
    age_months = age_months,
    raw_barcodes = length(barcodes),
    retained_nuclei = sum(keep),
    median_genes_retained = median(detected_genes[keep]),
    median_counts_retained = median(total_counts[keep]),
    median_percent_mt_retained = median(percent_mt[keep]),
    stringsAsFactors = FALSE
  )
  metadata <- data.frame(
    gsm = gsm,
    library = library,
    age_months = age_months,
    nCount_RNA = total_counts[keep],
    nFeature_RNA = detected_genes[keep],
    percent_mt = percent_mt[keep],
    row.names = colnames(counts)[keep],
    stringsAsFactors = FALSE
  )
  list(counts = counts[, keep, drop = FALSE], metadata = metadata, qc = qc)
}

sample_data <- lapply(seq_len(nrow(samples)), function(index) {
  read_sample(
    samples$gsm[index],
    samples$library[index],
    samples$age_months[index]
  )
})
names(sample_data) <- samples$library

combined_counts <- do.call(
  cbind,
  lapply(sample_data, function(item) item$counts)
)
combined_metadata <- do.call(
  rbind,
  lapply(sample_data, function(item) item$metadata)
)
qc_summary <- do.call(
  rbind,
  lapply(sample_data, function(item) item$qc)
)
write.csv(
  qc_summary,
  file.path(output_dir, "GSE274279_four_matrix_QC_by_library.csv"),
  row.names = FALSE
)

seurat <- CreateSeuratObject(
  counts = combined_counts,
  meta.data = combined_metadata,
  project = "GSE274279_four_matrix",
  min.cells = 3,
  min.features = 200
)
seurat[["percent.mt"]] <- seurat$percent_mt
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize",
                        scale.factor = 10000, verbose = FALSE)
seurat <- FindVariableFeatures(
  seurat,
  selection.method = "vst",
  nfeatures = 2500,
  verbose = FALSE
)
seurat <- ScaleData(
  seurat,
  features = VariableFeatures(seurat),
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose = FALSE
)
seurat <- RunPCA(
  seurat,
  features = VariableFeatures(seurat),
  npcs = 40,
  seed.use = 2026,
  verbose = FALSE
)
seurat <- FindNeighbors(seurat, dims = 1:30, verbose = FALSE)
seurat <- FindClusters(seurat, resolution = 0.6, random.seed = 2026,
                       verbose = FALSE)
seurat <- RunUMAP(
  seurat,
  dims = 1:30,
  seed.use = 2026,
  verbose = FALSE
)

data_layer <- GetAssayData(seurat, assay = "RNA", layer = "data")
score_rows <- list()
for (cell_type in names(markers)) {
  available <- intersect(markers[[cell_type]], rownames(data_layer))
  if (length(available) == 0) {
    next
  }
  values <- Matrix::colMeans(data_layer[available, , drop = FALSE])
  score_rows[[cell_type]] <- data.frame(
    cell = colnames(data_layer),
    cluster = as.character(seurat$seurat_clusters),
    cell_type = cell_type,
    score = as.numeric(values),
    markers_available = length(available),
    stringsAsFactors = FALSE
  )
}
cell_scores <- do.call(rbind, score_rows)
cluster_marker_scores <- cell_scores %>%
  group_by(cluster, cell_type) %>%
  summarise(
    mean_marker_score = mean(score),
    markers_available = first(markers_available),
    .groups = "drop"
  )
write.csv(
  cluster_marker_scores,
  file.path(output_dir, "GSE274279_four_matrix_cluster_marker_scores.csv"),
  row.names = FALSE
)

cluster_labels <- cluster_marker_scores %>%
  group_by(cluster) %>%
  slice_max(mean_marker_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(cluster, cell_type)
seurat$cell_type <- unname(
  setNames(cluster_labels$cell_type, cluster_labels$cluster)[
    as.character(seurat$seurat_clusters)
  ]
)

cell_type_counts <- seurat@meta.data %>%
  count(gsm, library, age_months, cell_type, name = "nuclei") %>%
  arrange(age_months, library, cell_type)
write.csv(
  cell_type_counts,
  file.path(output_dir, "GSE274279_four_matrix_cell_type_counts.csv"),
  row.names = FALSE
)

expression_values <- function(gene, selected, matrix = data_layer) {
  if (!gene %in% rownames(matrix) || sum(selected) == 0) {
    return(numeric())
  }
  as.numeric(matrix[gene, selected, drop = FALSE])
}

expression_rows <- list()
for (gene in candidates) {
  if (!gene %in% rownames(data_layer)) {
    next
  }
  for (cell_type in sort(unique(seurat$cell_type))) {
    for (age in sort(unique(seurat$age_months))) {
      selected <- seurat$cell_type == cell_type & seurat$age_months == age
      values <- expression_values(gene, selected, data_layer)
      if (length(values) == 0) {
        next
      }
      expression_rows[[length(expression_rows) + 1L]] <- data.frame(
        gene = gene,
        cell_type = cell_type,
        age_months = as.integer(age),
        cells = length(values),
        mean_log_expression = mean(values),
        percent_expressing = mean(values > 0) * 100,
        stringsAsFactors = FALSE
      )
    }
  }
}
candidate_expression <- do.call(rbind, expression_rows)
write.csv(
  candidate_expression,
  file.path(
    output_dir,
    "GSE274279_four_matrix_candidate_expression_by_cell_type.csv"
  ),
  row.names = FALSE
)

primary_intersection <- c("Cd68", "Clec7a", "Fcgr3", "Lyz1", "Lyz2",
                          "Mpeg1", "Ptprc")
myeloid <- subset(seurat, subset = cell_type == "Macrophages/Microglia")
if (ncol(myeloid) < 10) {
  stop("Too few marker-defined myeloid nuclei for descriptive subclustering.")
}
myeloid_data_layer <- GetAssayData(
  myeloid,
  assay = "RNA",
  layer = "data"
)

myeloid_age_summary <- do.call(rbind, lapply(primary_intersection, function(gene) {
  if (!gene %in% rownames(myeloid_data_layer)) {
    return(NULL)
  }
  do.call(rbind, lapply(sort(unique(myeloid$age_months)), function(age) {
    selected <- myeloid$age_months == age
    values <- expression_values(gene, selected, myeloid_data_layer)
    data.frame(
      gene = gene,
      age_months = as.integer(age),
      cells_total = length(values),
      cells_detected = sum(values > 0),
      percent_expressing = mean(values > 0) * 100,
      mean_log_expression = mean(values),
      stringsAsFactors = FALSE
    )
  }))
}))
write.csv(
  myeloid_age_summary,
  file.path(output_dir, "GSE274279_myeloid_age_descriptive_summary.csv"),
  row.names = FALSE
)

myeloid_library_summary <- do.call(
  rbind,
  lapply(primary_intersection, function(gene) {
    if (!gene %in% rownames(myeloid_data_layer)) {
      return(NULL)
    }
    do.call(rbind, lapply(c(
      "3M", "12M", "24M", "24M_Cochlea_BU"
    ), function(library) {
      selected <- myeloid$library == library
      values <- expression_values(gene, selected, myeloid_data_layer)
      data.frame(
        gene = gene,
        library = library,
        age_months = unique(myeloid$age_months[selected])[[1]],
        cells_total = length(values),
        cells_detected = sum(values > 0),
        percent_expressing = mean(values > 0) * 100,
        mean_log_expression = mean(values),
        stringsAsFactors = FALSE
      )
    }))
  })
)
write.csv(
  myeloid_library_summary,
  file.path(
    output_dir,
    "GSE274279_myeloid_library_descriptive_summary.csv"
  ),
  row.names = FALSE
)

myeloid <- NormalizeData(myeloid, verbose = FALSE)
myeloid <- FindVariableFeatures(myeloid, nfeatures = 1000, verbose = FALSE)
myeloid <- ScaleData(
  myeloid,
  vars.to.regress = c("nCount_RNA", "percent.mt"),
  verbose = FALSE
)
myeloid <- RunPCA(myeloid, npcs = 20, seed.use = 2026, verbose = FALSE)
myeloid <- FindNeighbors(myeloid, dims = 1:15, verbose = FALSE)
myeloid <- FindClusters(myeloid, resolution = 0.4, random.seed = 2026,
                        verbose = FALSE)
myeloid <- RunUMAP(myeloid, dims = 1:15, seed.use = 2026, verbose = FALSE)

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
  genes <- intersect(module_sets[[module]], rownames(myeloid))
  if (length(genes) >= 2) {
    myeloid <- AddModuleScore(
      myeloid,
      features = list(genes),
      name = paste0(module, "_score"),
      assay = "RNA"
    )
  }
}
myeloid$myeloid_state <- paste0("Myeloid_C", myeloid$seurat_clusters)
module_columns <- paste0(names(module_sets), "_score1")
available_module_columns <- intersect(module_columns, colnames(myeloid@meta.data))

cluster_scores <- myeloid@meta.data %>%
  group_by(myeloid_state) %>%
  summarise(
    nuclei = n(),
    across(all_of(available_module_columns), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
write.csv(
  cluster_scores,
  file.path(output_dir, "myeloid_four_matrix_cluster_module_scores.csv"),
  row.names = FALSE
)

cluster_proportions <- myeloid@meta.data %>%
  count(age_months, myeloid_state, name = "nuclei") %>%
  group_by(age_months) %>%
  mutate(proportion = nuclei / sum(nuclei)) %>%
  ungroup()
write.csv(
  cluster_proportions,
  file.path(output_dir, "myeloid_four_matrix_cluster_proportions_by_age.csv"),
  row.names = FALSE
)

module_summary <- myeloid@meta.data %>%
  group_by(myeloid_state, age_months) %>%
  summarise(
    nuclei = n(),
    across(all_of(available_module_columns), ~ mean(.x, na.rm = TRUE)),
    .groups = "drop"
  )
write.csv(
  module_summary,
  file.path(output_dir, "myeloid_four_matrix_module_scores_by_state_age.csv"),
  row.names = FALSE
)

write.csv(
  myeloid@meta.data,
  file.path(output_dir, "myeloid_four_matrix_cell_metadata.csv"),
  row.names = TRUE
)

myeloid_umap <- as.data.frame(Embeddings(myeloid, "umap"))
myeloid_umap$myeloid_state <- myeloid$myeloid_state
myeloid_umap$age_months <- myeloid$age_months
myeloid_umap$library <- myeloid$library
write.csv(
  myeloid_umap,
  file.path(output_dir, "myeloid_four_matrix_umap_coordinates.csv"),
  row.names = TRUE
)

bootstrap_rows <- list()
for (bootstrap_index in seq_len(500)) {
  sampled_cells <- unlist(lapply(c(3, 12, 24), function(age) {
    cells <- rownames(myeloid@meta.data)[myeloid$age_months == age]
    sample(cells, min(25, length(cells)), replace = TRUE)
  }))
  sampled <- myeloid@meta.data[sampled_cells, , drop = FALSE]
  for (module in available_module_columns) {
    means <- tapply(sampled[[module]], sampled$age_months, mean, na.rm = TRUE)
    bootstrap_rows[[length(bootstrap_rows) + 1L]] <- data.frame(
      bootstrap = bootstrap_index,
      module = sub("_score1$", "", module),
      mean_3m = means[["3"]],
      mean_12m = means[["12"]],
      mean_24m = means[["24"]],
      delta_24m_minus_3m = means[["24"]] - means[["3"]]
    )
  }
}
write.csv(
  do.call(rbind, bootstrap_rows),
  file.path(output_dir, "myeloid_four_matrix_cell_bootstrap_samples.csv"),
  row.names = FALSE
)

full_umap <- as.data.frame(Embeddings(seurat, "umap"))
full_umap$cell_type <- seurat$cell_type
full_umap$age_months <- seurat$age_months
full_umap$library <- seurat$library
write.csv(
  full_umap,
  file.path(output_dir, "GSE274279_four_matrix_umap_coordinates.csv"),
  row.names = TRUE
)

saveRDS(
  seurat,
  file.path(output_dir, "GSE274279_four_matrix_annotated.rds"),
  compress = "gzip"
)
saveRDS(
  myeloid,
  file.path(output_dir, "GSE274279_four_matrix_myeloid.rds"),
  compress = "gzip"
)

summary <- list(
  source = "Xia et al., Cell Reports 2025;44:115781",
  strain = "C57BL/6J",
  gsm = samples$gsm,
  raw_barcodes = setNames(qc_summary$raw_barcodes, qc_summary$library),
  retained_nuclei = setNames(qc_summary$retained_nuclei, qc_summary$library),
  total_retained_nuclei = ncol(seurat),
  genes = nrow(seurat),
  myeloid_nuclei = ncol(myeloid),
  myeloid_by_age = as.list(table(myeloid$age_months)),
  myeloid_by_library = as.list(table(myeloid$library)),
  cell_types = as.list(table(seurat$cell_type)),
  analytic_boundary = paste(
    "Age, library, and 10x batch are not separable. The animal-level",
    "relationship of the 24M_Cochlea_BU library to the 24M library is not",
    "documented in GEO."
  )
)
write_json(
  summary,
  file.path(output_dir, "GSE274279_four_matrix_summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

cat(
  "Retained", ncol(seurat), "nuclei;",
  "marker-defined myeloid:", ncol(myeloid), "nuclei.\n"
)
print(summary)
