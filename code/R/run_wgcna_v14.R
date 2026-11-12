#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run with Rscript so repository paths can be resolved.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
source(file.path(dirname(script_file), "common_v14.R"))
configure_v14_library()
suppressPackageStartupMessages({
  library(digest)
  library(jsonlite)
  library(limma)
  library(WGCNA)
})
REPO_ROOT <- v14_repo_root()
ANALYSIS_ROOT <- v14_prepared_dir(REPO_ROOT)
SAMPLE_MASTER <- file.path(ANALYSIS_ROOT, "sample_master.csv")
if (!file.exists(SAMPLE_MASTER)) {
  stop("Prepared sample metadata not found: ", SAMPLE_MASTER)
}
RESULTS <- v14_primary_results_dir(REPO_ROOT)
OUT <- file.path(RESULTS, "wgcna")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

GROUPS <- c("YC", "MA", "MP", "SP")
SEXES <- c("Female", "Male")
CONTRASTS <- data.frame(
  contrast_id = c("MA_vs_YC", "SP_vs_MP", "SP_vs_YC"),
  numerator = c("MA", "SP", "SP"),
  denominator = c("YC", "MP", "YC"),
  stringsAsFactors = FALSE
)

TISSUES <- list(
  cochlea = list(
    dataset = "GSE49543",
    expression = file.path(ANALYSIS_ROOT, "expression_GSE49543_log2.csv")
  ),
  inferior_colliculus = list(
    dataset = "GSE49522",
    expression = file.path(ANALYSIS_ROOT, "expression_GSE49522_log2.csv")
  )
)

clean_symbols <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | x == "---"] <- NA_character_
  x
}

read_expression <- function(path, genes) {
  data <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  symbols <- clean_symbols(data[[1]])
  keep <- !is.na(symbols) & symbols %in% genes
  data <- data[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (anyDuplicated(symbols)) {
    stop("Duplicate selected gene symbols in: ", path)
  }
  matrix <- as.matrix(data[, -1, drop = FALSE])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- symbols
  if (any(!is.finite(matrix))) {
    stop("Non-finite expression value found in: ", path)
  }
  matrix
}

load_metadata <- function(matrix, dataset) {
  metadata <- read.csv(
    SAMPLE_MASTER,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "")
  )
  metadata <- metadata[metadata$dataset == dataset, , drop = FALSE]
  index <- match(colnames(matrix), metadata$expression_column)
  if (anyNA(index)) {
    stop("Expression columns missing from sample_master.csv for ", dataset)
  }
  metadata <- metadata[index, , drop = FALSE]
  if (anyNA(metadata[, c("sample_id", "cohort_group", "sex")])) {
    stop("Required metadata contain NA for ", dataset)
  }
  metadata$cohort_group <- factor(metadata$cohort_group, levels = GROUPS)
  metadata$sex <- factor(metadata$sex, levels = SEXES)
  if (anyNA(metadata$cohort_group) || anyNA(metadata$sex)) {
    stop("Unexpected cohort_group or sex level for ", dataset)
  }
  metadata
}

design_for <- function(metadata, dataset) {
  design <- model.matrix(~ 0 + cohort_group + sex, data = metadata)
  rank <- qr(design)$rank
  if (rank != ncol(design)) {
    stop("Rank-deficient WGCNA design for ", dataset)
  }
  design
}

contrast_matrix_for <- function(design) {
  makeContrasts(
    MA_vs_YC = cohort_groupMA - cohort_groupYC,
    SP_vs_MP = cohort_groupSP - cohort_groupMP,
    SP_vs_YC = cohort_groupSP - cohort_groupYC,
    levels = design
  )
}

universe <- read.csv(
  file.path(RESULTS, "analysis_universe.csv"),
  stringsAsFactors = FALSE
)
if (anyDuplicated(universe$gene)) {
  stop("Duplicate genes in analysis_universe.csv")
}

soft_threshold_all <- list()
module_assignments <- list()
module_sizes <- list()
module_eigengenes_all <- list()
module_associations <- list()
wgcna_summary <- list()

set.seed(20261028)
disableWGCNAThreads()

for (tissue in names(TISSUES)) {
  config <- TISSUES[[tissue]]
  variance_column <- paste0(tissue, "_variance")
  if (!variance_column %in% colnames(universe)) {
    stop("Missing variance column in analysis universe: ", variance_column)
  }
  ranked <- universe[order(-universe[[variance_column]]), , drop = FALSE]
  selected_genes <- head(ranked$gene, 5000)
  if (length(selected_genes) < 5000) {
    stop("Fewer than 5,000 valid genes are available for ", tissue)
  }
  expression <- read_expression(config$expression, selected_genes)
  expression <- expression[selected_genes, , drop = FALSE]
  metadata <- load_metadata(expression, config$dataset)
  data_expression <- t(expression)

  gsg <- goodSamplesGenes(data_expression, verbose = 3)
  if (!gsg$allOK) {
    stop("WGCNA input contains samples or genes failing goodSamplesGenes for ", tissue)
  }

  soft <- pickSoftThreshold(
    data_expression,
    powerVector = 1:20,
    networkType = "signed hybrid",
    corFnc = "cor",
    corOptions = list(use = "p"),
    verbose = 5
  )
  soft_table <- as.data.frame(soft$fitIndices)
  soft_table$selected <- FALSE
  eligible <- which(
    is.finite(soft_table$SFT.R.sq) & soft_table$SFT.R.sq >= 0.80
  )
  if (length(eligible) == 0) {
    write.csv(
      data.frame(tissue = tissue, soft_table),
      file.path(OUT, paste0(tissue, "_soft_thresholds.csv")),
      row.names = FALSE
    )
    stop("No prespecified soft-threshold power reached scale-free fit >= 0.80 for ", tissue)
  }
  selected_index <- min(eligible)
  selected_power <- soft_table$Power[selected_index]
  soft_table$selected[selected_index] <- TRUE
  soft_threshold_all[[tissue]] <- data.frame(tissue = tissue, soft_table)

  network <- blockwiseModules(
    data_expression,
    power = selected_power,
    networkType = "signed hybrid",
    TOMType = "signed",
    minModuleSize = 30,
    deepSplit = 2,
    mergeCutHeight = 0.25,
    maxBlockSize = 5000,
    reassignThreshold = 0,
    pamRespectsDendro = FALSE,
    numericLabels = TRUE,
    saveTOMs = FALSE,
    randomSeed = 20261028,
    verbose = 3
  )

  numeric_modules <- network$colors
  module_colors <- labels2colors(numeric_modules)
  module_assignments[[tissue]] <- data.frame(
    tissue = tissue,
    dataset = config$dataset,
    gene = colnames(data_expression),
    module_label = numeric_modules,
    module_color = module_colors,
    is_grey_unassigned = numeric_modules == 0,
    stringsAsFactors = FALSE
  )

  module_table <- as.data.frame(table(module_label = numeric_modules))
  color_by_label <- tapply(
    module_colors,
    numeric_modules,
    function(x) unique(x)[1]
  )
  module_table$module_color <- unname(color_by_label[as.character(module_table$module_label)])
  module_table$tissue <- tissue
  module_sizes[[tissue]] <- module_table[
    , c("tissue", "module_label", "module_color", "Freq")
  ]
  names(module_sizes[[tissue]])[names(module_sizes[[tissue]]) == "Freq"] <- "n_genes"

  module_eigengenes <- network$MEs
  numeric_names <- as.integer(sub("^ME", "", colnames(module_eigengenes)))
  color_names <- labels2colors(numeric_names)
  colnames(module_eigengenes) <- paste0("ME", color_names)
  module_eigengenes_all[[tissue]] <- data.frame(
    tissue = tissue,
    sample_id = metadata$sample_id,
    expression_column = metadata$expression_column,
    cohort_group = as.character(metadata$cohort_group),
    sex = as.character(metadata$sex),
    module_eigengenes,
    check.names = FALSE
  )

  design <- design_for(metadata, config$dataset)
  contrasts <- contrast_matrix_for(design)
  me_matrix <- t(as.matrix(module_eigengenes))
  fit <- lmFit(me_matrix, design)
  fit <- contrasts.fit(fit, contrasts)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  for (contrast_id in CONTRASTS$contrast_id) {
    table <- topTable(
      fit,
      coef = contrast_id,
      number = Inf,
      adjust.method = "BH",
      sort.by = "P"
    )
    table$module <- rownames(table)
    table <- table[, c("module", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
    module_associations[[paste(tissue, contrast_id, sep = "__")]] <- data.frame(
      tissue = tissue,
      contrast = contrast_id,
      table,
      row.names = NULL
    )
  }

  wgcna_summary[[tissue]] <- data.frame(
    dataset = config$dataset,
    tissue = tissue,
    n_samples = nrow(data_expression),
    n_genes = ncol(data_expression),
    selected_power = selected_power,
    scale_free_fit = soft_table$SFT.R.sq[selected_index],
    mean_connectivity = soft_table$mean.k.[selected_index],
    n_modules_excluding_grey = length(unique(numeric_modules[numeric_modules != 0])),
    n_grey_genes = sum(numeric_modules == 0),
    max_block_size = 5000
  )
}

soft_all <- do.call(rbind, soft_threshold_all)
assignments <- do.call(rbind, module_assignments)
sizes <- do.call(rbind, module_sizes)
module_columns <- sort(unique(unlist(lapply(
  module_eigengenes_all,
  function(x) grep("^ME", colnames(x), value = TRUE)
))))
id_columns <- c(
  "tissue", "sample_id", "expression_column", "cohort_group", "sex"
)
eigengenes <- do.call(rbind, lapply(module_eigengenes_all, function(x) {
  missing_columns <- setdiff(module_columns, colnames(x))
  for (column in missing_columns) {
    x[[column]] <- NA_real_
  }
  x[, c(id_columns, module_columns), drop = FALSE]
}))
associations <- do.call(rbind, module_associations)
summary <- do.call(rbind, wgcna_summary)
module_sizes_by_tissue_color <- sizes$n_genes
names(module_sizes_by_tissue_color) <- paste(
  sizes$tissue,
  sizes$module_color,
  sep = "__"
)
association_key <- paste(
  associations$tissue,
  sub("^ME", "", associations$module),
  sep = "__"
)
associations$module_color <- sub("^ME", "", associations$module)
associations$is_grey_unassigned <- associations$module_color == "grey"
associations$module_type <- ifelse(
  associations$is_grey_unassigned,
  "unassigned_grey",
  "biological_module"
)
associations$n_genes <- unname(module_sizes_by_tissue_color[association_key])
associations <- associations[
  ,
  c(
    "tissue", "contrast", "module", "module_color", "module_type",
    "is_grey_unassigned", "n_genes", "logFC", "AveExpr", "t",
    "P.Value", "adj.P.Val", "B"
  ),
  drop = FALSE
]
association_summary <- do.call(rbind, lapply(
  split(associations, list(associations$tissue, associations$contrast), drop = TRUE),
  function(x) {
    data.frame(
      tissue = x$tissue[1],
      contrast = x$contrast[1],
      n_eigengenes_including_grey = nrow(x),
      n_biological_modules = sum(!x$is_grey_unassigned),
      n_biological_modules_fdr_lt_0_05 = sum(
        !x$is_grey_unassigned & x$adj.P.Val < 0.05
      ),
      n_grey_eigengene_fdr_lt_0_05 = sum(
        x$is_grey_unassigned & x$adj.P.Val < 0.05
      ),
      module_fdr_rule = paste(
        "BH-adjusted P < 0.05; descriptive module-level rule;",
        "grey denotes unassigned genes and is not a biological module"
      )
    )
  }
))
rownames(association_summary) <- NULL
fdr_modules <- associations[associations$adj.P.Val < 0.05, , drop = FALSE]
fdr_modules_non_grey <- fdr_modules[!fdr_modules$is_grey_unassigned, , drop = FALSE]
fdr_modules_grey <- fdr_modules[fdr_modules$is_grey_unassigned, , drop = FALSE]

write.csv(soft_all, file.path(OUT, "soft_thresholds.csv"), row.names = FALSE)
write.csv(assignments, file.path(OUT, "module_assignments.csv"), row.names = FALSE)
write.csv(sizes, file.path(OUT, "module_sizes.csv"), row.names = FALSE)
write.csv(eigengenes, file.path(OUT, "module_eigengenes.csv"), row.names = FALSE)
write.csv(
  associations,
  file.path(OUT, "module_contrast_associations.csv"),
  row.names = FALSE
)
write.csv(
  association_summary,
  file.path(OUT, "module_contrast_summary.csv"),
  row.names = FALSE
)
write.csv(
  fdr_modules,
  file.path(OUT, "module_contrast_fdr_lt_0_05.csv"),
  row.names = FALSE
)
write.csv(
  fdr_modules_non_grey,
  file.path(OUT, "module_contrast_fdr_lt_0_05_non_grey.csv"),
  row.names = FALSE
)
write.csv(
  fdr_modules_grey,
  file.path(OUT, "module_contrast_fdr_lt_0_05_grey.csv"),
  row.names = FALSE
)
write.csv(summary, file.path(OUT, "wgcna_summary.csv"), row.names = FALSE)

saveRDS(
  list(
    version = "V14 WGCNA",
    parameters = list(
      n_top_variable_genes = 5000,
      network_type = "signed hybrid",
      tom_type = "signed",
      power_rule = "first power in 1:20 with scale-free fit >= 0.80",
      min_module_size = 30,
      deep_split = 2,
      merge_cut_height = 0.25,
      max_block_size = 5000
    ),
    session_info = sessionInfo()
  ),
  file.path(OUT, "wgcna_run.rds")
)

input_files <- c(
  file.path(ANALYSIS_ROOT, "expression_GSE49543_log2.csv"),
  file.path(ANALYSIS_ROOT, "expression_GSE49522_log2.csv"),
  SAMPLE_MASTER,
  file.path(RESULTS, "analysis_universe.csv")
)
input_checksums <- data.frame(
  path = vapply(
    unname(input_files),
    v14_relative_path,
    character(1),
    root = REPO_ROOT
  ),
  sha256 = vapply(
    unname(input_files),
    digest,
    character(1),
    algo = "sha256",
    file = TRUE
  ),
  stringsAsFactors = FALSE,
  row.names = NULL
)
write.csv(
  input_checksums,
  file.path(OUT, "input_checksums.csv"),
  row.names = FALSE
)

metadata <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  wgcna_version = as.character(packageVersion("WGCNA")),
  model_formula = "~ 0 + cohort_group + sex",
  contrasts = CONTRASTS,
  parameters = list(
    gene_selection = "top 5,000 by variance within each tissue from the common valid analysis universe",
    network_type = "signed hybrid",
    tom_type = "signed",
    power_selection = "lowest power in 1:20 with scale-free topology fit >= 0.80",
    min_module_size = 30,
    deep_split = 2,
    merge_cut_height = 0.25,
    max_block_size = 5000
  ),
  input_checksums = input_checksums,
  outputs = list(
    soft_thresholds = "results/wgcna/soft_thresholds.csv",
    module_assignments = "results/wgcna/module_assignments.csv",
    module_sizes = "results/wgcna/module_sizes.csv",
    module_eigengenes = "results/wgcna/module_eigengenes.csv",
    module_contrast_associations = "results/wgcna/module_contrast_associations.csv",
    module_contrast_summary = "results/wgcna/module_contrast_summary.csv",
    module_contrast_fdr_lt_0_05 = "results/wgcna/module_contrast_fdr_lt_0_05.csv",
    module_contrast_fdr_lt_0_05_non_grey = "results/wgcna/module_contrast_fdr_lt_0_05_non_grey.csv",
    module_contrast_fdr_lt_0_05_grey = "results/wgcna/module_contrast_fdr_lt_0_05_grey.csv",
    summary = "results/wgcna/wgcna_summary.csv"
  ),
  session_info = capture.output(sessionInfo())
)
write_json(
  metadata,
  file.path(OUT, "wgcna_run_metadata.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("V14 WGCNA analysis complete.\n")
print(summary, row.names = FALSE)
