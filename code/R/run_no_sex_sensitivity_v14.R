#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

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
})

REPO_ROOT <- v14_repo_root()
GEO_REBUILD_DIR <- v14_path(
  "V14_GEO_REBUILD_DIR",
  file.path(REPO_ROOT, "validation", "geo_rebuild")
)
SAMPLE_MASTER <- v14_path(
  "V14_SAMPLE_MASTER",
  file.path(
    REPO_ROOT, "..", "S1_sample_master", "sample_master.csv"
  )
)
UNIVERSE_FILE <- v14_path(
  "V14_ANALYSIS_UNIVERSE",
  file.path(
    REPO_ROOT, "..", "S2_model_and_contrasts", "analysis_universe.csv"
  )
)
SEX_ADJUSTED_DIR <- v14_path(
  "V14_SEX_ADJUSTED_LIMMA_DIR",
  file.path(
    REPO_ROOT, "..", "S9_significant_and_intersection",
    "primary", "limma"
  )
)
RESULTS_DIR <- v14_path(
  "V14_NO_SEX_RESULTS_DIR",
  file.path(
    REPO_ROOT, "..", "S9_significant_and_intersection", "sensitivity_no_sex"
  )
)
dir.create(RESULTS_DIR, recursive = TRUE, showWarnings = FALSE)

GROUPS <- c("YC", "MA", "MP", "SP")
FDR_CUTOFF <- 0.05
ABS_LOG2FC_CUTOFF <- 1
TISSUES <- list(
  cochlea = list(
    dataset = "GSE49543",
    expression = file.path(GEO_REBUILD_DIR, "expression_GSE49543_log2.csv")
  ),
  inferior_colliculus = list(
    dataset = "GSE49522",
    expression = file.path(GEO_REBUILD_DIR, "expression_GSE49522_log2.csv")
  )
)
CONTRASTS <- data.frame(
  contrast_id = c("MA_vs_YC", "SP_vs_MP", "SP_vs_YC"),
  numerator = c("MA", "SP", "SP"),
  denominator = c("YC", "MP", "YC"),
  stringsAsFactors = FALSE
)

clean_symbols <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | x == "---"] <- NA_character_
  x
}

read_expression <- function(path) {
  if (!file.exists(path)) {
    stop("Missing expression input: ", path)
  }
  data <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    fileEncoding = "UTF-8-BOM"
  )
  symbols <- clean_symbols(data[[1]])
  keep <- !is.na(symbols)
  data <- data[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (anyDuplicated(symbols)) {
    stop("Duplicate gene symbols after cleaning: ", path)
  }
  matrix <- as.matrix(data[, -1, drop = FALSE])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- symbols
  if (any(!is.finite(matrix))) {
    stop("Non-finite expression value found: ", path)
  }
  matrix
}

sample_master <- read.csv(
  SAMPLE_MASTER,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)
frozen_genes <- read.csv(
  UNIVERSE_FILE,
  check.names = FALSE,
  stringsAsFactors = FALSE
)$gene
if (length(frozen_genes) != 12437 || anyDuplicated(frozen_genes)) {
  stop("Unexpected frozen analysis universe.")
}

tissue_matrices <- lapply(TISSUES, function(config) {
  matrix <- read_expression(config$expression)
  if (any(!frozen_genes %in% rownames(matrix))) {
    stop("Frozen genes unavailable in ", config$dataset)
  }
  matrix[frozen_genes, , drop = FALSE]
})

load_metadata <- function(matrix, dataset) {
  metadata <- sample_master[sample_master$dataset == dataset, , drop = FALSE]
  index <- match(colnames(matrix), metadata$expression_column)
  if (anyNA(index) || anyDuplicated(metadata$expression_column)) {
    stop("Expression metadata could not be aligned for ", dataset)
  }
  metadata <- metadata[index, , drop = FALSE]
  metadata$cohort_group <- factor(metadata$cohort_group, levels = GROUPS)
  if (anyNA(metadata$cohort_group)) {
    stop("Unexpected cohort level in ", dataset)
  }
  metadata
}

fit_no_sex <- function(matrix, metadata, tissue, dataset) {
  design <- model.matrix(~ 0 + cohort_group, data = metadata)
  rank <- qr(design)$rank
  if (rank != ncol(design)) {
    stop(
      "Rank-deficient no-sex design for ", dataset,
      ": rank=", rank, ", columns=", ncol(design)
    )
  }
  write.csv(
    data.frame(
      sample_id = metadata$sample_id,
      expression_column = metadata$expression_column,
      cohort_group = as.character(metadata$cohort_group),
      design,
      check.names = FALSE
    ),
    file.path(RESULTS_DIR, paste0(tissue, "_no_sex_design_matrix.csv")),
    row.names = FALSE
  )

  contrast_matrix <- makeContrasts(
    MA_vs_YC = cohort_groupMA - cohort_groupYC,
    SP_vs_MP = cohort_groupSP - cohort_groupMP,
    SP_vs_YC = cohort_groupSP - cohort_groupYC,
    levels = design
  )
  fit <- lmFit(matrix, design)
  fit <- contrasts.fit(fit, contrast_matrix)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  lapply(CONTRASTS$contrast_id, function(contrast_id) {
    table <- topTable(
      fit,
      coef = contrast_id,
      number = Inf,
      adjust.method = "BH",
      sort.by = "P"
    )
    table$gene <- rownames(table)
    table[, c(
      "gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"
    )]
  }) |>
    setNames(CONTRASTS$contrast_id)
}

no_sex_tables <- list()
design_summary <- list()
for (tissue in names(TISSUES)) {
  config <- TISSUES[[tissue]]
  metadata <- load_metadata(tissue_matrices[[tissue]], config$dataset)
  no_sex_tables[[tissue]] <- fit_no_sex(
    tissue_matrices[[tissue]],
    metadata,
    tissue,
    config$dataset
  )
  design_summary[[tissue]] <- data.frame(
    tissue = tissue,
    dataset = config$dataset,
    n_samples = nrow(metadata),
    model = "~ 0 + cohort_group",
    design_columns = 4,
    design_rank = 4
  )
}
write.csv(
  do.call(rbind, design_summary),
  file.path(RESULTS_DIR, "no_sex_design_summary.csv"),
  row.names = FALSE
)

count_rows <- list()
shared_rows <- list()
comparison_rows <- list()
for (tissue in names(TISSUES)) {
  for (contrast_id in CONTRASTS$contrast_id) {
    table <- no_sex_tables[[tissue]][[contrast_id]]
    significant <- table$adj.P.Val < FDR_CUTOFF &
      abs(table$logFC) >= ABS_LOG2FC_CUTOFF
    table$significant <- significant
    table$direction <- ifelse(
      significant & table$logFC > 0,
      "up",
      ifelse(significant & table$logFC < 0, "down", "not_significant")
    )
    write.csv(
      table,
      file.path(
        RESULTS_DIR,
        paste0(tissue, "_", contrast_id, "_no_sex.csv")
      ),
      row.names = FALSE
    )
    count_rows[[paste(tissue, contrast_id, sep = "__")]] <- data.frame(
      tissue = tissue,
      dataset = TISSUES[[tissue]]$dataset,
      contrast = contrast_id,
      n_up = sum(table$direction == "up"),
      n_down = sum(table$direction == "down"),
      n_fdr_only = sum(table$adj.P.Val < FDR_CUTOFF),
      n_effect_size_only = sum(abs(table$logFC) >= ABS_LOG2FC_CUTOFF),
      n_significant = sum(significant),
      fdr_cutoff = FDR_CUTOFF,
      abs_log2fc_cutoff = ABS_LOG2FC_CUTOFF
    )
    selected_rows <- which(significant)
    shared_rows[[paste(tissue, contrast_id, sep = "__")]] <- data.frame(
      tissue = rep(tissue, length(selected_rows)),
      contrast = rep(contrast_id, length(selected_rows)),
      gene = table$gene[selected_rows],
      logFC_no_sex = table$logFC[selected_rows],
      adj.P.Val_no_sex = table$adj.P.Val[selected_rows],
      direction_no_sex = table$direction[selected_rows],
      stringsAsFactors = FALSE
    )

    sex_file <- file.path(
      SEX_ADJUSTED_DIR,
      paste0(tissue, "_", contrast_id, ".csv")
    )
    sex_adjusted <- read.csv(
      sex_file,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
    merged <- merge(
      sex_adjusted[, c("gene", "logFC", "adj.P.Val")],
      table[, c("gene", "logFC", "adj.P.Val")],
      by = "gene",
      suffixes = c("_sex_adjusted", "_no_sex"),
      sort = FALSE
    )
    merged$same_direction <- sign(merged$logFC_sex_adjusted) ==
      sign(merged$logFC_no_sex)
    merged$significant_sex_adjusted <- merged$adj.P.Val_sex_adjusted <
      FDR_CUTOFF & abs(merged$logFC_sex_adjusted) >= ABS_LOG2FC_CUTOFF
    merged$significant_no_sex <- merged$adj.P.Val_no_sex < FDR_CUTOFF &
      abs(merged$logFC_no_sex) >= ABS_LOG2FC_CUTOFF
    merged$tissue <- tissue
    merged$contrast <- contrast_id
    comparison_rows[[paste(tissue, contrast_id, sep = "__")]] <- merged
  }
}

counts <- do.call(rbind, count_rows)
shared <- do.call(rbind, shared_rows)
rownames(shared) <- NULL
comparison <- do.call(rbind, comparison_rows)
rownames(comparison) <- NULL
write.csv(
  counts,
  file.path(RESULTS_DIR, "no_sex_significant_counts.csv"),
  row.names = FALSE
)
write.csv(
  shared,
  file.path(RESULTS_DIR, "no_sex_significant_genes.csv"),
  row.names = FALSE
)
write.csv(
  comparison,
  file.path(RESULTS_DIR, "no_sex_vs_sex_adjusted_all_genes.csv"),
  row.names = FALSE
)

tissues <- names(TISSUES)
intersection_rows <- list()
summary_rows <- list()
for (contrast_id in CONTRASTS$contrast_id) {
  sex_sets <- lapply(tissues, function(tissue) {
    data <- read.csv(
      file.path(SEX_ADJUSTED_DIR, paste0(tissue, "_", contrast_id, ".csv")),
      stringsAsFactors = FALSE
    )
    data$gene[
      data$adj.P.Val < FDR_CUTOFF & abs(data$logFC) >= ABS_LOG2FC_CUTOFF
    ]
  })
  no_sex_sets <- lapply(tissues, function(tissue) {
    data <- no_sex_tables[[tissue]][[contrast_id]]
    data$gene[
      data$adj.P.Val < FDR_CUTOFF & abs(data$logFC) >= ABS_LOG2FC_CUTOFF
    ]
  })
  names(sex_sets) <- tissues
  names(no_sex_sets) <- tissues
  sex_intersection <- Reduce(intersect, sex_sets)
  no_sex_intersection <- Reduce(intersect, no_sex_sets)
  intersection_rows[[contrast_id]] <- data.frame(
    contrast = rep(
      contrast_id,
      length(sex_intersection) + length(no_sex_intersection)
    ),
    model = c(
      rep("sex_adjusted", length(sex_intersection)),
      rep("no_sex", length(no_sex_intersection))
    ),
    gene = c(sex_intersection, no_sex_intersection),
    stringsAsFactors = FALSE
  )
  summary_rows[[contrast_id]] <- data.frame(
    contrast = contrast_id,
    sex_adjusted_intersection_n = length(sex_intersection),
    no_sex_intersection_n = length(no_sex_intersection),
    retained_in_no_sex = sum(sex_intersection %in% no_sex_intersection),
    gained_in_no_sex = sum(!no_sex_intersection %in% sex_intersection),
    lost_in_no_sex = sum(!sex_intersection %in% no_sex_intersection)
  )
}
intersections <- do.call(rbind, intersection_rows)
summary <- do.call(rbind, summary_rows)
rownames(intersections) <- NULL
rownames(summary) <- NULL
write.csv(
  intersections,
  file.path(RESULTS_DIR, "no_sex_intersection_membership.csv"),
  row.names = FALSE
)
write.csv(
  summary,
  file.path(RESULTS_DIR, "no_sex_intersection_summary.csv"),
  row.names = FALSE
)

sp_yc_comparison <- comparison[comparison$contrast == "SP_vs_YC", ]
sp_yc_primary_genes <- sort(unique(c(
  read.csv(
    file.path(
      REPO_ROOT, "..", "S9_significant_and_intersection",
      "primary", "shared_genes.csv"
    ),
    stringsAsFactors = FALSE
  )$gene,
  intersections$gene[
    intersections$contrast == "SP_vs_YC" &
      intersections$model == "no_sex"
  ]
)))
primary_intersection_table <- sp_yc_comparison[
  sp_yc_comparison$gene %in% sp_yc_primary_genes,
  c(
    "tissue", "gene", "logFC_sex_adjusted", "logFC_no_sex",
    "adj.P.Val_sex_adjusted", "adj.P.Val_no_sex", "same_direction",
    "significant_sex_adjusted", "significant_no_sex"
  )
]
write.csv(
  primary_intersection_table,
  file.path(RESULTS_DIR, "no_sex_primary_intersection_comparison.csv"),
  row.names = FALSE
)

correlation_rows <- lapply(tissues, function(tissue) {
  data <- sp_yc_comparison[sp_yc_comparison$tissue == tissue, ]
  data.frame(
    tissue = tissue,
    n_genes = nrow(data),
    spearman_logFC = cor(
      data$logFC_sex_adjusted,
      data$logFC_no_sex,
      method = "spearman"
    ),
    sign_concordance = mean(data$same_direction),
    max_abs_logFC_difference = max(
      abs(data$logFC_sex_adjusted - data$logFC_no_sex)
    )
  )
})
correlations <- do.call(rbind, correlation_rows)
write.csv(
  correlations,
  file.path(RESULTS_DIR, "no_sex_sp_vs_yc_effect_concordance.csv"),
  row.names = FALSE
)

input_files <- c(
  TISSUES$cochlea$expression,
  TISSUES$inferior_colliculus$expression,
  SAMPLE_MASTER,
  UNIVERSE_FILE
)
checksums <- data.frame(
  path = vapply(
    input_files,
    v14_relative_path,
    character(1),
    root = REPO_ROOT
  ),
  sha256 = vapply(
    input_files,
    digest,
    character(1),
    algo = "sha256",
    file = TRUE
  ),
  stringsAsFactors = FALSE
)
write.csv(
  checksums,
  file.path(RESULTS_DIR, "no_sex_input_checksums.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(RESULTS_DIR, "sessionInfo.txt")
)
write_json(
  list(
    generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
    model = "~ 0 + cohort_group",
    rationale = paste(
      "Sensitivity model omitting sex adjustment because the inferior",
      "colliculus SP group contained only one male sample."
    ),
    frozen_gene_universe = length(frozen_genes),
    fdr_method = "Benjamini-Hochberg within tissue and contrast",
    fdr_cutoff = FDR_CUTOFF,
    abs_log2fc_cutoff = ABS_LOG2FC_CUTOFF,
    limma = list(
      lmFit = "limma::lmFit",
      eBayes = "limma::eBayes(fit, trend = TRUE, robust = TRUE)"
    ),
    input_checksums = checksums,
    session_info = capture.output(sessionInfo())
  ),
  file.path(RESULTS_DIR, "no_sex_sensitivity_metadata.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("No-sex sensitivity outputs written to", RESULTS_DIR, "\n")
print(summary)
print(correlations)
