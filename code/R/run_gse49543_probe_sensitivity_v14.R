#!/usr/bin/env Rscript

# Sensitivity of the SP-versus-YC intersection to GSE49543 probe-to-symbol
# aggregation. The inferior-colliculus model and frozen universe are unchanged.

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run with Rscript so repository paths can be resolved.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
source(file.path(dirname(script_file), "common_v14.R"))
configure_v14_library()

suppressPackageStartupMessages({
  library(limma)
  library(jsonlite)
})

repo_root <- v14_repo_root()
input_dir <- v14_path(
  "V14_PROBE_SENSITIVITY_DIR",
  file.path(repo_root, "inputs", "prepared", "probe_sensitivity")
)
output_dir <- v14_path(
  "V14_PROBE_SENSITIVITY_OUTPUT_DIR",
  file.path(
    repo_root, "results", "sensitivity_probe_aggregation"
  )
)
sample_master_path <- v14_path(
  "V14_SAMPLE_MASTER",
  file.path(repo_root, "inputs", "prepared", "sample_master.csv")
)
primary_dir <- v14_path(
  "V14_PRIMARY_RESULTS_DIR",
  file.path(repo_root, "results", "primary")
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

fdr_cutoff <- 0.05
effect_cutoff <- 1
schemes <- c("first_probe", "mean_probe", "max_mean_probe")
primary_genes <- c("Cd68", "Clec7a", "Fcgr3", "Lyz1", "Lyz2", "Mpeg1", "Ptprc")

clean_symbols <- function(x) {
  x <- trimws(as.character(x))
  x[is.na(x) | x == "" | x == "---"] <- NA_character_
  x
}

read_expression <- function(path) {
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
    stop("Duplicate gene symbols in ", path)
  }
  matrix <- as.matrix(data[, -1, drop = FALSE])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- symbols
  matrix
}

sample_master <- read.csv(
  sample_master_path,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
cochlea_metadata <- sample_master[
  sample_master$dataset == "GSE49543" &
    sample_master$include_main %in% c(TRUE, "TRUE"),
  ,
  drop = FALSE
]
if (nrow(cochlea_metadata) != 41L) {
  stop("Expected 41 primary GSE49543 samples.")
}
cochlea_metadata$cohort_group <- factor(
  cochlea_metadata$cohort_group,
  levels = c("YC", "MA", "MP", "SP")
)
cochlea_metadata$sex <- factor(
  cochlea_metadata$sex,
  levels = c("Female", "Male")
)
design <- model.matrix(~ 0 + cohort_group + sex, data = cochlea_metadata)
contrast <- makeContrasts(
  SP_vs_YC = cohort_groupSP - cohort_groupYC,
  levels = design
)

primary_universe_path <- file.path(primary_dir, "analysis_universe.csv")
if (!file.exists(primary_universe_path)) {
  stop("Primary gene universe not found: ", primary_universe_path)
}
frozen_universe <- read.csv(
  primary_universe_path,
  stringsAsFactors = FALSE
)$gene

ic_path <- file.path(
  primary_dir,
  "limma",
  "inferior_colliculus_SP_vs_YC.csv"
)
if (!file.exists(ic_path)) {
  stop("Primary inferior-colliculus result not found: ", ic_path)
}
ic_table <- read.csv(ic_path, stringsAsFactors = FALSE)
ic_significant <- ic_table$gene[
  ic_table$adj.P.Val < fdr_cutoff & abs(ic_table$logFC) >= effect_cutoff
]

full_results <- list()
summary_rows <- list()
membership_rows <- list()
for (scheme in schemes) {
  expression_path <- file.path(
    input_dir,
    paste0("GSE49543_", scheme, "_log2.csv")
  )
  expression <- read_expression(expression_path)
  index <- match(cochlea_metadata$expression_column, colnames(expression))
  if (anyNA(index)) {
    stop(
      "Sample columns missing from ", scheme, ": ",
      paste(cochlea_metadata$expression_column[is.na(index)], collapse = ", ")
    )
  }
  expression <- expression[, index, drop = FALSE]
  if (!identical(colnames(expression), cochlea_metadata$expression_column)) {
    stop("Expression columns and sample metadata are not aligned for ", scheme)
  }
  missing_universe <- setdiff(frozen_universe, rownames(expression))
  if (length(missing_universe) > 0) {
    stop(
      scheme, " lacks ", length(missing_universe),
      " genes in the frozen primary universe."
    )
  }

  fit <- lmFit(expression[frozen_universe, , drop = FALSE], design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit, trend = TRUE, robust = TRUE)
  table <- topTable(
    fit,
    coef = "SP_vs_YC",
    number = Inf,
    adjust.method = "BH",
    sort.by = "P"
  )
  table$gene <- rownames(table)
  table <- table[, c(
    "gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"
  )]
  table$significant <- with(
    table,
    adj.P.Val < fdr_cutoff & abs(logFC) >= effect_cutoff
  )
  write.csv(
    table,
    file.path(output_dir, paste0("cochlea_SP_vs_YC_", scheme, ".csv")),
    row.names = FALSE
  )
  full_results[[scheme]] <- table

  cochlea_significant <- table$gene[table$significant]
  intersection <- intersect(cochlea_significant, ic_significant)
  retained <- intersect(primary_genes, intersection)
  gained <- setdiff(intersection, primary_genes)
  summary_rows[[scheme]] <- data.frame(
    probe_rule = scheme,
    genes_tested = nrow(table),
    cochlea_significant_up = sum(
      table$significant & table$logFC > 0
    ),
    cochlea_significant_down = sum(
      table$significant & table$logFC < 0
    ),
    inferior_colliculus_significant = length(ic_significant),
    intersection_size = length(intersection),
    primary_seven_retained = length(retained),
    additional_intersection_genes = paste(sort(gained), collapse = ";"),
    stringsAsFactors = FALSE
  )
  membership_rows[[scheme]] <- data.frame(
    probe_rule = scheme,
    gene = sort(union(primary_genes, intersection)),
    primary_seven = sort(union(primary_genes, intersection)) %in% primary_genes,
    in_intersection = sort(union(primary_genes, intersection)) %in% intersection,
    cochlea_logFC = table$logFC[
      match(sort(union(primary_genes, intersection)), table$gene)
    ],
    cochlea_adj_p = table$adj.P.Val[
      match(sort(union(primary_genes, intersection)), table$gene)
    ],
    stringsAsFactors = FALSE
  )
}

summary_table <- do.call(rbind, summary_rows)
rownames(summary_table) <- NULL
write.csv(
  summary_table,
  file.path(output_dir, "probe_scheme_summary.csv"),
  row.names = FALSE
)
membership_table <- do.call(rbind, membership_rows)
rownames(membership_table) <- NULL
write.csv(
  membership_table,
  file.path(output_dir, "probe_scheme_intersection_membership.csv"),
  row.names = FALSE
)

pairwise <- list()
for (left in schemes) {
  for (right in schemes) {
    if (left >= right) {
      next
    }
    left_table <- full_results[[left]]
    right_table <- full_results[[right]]
    aligned <- merge(
      left_table[, c("gene", "logFC")],
      right_table[, c("gene", "logFC")],
      by = "gene",
      suffixes = c("_left", "_right")
    )
    pairwise[[paste(left, right, sep = "__")]] <- data.frame(
      scheme_left = left,
      scheme_right = right,
      genes = nrow(aligned),
      spearman_logfc = cor(
        aligned$logFC_left,
        aligned$logFC_right,
        method = "spearman"
      ),
      max_abs_logfc_difference = max(
        abs(aligned$logFC_left - aligned$logFC_right)
      ),
      stringsAsFactors = FALSE
    )
  }
}
pairwise_table <- do.call(rbind, pairwise)
rownames(pairwise_table) <- NULL
write.csv(
  pairwise_table,
  file.path(output_dir, "probe_scheme_pairwise_concordance.csv"),
  row.names = FALSE
)

metadata <- list(
  source_preparation = "code/python/prepare_gse49543_probe_sensitivity.py",
  input_directory = v14_relative_path(input_dir, repo_root),
  sample_master = v14_relative_path(sample_master_path, repo_root),
  inferior_colliculus_result = v14_relative_path(ic_path, repo_root),
  frozen_universe = v14_relative_path(primary_universe_path, repo_root),
  model = "~ 0 + cohort_group + sex",
  contrast = "SP - YC",
  thresholds = list(
    bh_adjusted_p = fdr_cutoff,
    absolute_log2_fold_change = effect_cutoff
  ),
  rules = list(
    first_probe = "First valid source row per exact gene symbol.",
    mean_probe = "Mean linear signal across rows sharing an exact gene symbol.",
    max_mean_probe = paste(
      "Row with the highest mean linear signal across the 41 samples",
      "for each exact gene symbol."
    )
  ),
  interpretation = paste(
    "Probe aggregation changes the cochlear transcript set tested;",
    "the inferior-colliculus model and frozen gene universe are unchanged."
  )
)
write_json(
  metadata,
  file.path(output_dir, "probe_scheme_analysis_metadata.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(output_dir, "sessionInfo.txt")
)

print(summary_table, row.names = FALSE)
print(pairwise_table, row.names = FALSE)
