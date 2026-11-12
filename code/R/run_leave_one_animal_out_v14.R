#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run with Rscript so repository paths can be resolved.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
source(file.path(dirname(script_file), "common_v14.R"))
configure_v14_library()
suppressPackageStartupMessages({
  library(jsonlite)
  library(limma)
})

REPO_ROOT <- v14_repo_root()
WORK <- v14_prepared_dir(REPO_ROOT)
PRIMARY_RESULTS <- v14_primary_results_dir(REPO_ROOT)
RESULTS <- file.path(REPO_ROOT, "results", "leave_one_animal_out")
dir.create(RESULTS, recursive = TRUE, showWarnings = FALSE)

GROUPS <- c("YC", "MA", "MP", "SP")
SEXES <- c("Female", "Male")
FDR_CUTOFF <- 0.05
ABS_LOG2FC_CUTOFF <- 1
CANDIDATE_GENES <- c(
  "Cd68", "Clec7a", "Fcgr3", "Lyz1", "Lyz2", "Mpeg1", "Ptprc", "Tyrobp"
)

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
    stop("Duplicate gene symbols in: ", path)
  }
  matrix <- as.matrix(data[, -1, drop = FALSE])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- symbols
  matrix
}

sample_master <- read.csv(
  file.path(WORK, "sample_master.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)
common_genes <- read.csv(
  file.path(PRIMARY_RESULTS, "analysis_universe.csv"),
  stringsAsFactors = FALSE
)$gene

tissues <- list(
  cochlea = list(
    dataset = "GSE49543",
    expression = file.path(WORK, "expression_GSE49543_log2.csv")
  ),
  inferior_colliculus = list(
    dataset = "GSE49522",
    expression = file.path(WORK, "expression_GSE49522_log2.csv")
  )
)

fit_shared_genes <- function(animal_id) {
  shared <- character()
  for (tissue in names(tissues)) {
    config <- tissues[[tissue]]
    matrix <- read_expression(config$expression)[common_genes, , drop = FALSE]
    metadata <- sample_master[sample_master$dataset == config$dataset, , drop = FALSE]
    metadata <- metadata[
      match(colnames(matrix), metadata$expression_column),
      ,
      drop = FALSE
    ]
    keep <- metadata$animal_id != animal_id
    matrix <- matrix[, keep, drop = FALSE]
    metadata <- metadata[keep, , drop = FALSE]
    metadata$cohort_group <- factor(
      metadata$cohort_group,
      levels = GROUPS
    )
    metadata$sex <- factor(metadata$sex, levels = SEXES)
    design <- model.matrix(~ 0 + cohort_group + sex, data = metadata)
    if (qr(design)$rank != ncol(design)) {
      stop("Rank-deficient design after removing animal ", animal_id)
    }
    fit <- lmFit(matrix, design)
    contrast <- makeContrasts(
      SP_vs_YC = cohort_groupSP - cohort_groupYC,
      levels = design
    )
    fit <- contrasts.fit(fit, contrast)
    fit <- eBayes(fit, trend = TRUE, robust = TRUE)
    table <- topTable(fit, coef = "SP_vs_YC", number = Inf)
    significant <- rownames(table)[
      table$adj.P.Val < FDR_CUTOFF &
        abs(table$logFC) >= ABS_LOG2FC_CUTOFF
    ]
    if (tissue == "cochlea") {
      shared <- significant
    } else {
      shared <- intersect(shared, significant)
    }
  }
  sort(shared)
}

animals <- sort(unique(sample_master$animal_id))
rows <- vector("list", length(animals))
gene_counts <- setNames(integer(length(CANDIDATE_GENES)), CANDIDATE_GENES)

for (index in seq_along(animals)) {
  animal_id <- animals[index]
  shared <- fit_shared_genes(animal_id)
  candidate_hits <- CANDIDATE_GENES[CANDIDATE_GENES %in% shared]
  gene_counts[candidate_hits] <- gene_counts[candidate_hits] + 1L
  rows[[index]] <- data.frame(
    animal_id = animal_id,
    shared_gene_count = length(shared),
    shared_genes = paste(shared, collapse = ";"),
    candidate_gene_hits = paste(candidate_hits, collapse = ";"),
    stringsAsFactors = FALSE
  )
  cat(
    sprintf(
      "Leave-one-animal-out %d/%d: animal %s, shared genes %d\n",
      index,
      length(animals),
      animal_id,
      length(shared)
    )
  )
}

counts <- do.call(rbind, rows)
write.csv(
  counts,
  file.path(RESULTS, "leave_one_animal_out_counts.csv"),
  row.names = FALSE
)

frequency <- data.frame(
  gene = CANDIDATE_GENES,
  n_models = unname(gene_counts),
  frequency = unname(gene_counts) / nrow(counts),
  stringsAsFactors = FALSE
)
write.csv(
  frequency,
  file.path(RESULTS, "leave_one_animal_out_gene_frequency.csv"),
  row.names = FALSE
)

count_distribution <- as.data.frame(
  table(counts$shared_gene_count),
  stringsAsFactors = FALSE
)
names(count_distribution) <- c("shared_gene_count", "n_models")
count_distribution$shared_gene_count <- as.integer(
  count_distribution$shared_gene_count
)
count_distribution$n_models <- as.integer(count_distribution$n_models)

summary <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  method = "For each of the 41 cochlear animal IDs, we removed that animal from the cochlear model and, when a paired inferior-colliculus sample existed, from the inferior-colliculus model. Both tissue models were then refitted.",
  n_animals = length(animals),
  shared_gene_count = list(
    min = min(counts$shared_gene_count),
    median = median(counts$shared_gene_count),
    max = max(counts$shared_gene_count),
    distribution = count_distribution
  ),
  candidate_gene_frequency = frequency,
  fdr_cutoff = FDR_CUTOFF,
  abs_log2fc_cutoff = ABS_LOG2FC_CUTOFF,
  model_formula = "~ 0 + cohort_group + sex",
  session_info = capture.output(sessionInfo())
)
write_json(
  summary,
  file.path(RESULTS, "leave_one_animal_out_summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("Leave-one-animal-out analysis complete.\n")
print(count_distribution, row.names = FALSE)
