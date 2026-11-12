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
})
REPO_ROOT <- v14_repo_root()
PREPARED_DIR <- v14_prepared_dir(REPO_ROOT)
PRIMARY_RESULTS <- v14_primary_results_dir(REPO_ROOT)
SAMPLE_MASTER <- file.path(PREPARED_DIR, "sample_master.csv")
if (!file.exists(SAMPLE_MASTER)) {
  stop("Prepared sample metadata not found: ", SAMPLE_MASTER)
}
ANALYSIS_ROOT <- PREPARED_DIR
RUN_LABEL <- Sys.getenv("V14_RUN_LABEL", unset = "primary")
EXCLUDE_COCHLEA_GSMS <- trimws(
  Sys.getenv("V14_EXCLUDE_COCHLEA_GSM", unset = "")
)
FREEZE_GENE_UNIVERSE <- identical(
  tolower(Sys.getenv("V14_FREEZE_GENE_UNIVERSE", unset = "false")),
  "true"
)
DEFAULT_RESULTS <- if (RUN_LABEL == "primary") {
  PRIMARY_RESULTS
} else {
  file.path(REPO_ROOT, "results", RUN_LABEL)
}
RESULTS <- normalizePath(
  Sys.getenv("V14_RESULTS_DIR", unset = DEFAULT_RESULTS),
  mustWork = FALSE
)
LIMMA_DIR <- file.path(RESULTS, "limma")
DESIGN_DIR <- file.path(RESULTS, "design")
dir.create(LIMMA_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(DESIGN_DIR, recursive = TRUE, showWarnings = FALSE)

if (RUN_LABEL == "primary" && nzchar(EXCLUDE_COCHLEA_GSMS)) {
  stop("A sample exclusion requires a non-primary V14_RUN_LABEL.")
}
if (RUN_LABEL != "primary" && !nzchar(EXCLUDE_COCHLEA_GSMS)) {
  stop("A non-primary sensitivity run requires V14_EXCLUDE_COCHLEA_GSM.")
}

GROUPS <- c("YC", "MA", "MP", "SP")
SEXES <- c("Female", "Male")
FDR_CUTOFF <- 0.05
ABS_LOG2FC_CUTOFF <- 1

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

CONTRASTS <- data.frame(
  contrast_id = c("MA_vs_YC", "SP_vs_MP", "SP_vs_YC"),
  numerator = c("MA", "SP", "SP"),
  denominator = c("YC", "MP", "YC"),
  interpretation = c(
    "Middle-aged normal-hearing cohort minus young normal-hearing cohort.",
    "SP cohort minus MP cohort according to the source labels; age remains partly confounded.",
    "SP cohort minus YC cohort according to the source labels; descriptive combined age and hearing-status difference."
  ),
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
  if (ncol(data) < 2) {
    stop("Expression input has fewer than two columns: ", path)
  }
  symbols <- clean_symbols(data[[1]])
  keep <- !is.na(symbols)
  data <- data[keep, , drop = FALSE]
  symbols <- symbols[keep]
  if (anyDuplicated(symbols)) {
    stop("Duplicate gene symbols after dropping invalid identifiers in: ", path)
  }
  matrix <- as.matrix(data[, -1, drop = FALSE])
  storage.mode(matrix) <- "double"
  rownames(matrix) <- symbols
  if (anyDuplicated(colnames(matrix))) {
    stop("Duplicate expression column names in: ", path)
  }
  if (any(!is.finite(matrix))) {
    stop("Non-finite expression value found in: ", path)
  }
  matrix
}

load_metadata <- function(matrix, dataset) {
  path <- SAMPLE_MASTER
  metadata <- read.csv(
    path,
    check.names = FALSE,
    stringsAsFactors = FALSE,
    na.strings = c("NA", "")
  )
  metadata <- metadata[metadata$dataset == dataset, , drop = FALSE]
  if (anyDuplicated(metadata$expression_column)) {
    stop("Duplicate expression columns in sample_master.csv for ", dataset)
  }
  index <- match(colnames(matrix), metadata$expression_column)
  if (anyNA(index)) {
    stop("Expression columns missing from sample_master.csv for ", dataset)
  }
  metadata <- metadata[index, , drop = FALSE]
  if (!identical(metadata$expression_column, colnames(matrix))) {
    stop("Expression columns and metadata are not aligned for ", dataset)
  }
  required <- c("sample_id", "cohort_group", "sex")
  if (anyNA(metadata[, required])) {
    stop("Required metadata contain NA for ", dataset)
  }
  metadata$cohort_group <- factor(metadata$cohort_group, levels = GROUPS)
  metadata$sex <- factor(metadata$sex, levels = SEXES)
  if (anyNA(metadata$cohort_group) || anyNA(metadata$sex)) {
    stop("Unexpected cohort_group or sex level for ", dataset)
  }
  metadata
}

build_design <- function(metadata, tissue, dataset) {
  if (any(!complete.cases(metadata[, c("cohort_group", "sex")]))) {
    stop("Incomplete cohort_group or sex in ", dataset)
  }
  group_sex_counts <- as.data.frame(
    table(
      cohort_group = metadata$cohort_group,
      sex = metadata$sex
    )
  )
  write.csv(
    group_sex_counts,
    file.path(DESIGN_DIR, paste0(tissue, "_group_sex_counts.csv")),
    row.names = FALSE
  )
  design <- model.matrix(~ 0 + cohort_group + sex, data = metadata)
  rank <- qr(design)$rank
  if (rank != ncol(design)) {
    stop(
      "Rank-deficient design for ", dataset, ": rank=", rank,
      ", columns=", ncol(design)
    )
  }
  write.csv(
    data.frame(
      sample_id = metadata$sample_id,
      expression_column = metadata$expression_column,
      cohort_group = as.character(metadata$cohort_group),
      sex = as.character(metadata$sex),
      design,
      check.names = FALSE
    ),
    file.path(DESIGN_DIR, paste0(tissue, "_design_matrix.csv")),
    row.names = FALSE
  )
  writeLines(
    c(
      paste0("dataset: ", dataset),
      paste0("tissue: ", tissue),
      paste0("model: ~ 0 + cohort_group + sex"),
      paste0("n_samples: ", nrow(design)),
      paste0("n_design_columns: ", ncol(design)),
      paste0("design_rank: ", rank),
      paste0("full_rank: ", rank == ncol(design))
    ),
    file.path(DESIGN_DIR, paste0(tissue, "_design_rank.txt"))
  )
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

tissue_matrices <- lapply(TISSUES, function(x) read_expression(x$expression))

sample_master <- read.csv(
  SAMPLE_MASTER,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "")
)

exclusion_record <- data.frame(
  run_label = character(),
  dataset = character(),
  tissue = character(),
  gsm = character(),
  animal_id = character(),
  adjudicated_group = character(),
  expression_column = character(),
  reason = character(),
  stringsAsFactors = FALSE
)
if (nzchar(EXCLUDE_COCHLEA_GSMS)) {
  exclude_gsms <- trimws(unlist(strsplit(
    EXCLUDE_COCHLEA_GSMS,
    split = ",",
    fixed = TRUE
  )))
  exclude_gsms <- unique(exclude_gsms[nzchar(exclude_gsms)])
  selected <- sample_master[
    sample_master$dataset == "GSE49543" &
      sample_master$gsm %in% exclude_gsms,
    ,
    drop = FALSE
  ]
  if (nrow(selected) != length(exclude_gsms) ||
      anyDuplicated(selected$gsm)) {
    stop(
      "Expected one GSE49543 sample for each requested GSM: ",
      paste(exclude_gsms, collapse = ", ")
    )
  }
  selected <- selected[match(exclude_gsms, selected$gsm), , drop = FALSE]
  exclusion_record <- data.frame(
    run_label = rep(RUN_LABEL, nrow(selected)),
    dataset = selected$dataset,
    tissue = selected$tissue,
    gsm = selected$gsm,
    animal_id = selected$animal_id,
    adjudicated_group = selected$adjudicated_group,
    expression_column = selected$expression_column,
    reason = paste0(
      "Cochlear sensitivity analysis excluding ",
      nrow(selected),
      " observed GSM record(s): ",
      paste(selected$gsm, collapse = ", "),
      "."
    ),
    stringsAsFactors = FALSE
  )
  tissue_matrices$cochlea <- tissue_matrices$cochlea[
    ,
    !colnames(tissue_matrices$cochlea) %in% selected$expression_column,
    drop = FALSE
  ]
  expected_cochlea_n <- 41 - nrow(selected)
  if (ncol(tissue_matrices$cochlea) != expected_cochlea_n) {
    stop(
      "The cochlear sensitivity matrix should contain ",
      expected_cochlea_n,
      " samples after exclusion; observed ",
      ncol(tissue_matrices$cochlea)
    )
  }
}
write.csv(
  exclusion_record,
  file.path(RESULTS, "sample_exclusion.csv"),
  row.names = FALSE
)

if (FREEZE_GENE_UNIVERSE) {
  frozen_universe_path <- file.path(PRIMARY_RESULTS, "analysis_universe.csv")
  if (!file.exists(frozen_universe_path)) {
    stop("Frozen primary gene universe not found: ", frozen_universe_path)
  }
  common_genes <- read.csv(
    frozen_universe_path,
    stringsAsFactors = FALSE
  )$gene
  unavailable_frozen <- vapply(
    tissue_matrices,
    function(matrix) any(!common_genes %in% rownames(matrix)),
    logical(1)
  )
  if (any(unavailable_frozen)) {
    stop("A frozen analysis gene is unavailable in a sensitivity matrix.")
  }
  invariant_genes <- do.call(
    rbind,
    lapply(names(tissue_matrices), function(tissue) {
      matrix <- tissue_matrices[[tissue]][common_genes, , drop = FALSE]
      invariant <- apply(matrix, 1, var) == 0
      if (!any(invariant)) {
        return(data.frame(
          tissue = character(),
          gene = character(),
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        tissue = rep(tissue, sum(invariant)),
        gene = names(invariant)[invariant],
        stringsAsFactors = FALSE
      )
    })
  )
  write.csv(
    invariant_genes,
    file.path(RESULTS, "frozen_universe_invariant_genes.csv"),
    row.names = FALSE
  )
} else {
  valid_by_tissue <- lapply(tissue_matrices, function(matrix) {
    apply(matrix, 1, function(x) all(is.finite(x)) && var(x) > 0)
  })
  common_genes <- Reduce(
    intersect,
    lapply(valid_by_tissue, function(valid) names(valid)[valid])
  )
  common_genes <- sort(common_genes)
}
if (length(common_genes) == 0) {
  stop("No common valid genes across tissues.")
}

universe <- data.frame(gene = common_genes, stringsAsFactors = FALSE)
for (tissue in names(tissue_matrices)) {
  matrix <- tissue_matrices[[tissue]][common_genes, , drop = FALSE]
  universe[[paste0(tissue, "_variance")]] <- apply(matrix, 1, var)
  universe[[paste0(tissue, "_min")]] <- apply(matrix, 1, min)
  universe[[paste0(tissue, "_max")]] <- apply(matrix, 1, max)
}
write.csv(
  universe,
  file.path(RESULTS, "analysis_universe.csv"),
  row.names = FALSE
)

counts <- list()
significant_long <- list()
shared_list <- list()

for (tissue in names(TISSUES)) {
  config <- TISSUES[[tissue]]
  matrix <- tissue_matrices[[tissue]][common_genes, , drop = FALSE]
  metadata <- load_metadata(matrix, config$dataset)
  design <- build_design(metadata, tissue, config$dataset)
  contrasts <- contrast_matrix_for(design)
  fit <- lmFit(matrix, design)
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
    table$gene <- rownames(table)
    table <- table[, c(
      "gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B"
    )]
    table$significant <- with(
      table,
      adj.P.Val < FDR_CUTOFF & abs(logFC) >= ABS_LOG2FC_CUTOFF
    )
    table$direction <- ifelse(
      table$significant & table$logFC > 0,
      "up",
      ifelse(table$significant & table$logFC < 0, "down", "not_significant")
    )
    write.csv(
      table,
      file.path(LIMMA_DIR, paste0(tissue, "_", contrast_id, ".csv")),
      row.names = FALSE
    )

    counts[[paste(tissue, contrast_id, sep = "__")]] <- data.frame(
      dataset = config$dataset,
      tissue = tissue,
      contrast = contrast_id,
      numerator = CONTRASTS$numerator[CONTRASTS$contrast_id == contrast_id],
      denominator = CONTRASTS$denominator[CONTRASTS$contrast_id == contrast_id],
      n_samples = nrow(metadata),
      genes_tested = nrow(table),
      n_up = sum(table$direction == "up"),
      n_down = sum(table$direction == "down"),
      n_fdr_only = sum(table$adj.P.Val < FDR_CUTOFF),
      n_effect_size_only = sum(abs(table$logFC) >= ABS_LOG2FC_CUTOFF),
      n_significant = sum(table$significant),
      fdr_cutoff = FDR_CUTOFF,
      abs_log2fc_cutoff = ABS_LOG2FC_CUTOFF
    )

    significant <- table[table$significant, , drop = FALSE]
    if (nrow(significant) > 0) {
      significant_long[[paste(tissue, contrast_id, sep = "__")]] <- data.frame(
        dataset = config$dataset,
        tissue = tissue,
        contrast = contrast_id,
        significant[, c("gene", "logFC", "adj.P.Val", "direction")],
        row.names = NULL
      )
    } else {
      significant_long[[paste(tissue, contrast_id, sep = "__")]] <- data.frame(
        dataset = character(),
        tissue = character(),
        contrast = character(),
        gene = character(),
        logFC = numeric(),
        adj.P.Val = numeric(),
        direction = character()
      )
    }
  }
}

counts_table <- do.call(rbind, counts)
rownames(counts_table) <- NULL
write.csv(
  counts_table,
  file.path(RESULTS, "significant_counts.csv"),
  row.names = FALSE
)
write.csv(
  do.call(rbind, significant_long),
  file.path(RESULTS, "significant_genes_long.csv"),
  row.names = FALSE
)

for (contrast_id in CONTRASTS$contrast_id) {
  cochlea_file <- file.path(
    LIMMA_DIR, paste0("cochlea_", contrast_id, ".csv")
  )
  ic_file <- file.path(
    LIMMA_DIR, paste0("inferior_colliculus_", contrast_id, ".csv")
  )
  cochlea_table <- read.csv(cochlea_file, stringsAsFactors = FALSE)
  ic_table <- read.csv(ic_file, stringsAsFactors = FALSE)
  shared <- merge(
    cochlea_table[
      cochlea_table$significant,
      c("gene", "logFC", "adj.P.Val", "direction")
    ],
    ic_table[
      ic_table$significant,
      c("gene", "logFC", "adj.P.Val", "direction")
    ],
    by = "gene",
    suffixes = c("_cochlea", "_inferior_colliculus")
  )
  if (nrow(shared) > 0) {
    shared <- data.frame(
      contrast = contrast_id,
      shared,
      same_direction = shared$direction_cochlea == shared$direction_inferior_colliculus,
      row.names = NULL
    )
  } else {
    shared <- data.frame(
      contrast = character(),
      gene = character(),
      logFC_cochlea = numeric(),
      adj.P.Val_cochlea = numeric(),
      direction_cochlea = character(),
      logFC_inferior_colliculus = numeric(),
      adj.P.Val_inferior_colliculus = numeric(),
      direction_inferior_colliculus = character(),
      same_direction = logical()
    )
  }
  shared_list[[contrast_id]] <- shared
}
shared_genes <- do.call(rbind, shared_list)
rownames(shared_genes) <- NULL
write.csv(
  shared_genes,
  file.path(RESULTS, "shared_genes.csv"),
  row.names = FALSE
)

shared_counts <- do.call(
  rbind,
  lapply(CONTRASTS$contrast_id, function(contrast_id) {
    selected <- shared_genes[shared_genes$contrast == contrast_id, , drop = FALSE]
    data.frame(
      contrast = contrast_id,
      n_shared = nrow(selected),
      n_same_direction = sum(selected$same_direction),
      n_opposite_direction = sum(!selected$same_direction)
    )
  })
)
write.csv(
  shared_counts,
  file.path(RESULTS, "shared_gene_counts.csv"),
  row.names = FALSE
)

writeLines(
  "~ 0 + cohort_group + sex",
  file.path(RESULTS, "model_formula.txt")
)
write.csv(
  CONTRASTS,
  file.path(RESULTS, "contrast_table.csv"),
  row.names = FALSE
)

input_files <- c(
  file.path(ANALYSIS_ROOT, "expression_GSE49543_log2.csv"),
  file.path(ANALYSIS_ROOT, "expression_GSE49522_log2.csv"),
  SAMPLE_MASTER
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
  file.path(RESULTS, "input_checksums.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    version = paste("V14 limma", RUN_LABEL),
    model = "~ 0 + cohort_group + sex",
    contrasts = CONTRASTS,
    run_label = RUN_LABEL,
    exclusion = exclusion_record,
    frozen_gene_universe = FREEZE_GENE_UNIVERSE,
    fdr_cutoff = FDR_CUTOFF,
    abs_log2fc_cutoff = ABS_LOG2FC_CUTOFF,
    common_genes = common_genes,
    session_info = sessionInfo()
  ),
  file.path(RESULTS, "limma_run.rds")
)

writeLines(
  capture.output(sessionInfo()),
  file.path(RESULTS, "sessionInfo.txt")
)

limma_metadata <- list(
  generated_at = format(Sys.time(), "%Y-%m-%dT%H:%M:%S%z"),
  run_label = RUN_LABEL,
  model_formula = "~ 0 + cohort_group + sex",
  sample_exclusion = exclusion_record,
  frozen_gene_universe = FREEZE_GENE_UNIVERSE,
  cohort_levels = GROUPS,
  sex_levels = SEXES,
  contrasts = CONTRASTS,
  expression_scale = "log2 values supplied by input matrices; no additional log transform",
  fdr_method = "Benjamini-Hochberg within tissue and contrast",
  fdr_cutoff = FDR_CUTOFF,
  abs_log2fc_cutoff = ABS_LOG2FC_CUTOFF,
  limma = list(
    lmFit = "limma::lmFit",
    eBayes = "limma::eBayes(fit, trend = TRUE, robust = TRUE)"
  ),
  common_analysis_genes = length(common_genes),
  supported_by_continuous_abr = FALSE,
  continuous_abr_reason = "No per-animal ABR values were deposited or recovered.",
  dpoae_used = FALSE,
  dpoae_reason = "No per-animal DPOAE values were deposited or recovered; no composite phenotype was created.",
  input_checksums = input_checksums,
  outputs = list(
    significant_counts = "results/significant_counts.csv",
    significant_genes_long = "results/significant_genes_long.csv",
    shared_genes = "results/shared_genes.csv",
    shared_gene_counts = "results/shared_gene_counts.csv",
    contrast_table = "results/contrast_table.csv",
    model_formula = "results/model_formula.txt",
    analysis_universe = "results/analysis_universe.csv",
    design_matrices = "results/design/",
    limma_tables = "results/limma/"
  ),
  session_info = capture.output(sessionInfo())
)
write_json(
  limma_metadata,
  file.path(RESULTS, "limma_run_metadata.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat("V14 limma analysis complete:", RUN_LABEL, "\n")
cat("Common analysis genes:", length(common_genes), "\n")
print(counts_table[, c(
  "tissue", "contrast", "n_fdr_only", "n_effect_size_only",
  "n_up", "n_down", "n_significant"
)], row.names = FALSE)
print(shared_counts, row.names = FALSE)
