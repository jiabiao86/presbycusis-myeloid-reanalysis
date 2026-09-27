#!/usr/bin/env Rscript

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
lib <- "/Users/jijiabiao/Desktop/耳聋课题/Rlibs"
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages(library(limma))

out <- file.path(root, "sensitivity_analysis")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

read_expression <- function(path) {
  frame <- read.csv(path, check.names = FALSE, row.names = 1)
  matrix <- as.matrix(frame)
  storage.mode(matrix) <- "double"
  matrix
}

read_series_metadata <- function(path) {
  lines <- readLines(gzfile(path), warn = FALSE)
  get_field <- function(prefix) {
    line <- lines[startsWith(lines, prefix)][1]
    if (is.na(line)) {
      return(character())
    }
    values <- strsplit(line, "\t", fixed = TRUE)[[1]][-1]
    gsub('^"|"$', "", values)
  }
  titles <- get_field("!Sample_title")
  accessions <- get_field("!Sample_geo_accession")
  characteristics <- lines[startsWith(lines, "!Sample_characteristics_ch1")]
  values <- lapply(characteristics, function(line) {
    gsub('^"|"$', "", strsplit(line, "\t", fixed = TRUE)[[1]][-1])
  })
  sex <- rep(NA_character_, length(accessions))
  hearing <- rep(NA_character_, length(accessions))
  for (value in values) {
    sex[grepl("^gender:", value)] <- sub("^gender: ", "", value[grepl("^gender:", value)])
    hearing[grepl("^hearing status:", value)] <- sub(
      "^hearing status: ",
      "",
      value[grepl("^hearing status:", value)]
    )
  }
  data.frame(
    gsm = accessions,
    title = titles,
    sex = sex,
    hearing_status = hearing,
    stringsAsFactors = FALSE
  )
}

cochlea <- read_expression(file.path(root, "expression_GSE49543_log2.csv"))
central <- read_expression(file.path(root, "expression_GSE49522_log2.csv"))

cochlea_metadata <- read.csv(
  file.path(root, "GSE49543_sample_metadata.csv"),
  encoding = "utf-8-sig",
  stringsAsFactors = FALSE
)
central_metadata <- read_series_metadata(
  "/tmp/ear_study_learning/GSE49522_series_matrix.txt.gz"
)

severity_map <- c("YC" = 0, "MA" = 1, "MP" = 2, "SP" = 3)
cochlea_group <- sub("-.*$", "", colnames(cochlea))
cochlea_severity <- severity_map[cochlea_group]
cochlea_sex <- cochlea_metadata$gender
if (length(cochlea_sex) != ncol(cochlea)) {
  stop("Cochlear metadata and expression columns do not match")
}

central_metadata <- central_metadata[match(colnames(central), central_metadata$gsm), ]
central_sex <- central_metadata$sex
central_group_code <- c(
  "Young Control" = "YC",
  "Middle-Aged" = "MA",
  "Mild Presbycusis" = "MP",
  "Severe Presbycusis" = "SP"
)[central_metadata$hearing_status]
central_group <- central_group_code
central_severity <- severity_map[central_group]
cat("Central metadata rows:", nrow(central_metadata), "NA sex:", sum(is.na(central_sex)), "\n")
print(head(data.frame(sample = colnames(central), gsm = central_metadata$gsm, sex = central_sex)))

fit_model <- function(expression, design, coefficient) {
  fit <- eBayes(lmFit(expression, design))
  result <- topTable(fit, coef = coefficient, number = Inf, sort.by = "P")
  result$gene <- rownames(result)
  result[, c("gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "B")]
}

cochlea_sex_design <- model.matrix(~ cochlea_severity + factor(cochlea_sex))
cochlea_age_design <- model.matrix(
  ~ cochlea_severity + I(cochlea_group %in% c("MP", "SP"))
)
central_sex_design <- model.matrix(~ central_severity + factor(central_sex))
central_age_design <- model.matrix(
  ~ central_severity + I(central_group %in% c("MP", "SP"))
)

cat("Cochlear sex design:", paste(dim(cochlea_sex_design), collapse = "x"), "\n")
cat("Cochlear age design:", paste(dim(cochlea_age_design), collapse = "x"), "\n")
cat("Central sex design:", paste(dim(central_sex_design), collapse = "x"), "\n")
cat("Central age design:", paste(dim(central_age_design), collapse = "x"), "\n")

cochlea_sex_result <- fit_model(cochlea, cochlea_sex_design, "cochlea_severity")
cochlea_age_result <- fit_model(cochlea, cochlea_age_design, "cochlea_severity")
central_sex_result <- fit_model(central, central_sex_design, "central_severity")
central_age_result <- fit_model(central, central_age_design, "central_severity")

write.csv(
  cochlea_sex_result,
  file.path(out, "cochlea_severity_sex_adjusted.csv"),
  row.names = FALSE
)
write.csv(
  cochlea_age_result,
  file.path(out, "cochlea_severity_age_adjusted.csv"),
  row.names = FALSE
)
write.csv(
  central_sex_result,
  file.path(out, "central_severity_sex_adjusted.csv"),
  row.names = FALSE
)
write.csv(
  central_age_result,
  file.path(out, "central_severity_age_adjusted.csv"),
  row.names = FALSE
)

candidate_genes <- c(
  "H2-Aa", "H2-Eb1", "Cd74", "Setd1a", "Ctss", "Fcgr3",
  "Cd68", "Tyrobp", "Mpeg1", "Lgals3", "C1qa", "C1qb",
  "C1qc", "C3ar1", "Ms4a7", "Csf1r", "Clec7a", "Clec4d"
)

cochlea_baseline <- read.csv(
  file.path(root, "limma_severity_trend.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
central_baseline <- read.csv(
  file.path(root, "gse49522_severity_trend.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)

cochlea_baseline_candidate <- cochlea_baseline[
  match(candidate_genes, cochlea_baseline$gene),
  c("gene", "logFC", "adj.P.Val")
]
cochlea_sex_candidate <- cochlea_sex_result[
  match(candidate_genes, cochlea_sex_result$gene),
  c("gene", "logFC", "adj.P.Val")
]
cochlea_age_candidate <- cochlea_age_result[
  match(candidate_genes, cochlea_age_result$gene),
  c("gene", "logFC", "adj.P.Val")
]
central_baseline_candidate <- central_baseline[
  match(candidate_genes, central_baseline$gene),
  c("gene", "central_severity_correlation", "central_fdr_bh")
]
central_sex_candidate <- central_sex_result[
  match(candidate_genes, central_sex_result$gene),
  c("gene", "logFC", "adj.P.Val")
]
central_age_candidate <- central_age_result[
  match(candidate_genes, central_age_result$gene),
  c("gene", "logFC", "adj.P.Val")
]

candidate_table <- data.frame(
  gene = candidate_genes,
  cochlea_baseline_logFC = cochlea_baseline_candidate$logFC,
  cochlea_baseline_fdr = cochlea_baseline_candidate$adj.P.Val,
  cochlea_sex_adjusted_logFC = cochlea_sex_candidate$logFC,
  cochlea_sex_adjusted_fdr = cochlea_sex_candidate$adj.P.Val,
  cochlea_age_adjusted_logFC = cochlea_age_candidate$logFC,
  cochlea_age_adjusted_fdr = cochlea_age_candidate$adj.P.Val,
  central_baseline_r = central_baseline_candidate$central_severity_correlation,
  central_baseline_fdr = central_baseline_candidate$central_fdr_bh,
  central_sex_adjusted_logFC = central_sex_candidate$logFC,
  central_sex_adjusted_fdr = central_sex_candidate$adj.P.Val,
  central_age_adjusted_logFC = central_age_candidate$logFC,
  central_age_adjusted_fdr = central_age_candidate$adj.P.Val,
  stringsAsFactors = FALSE
)
candidate_table$cochlea_sex_direction_concordant <-
  sign(candidate_table$cochlea_baseline_logFC) ==
    sign(candidate_table$cochlea_sex_adjusted_logFC)
candidate_table$cochlea_age_direction_concordant <-
  sign(candidate_table$cochlea_baseline_logFC) ==
    sign(candidate_table$cochlea_age_adjusted_logFC)
candidate_table$central_sex_direction_concordant <-
  sign(candidate_table$central_baseline_r) ==
    sign(candidate_table$central_sex_adjusted_logFC)
candidate_table$central_age_direction_concordant <-
  sign(candidate_table$central_baseline_r) ==
    sign(candidate_table$central_age_adjusted_logFC)

write.csv(
  candidate_table,
  file.path(out, "candidate_models_sensitivity.csv"),
  row.names = FALSE
)

set.seed(2026)
n_bootstrap <- 500
bootstrap_candidates <- function(expression, severity, sex, genes) {
  result <- matrix(NA_real_, nrow = length(genes), ncol = 4)
  rownames(result) <- genes
  colnames(result) <- c("lower", "median", "upper", "positive_fraction")
  for (gene_index in seq_along(genes)) {
    gene <- genes[gene_index]
    if (!gene %in% rownames(expression)) {
      next
    }
    values <- as.numeric(expression[gene, ])
    estimates <- numeric(n_bootstrap)
    for (bootstrap_index in seq_len(n_bootstrap)) {
      sampled <- unlist(lapply(sort(unique(severity)), function(level) {
        indices <- which(severity == level)
        sample(indices, length(indices), replace = TRUE)
      }))
      data <- data.frame(
        expression = values[sampled],
        severity = severity[sampled],
        sex = factor(sex[sampled])
      )
      model <- lm(expression ~ severity + sex, data = data)
      estimates[bootstrap_index] <- coef(model)[["severity"]]
    }
    result[gene, ] <- c(
      quantile(estimates, 0.025, na.rm = TRUE),
      median(estimates, na.rm = TRUE),
      quantile(estimates, 0.975, na.rm = TRUE),
      mean(estimates > 0, na.rm = TRUE)
    )
  }
  as.data.frame(result)
}

bootstrap_cochlea <- bootstrap_candidates(
  cochlea,
  cochlea_severity,
  cochlea_sex,
  candidate_genes
)
bootstrap_central <- bootstrap_candidates(
  central,
  central_severity,
  central_sex,
  candidate_genes
)
bootstrap_cochlea$tissue <- "Cochlea"
bootstrap_cochlea$gene <- rownames(bootstrap_cochlea)
bootstrap_central$tissue <- "Inferior colliculus"
bootstrap_central$gene <- rownames(bootstrap_central)
bootstrap_table <- rbind(bootstrap_cochlea, bootstrap_central)
write.csv(
  bootstrap_table,
  file.path(out, "candidate_severity_bootstrap.csv"),
  row.names = FALSE
)

summary <- list(
  cochlea = list(
    severity_only_fdr005 = sum(cochlea_baseline$adj.P.Val < 0.05),
    sex_adjusted_fdr005 = sum(cochlea_sex_result$adj.P.Val < 0.05),
    age_adjusted_fdr005 = sum(cochlea_age_result$adj.P.Val < 0.05),
    sex_adjusted_effect_correlation = cor(
      cochlea_baseline$logFC[match(cochlea_sex_result$gene, cochlea_baseline$gene)],
      cochlea_sex_result$logFC,
      use = "complete.obs"
    ),
    age_adjusted_effect_correlation = cor(
      cochlea_baseline$logFC[match(cochlea_age_result$gene, cochlea_baseline$gene)],
      cochlea_age_result$logFC,
      use = "complete.obs"
    )
  ),
  central = list(
    severity_only_fdr005 = sum(central_baseline$central_fdr_bh < 0.05),
    sex_adjusted_fdr005 = sum(central_sex_result$adj.P.Val < 0.05),
    age_adjusted_fdr005 = sum(central_age_result$adj.P.Val < 0.05),
    sex_adjusted_effect_correlation = cor(
      central_baseline$central_severity_correlation[match(central_sex_result$gene, central_baseline$gene)],
      central_sex_result$logFC,
      use = "complete.obs"
    ),
    age_adjusted_effect_correlation = cor(
      central_baseline$central_severity_correlation[match(central_age_result$gene, central_baseline$gene)],
      central_age_result$logFC,
      use = "complete.obs"
    )
  )
)

writeLines(
  jsonlite::toJSON(summary, auto_unbox = TRUE, pretty = TRUE),
  file.path(out, "sensitivity_summary.json")
)

cat("Sensitivity analysis complete:", out, "\n")
