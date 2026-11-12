#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run with Rscript so repository paths can be resolved.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."))
prepared_dir <- file.path(repo_root, "inputs", "external", "prepared")
outdir <- file.path(repo_root, "results", "external")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

candidate_genes <- c(
  "Cd68", "Clec7a", "Fcgr3", "Lyz1", "Lyz2", "Mpeg1", "Ptprc"
)

read_expression <- function(accession) {
  path <- file.path(
    prepared_dir, paste0(accession, "_expression_log2.csv")
  )
  data <- read.csv(path, check.names = FALSE, stringsAsFactors = FALSE)
  if (!"gene" %in% names(data)) {
    stop(accession, ": missing gene column")
  }
  data$gene <- as.character(data$gene)
  if (anyDuplicated(data$gene)) {
    stop(accession, ": duplicate gene symbols in prepared matrix")
  }
  data
}

format_p <- function(values) {
  ifelse(
    is.na(values),
    NA_character_,
    ifelse(values == 0, "<2.2e-16", formatC(values, digits = 15, format = "g"))
  )
}

p_status <- function(values) {
  ifelse(
    is.na(values),
    "not_estimable",
    ifelse(values == 0, "numeric_underflow", "finite")
  )
}

welch_scan <- function(
  data,
  target_columns,
  control_columns,
  dataset,
  context,
  cell_type,
  comparison
) {
  missing_columns <- setdiff(c(target_columns, control_columns), names(data))
  if (length(missing_columns) > 0) {
    stop(context, ": missing columns ", paste(missing_columns, collapse = ", "))
  }
  target <- as.matrix(data[, target_columns, drop = FALSE])
  storage.mode(target) <- "double"
  control <- as.matrix(data[, control_columns, drop = FALSE])
  storage.mode(control) <- "double"
  if (any(!is.finite(target)) || any(!is.finite(control))) {
    stop(context, ": prepared expression matrix contains non-finite values")
  }

  n_target <- ncol(target)
  n_control <- ncol(control)
  mean_target <- rowMeans(target)
  mean_control <- rowMeans(control)
  variance_target <- apply(target, 1, var)
  variance_control <- apply(control, 1, var)
  squared_error <- variance_target / n_target + variance_control / n_control
  delta <- mean_target - mean_control
  t_value <- ifelse(
    squared_error > 0,
    delta / sqrt(squared_error),
    ifelse(delta == 0, 0, sign(delta) * Inf)
  )
  numerator <- squared_error^2
  denominator <-
    (variance_target / n_target)^2 / (n_target - 1) +
    (variance_control / n_control)^2 / (n_control - 1)
  degrees_freedom <- ifelse(
    squared_error > 0 & denominator > 0,
    numerator / denominator,
    Inf
  )
  p_value <- ifelse(
    is.infinite(t_value),
    0,
    2 * pt(-abs(t_value), df = degrees_freedom)
  )
  p_value <- pmin(pmax(p_value, 0), 1)
  fdr_bh <- p.adjust(p_value, method = "BH")

  result <- data.frame(
    dataset = dataset,
    context = context,
    cell_type = cell_type,
    comparison = comparison,
    gene = data$gene,
    n_target = n_target,
    n_control = n_control,
    mean_target_log2 = mean_target,
    mean_control_log2 = mean_control,
    log2FC = delta,
    fold_change = 2^delta,
    welch_t = t_value,
    welch_df = degrees_freedom,
    p_value = p_value,
    p_value_status = p_status(p_value),
    p_value_display = format_p(p_value),
    fdr_bh = fdr_bh,
    fdr_bh_status = p_status(fdr_bh),
    fdr_bh_display = format_p(fdr_bh),
    tested_transcripts = nrow(data),
    normalization = "submitter-normalized values; log2(x + 1)",
    stringsAsFactors = FALSE
  )
  result[
    order(result$fdr_bh, result$p_value, result$log2FC, result$gene),
    ,
    drop = FALSE
  ]
}

write_context <- function(result, stem) {
  write.csv(
    result,
    file.path(outdir, paste0(stem, "_all_genes.csv")),
    row.names = FALSE
  )
  candidate <- result[match(candidate_genes, result$gene), , drop = FALSE]
  missing_index <- is.na(candidate$gene)
  candidate$gene <- candidate_genes
  candidate$dataset[missing_index] <- result$dataset[1]
  candidate$context[missing_index] <- result$context[1]
  candidate$cell_type[missing_index] <- result$cell_type[1]
  candidate$comparison[missing_index] <- result$comparison[1]
  candidate$n_target[missing_index] <- result$n_target[1]
  candidate$n_control[missing_index] <- result$n_control[1]
  candidate$tested_transcripts[missing_index] <- result$tested_transcripts[1]
  candidate$normalization[missing_index] <- result$normalization[1]
  candidate$p_value_status[missing_index] <- "gene_not_available"
  candidate$fdr_bh_status[missing_index] <- "gene_not_available"
  write.csv(
    candidate,
    file.path(outdir, paste0(stem, "_candidate_genes.csv")),
    row.names = FALSE
  )
  candidate
}

gse233798 <- read_expression("GSE233798")
res_233798 <- welch_scan(
  gse233798,
  target_columns = grep("^fpkm_Base_old", names(gse233798), value = TRUE),
  control_columns = grep("^fpkm_Base_young", names(gse233798), value = TRUE),
  dataset = "GSE233798",
  context = "GSE233798 cochlea",
  cell_type = "whole cochlea",
  comparison = "old_vs_young"
)
cand_233798 <- write_context(res_233798, "GSE233798")

gse153882 <- read_expression("GSE153882")
ihc_target <- paste0("26m_IHC_", 1:4)
ihc_control <- paste0("9m_IHC_", 1:4)
ohc_target <- paste0("26m_OHC_", 1:4)
ohc_control <- paste0("9m_OHC_", 1:4)
res_ihc <- welch_scan(
  gse153882,
  target_columns = ihc_target,
  control_columns = ihc_control,
  dataset = "GSE153882",
  context = "GSE153882 IHC",
  cell_type = "inner hair cell",
  comparison = "26m_vs_9m"
)
cand_ihc <- write_context(res_ihc, "GSE153882_IHC")
res_ohc <- welch_scan(
  gse153882,
  target_columns = ohc_target,
  control_columns = ohc_control,
  dataset = "GSE153882",
  context = "GSE153882 OHC",
  cell_type = "outer hair cell",
  comparison = "26m_vs_9m"
)
cand_ohc <- write_context(res_ohc, "GSE153882_OHC")

gse154833 <- read_expression("GSE154833")
res_154833 <- welch_scan(
  gse154833,
  target_columns = paste0("26mo_SV_", c("A", "B", "C")),
  control_columns = paste0("1mo_SV_", c("A", "B", "C")),
  dataset = "GSE154833",
  context = "GSE154833 stria vascularis",
  cell_type = "stria vascularis",
  comparison = "26m_vs_1m"
)
cand_154833 <- write_context(res_154833, "GSE154833")

all_results <- rbind(res_233798, res_ihc, res_ohc, res_154833)
candidate_results <- rbind(cand_233798, cand_ihc, cand_ohc, cand_154833)
write.csv(
  candidate_results,
  file.path(outdir, "external_candidate_genes.csv"),
  row.names = FALSE
)

contexts <- unique(all_results[, c("dataset", "context", "comparison")])
summary_rows <- lapply(seq_len(nrow(contexts)), function(index) {
  current <- all_results[all_results$context == contexts$context[index], ]
  data.frame(
    dataset = contexts$dataset[index],
    context = contexts$context[index],
    comparison = contexts$comparison[index],
    tested_transcripts = nrow(current),
    nominal_p_lt_0_05 = sum(current$p_value < 0.05),
    bh_fdr_lt_0_05 = sum(current$fdr_bh < 0.05),
    p_value_numeric_underflow = sum(current$p_value_status == "numeric_underflow"),
    fdr_numeric_underflow = sum(current$fdr_bh_status == "numeric_underflow"),
    method = "Welch t test; Welch-Satterthwaite df; BH across tested transcripts",
    stringsAsFactors = FALSE
  )
})
summary <- do.call(rbind, summary_rows)
write.csv(
  summary,
  file.path(outdir, "external_welch_summary.csv"),
  row.names = FALSE
)

underflow <- all_results[all_results$p_value_status == "numeric_underflow", ]
underflow <- underflow[, c(
  "dataset", "context", "gene", "log2FC", "welch_t", "welch_df",
  "p_value", "p_value_status", "p_value_display",
  "fdr_bh", "fdr_bh_status", "fdr_bh_display"
)]
write.csv(
  underflow,
  file.path(outdir, "external_p_value_underflow.csv"),
  row.names = FALSE
)

readme <- c(
  "# External Welch and BH results",
  "",
  "This directory contains the complete transcriptome-wide Welch tests and",
  "Benjamini-Hochberg adjusted values for the three external GEO datasets.",
  "",
  "- GSE233798: old versus young whole cochlea; 3 versus 3 biological samples.",
  "- GSE153882: 26-month versus 9-month IHC and OHC; 4 versus 4 biological samples each.",
  "- GSE154833: 26-month versus 1-month stria vascularis; 3 versus 3 biological samples.",
  "",
  "The submitters supplied normalized FPKM or RPKM values. The analysis applies",
  "log2(x + 1), resolves duplicate symbols by the highest mean submitted value,",
  "and performs Welch t tests with Welch-Satterthwaite degrees of freedom.",
  "BH correction is applied across every tested transcript within each stated",
  "dataset/context comparison. No external P value changes the primary analysis.",
  "",
  "GSE154833 has 36 GSMs representing four sequencing lanes for each of nine",
  "biological samples. The submitted All values worksheet already supplies nine",
  "biological-sample columns (A/B/C at each age); lanes are not treated as",
  "independent replicates. The 36-run mapping is retained in",
  "inputs/external/prepared/run_mapping_external.csv.",
  "",
  "When an extreme Welch statistic exceeds the floating-point resolution of the",
  "student t distribution, R returns p = 0. Such rows are retained without",
  "imputation, labelled p_value_status = numeric_underflow, and displayed as",
  "<2.2e-16. BH correction uses the numeric p values as returned by R. The exact",
  "affected rows are traceable in external_p_value_underflow.csv.",
  "",
  "Files:",
  "- *_all_genes.csv: complete results for each dataset/context.",
  "- *_candidate_genes.csv: the seven primary intersection genes; absent genes",
  "  are retained as gene_not_available rather than imputed.",
  "- external_candidate_genes.csv: combined candidate table for Figure S2.",
  "- external_welch_summary.csv: transcript counts and significance totals.",
  "- external_p_value_underflow.csv: traceable rows with numeric p = 0.",
  "- external_welch_validation*.csv and .json: independent R t.test checks for",
  "  all candidate genes and p.adjust checks across every tested transcript."
)
writeLines(readme, file.path(outdir, "EXTERNAL_RESULTS_README.md"))
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))

manifest <- list(
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  script = "code/R/run_external_welch_bh.R",
  method = "Welch t test with Welch-Satterthwaite df; BH correction",
  contexts = summary,
  candidate_genes = candidate_genes,
  p_value_zero_policy = paste(
    "retain numeric zero; label numeric_underflow; display <2.2e-16;",
    "BH uses the numeric values and zero rows are exported for traceability"
  )
)
jsonlite_available <- requireNamespace("jsonlite", quietly = TRUE)
if (jsonlite_available) {
  jsonlite::write_json(
    manifest,
    file.path(outdir, "external_welch_manifest.json"),
    pretty = TRUE,
    auto_unbox = TRUE,
    na = "null"
  )
}

print(summary)
