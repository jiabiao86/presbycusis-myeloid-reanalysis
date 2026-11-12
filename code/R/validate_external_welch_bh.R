#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run with Rscript so repository paths can be resolved.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
repo_root <- normalizePath(file.path(dirname(script_file), "..", ".."))
prepared_dir <- file.path(repo_root, "inputs", "external", "prepared")
results_dir <- file.path(repo_root, "results", "external")

candidate_genes <- c(
  "Cd68", "Clec7a", "Fcgr3", "Lyz1", "Lyz2", "Mpeg1", "Ptprc"
)

contexts <- list(
  list(
    accession = "GSE233798",
    stem = "GSE233798",
    target = grep("^fpkm_Base_old", names(
      read.csv(file.path(prepared_dir, "GSE233798_expression_log2.csv"),
               nrows = 1, check.names = FALSE)
    ), value = TRUE),
    control = grep("^fpkm_Base_young", names(
      read.csv(file.path(prepared_dir, "GSE233798_expression_log2.csv"),
               nrows = 1, check.names = FALSE)
    ), value = TRUE)
  ),
  list(
    accession = "GSE153882",
    stem = "GSE153882_IHC",
    target = paste0("26m_IHC_", 1:4),
    control = paste0("9m_IHC_", 1:4)
  ),
  list(
    accession = "GSE153882",
    stem = "GSE153882_OHC",
    target = paste0("26m_OHC_", 1:4),
    control = paste0("9m_OHC_", 1:4)
  ),
  list(
    accession = "GSE154833",
    stem = "GSE154833",
    target = paste0("26mo_SV_", c("A", "B", "C")),
    control = paste0("1mo_SV_", c("A", "B", "C"))
  )
)

same_number <- function(left, right, tolerance = 1e-10) {
  if (length(left) != 1 || length(right) != 1) {
    return(FALSE)
  }
  if (is.na(left) || is.na(right)) {
    return(is.na(left) && is.na(right))
  }
  isTRUE(all.equal(left, right, tolerance = tolerance, check.attributes = FALSE))
}

rows <- list()
summaries <- list()

for (context in contexts) {
  expression <- read.csv(
    file.path(prepared_dir, paste0(context$accession, "_expression_log2.csv")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  results <- read.csv(
    file.path(results_dir, paste0(context$stem, "_all_genes.csv")),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  missing_columns <- setdiff(
    c("gene", context$target, context$control),
    names(expression)
  )
  if (length(missing_columns) > 0) {
    stop(context$stem, ": missing expression columns ",
         paste(missing_columns, collapse = ", "))
  }
  if (anyDuplicated(expression$gene) || anyDuplicated(results$gene)) {
    stop(context$stem, ": duplicate gene symbols")
  }
  if (length(expression$gene) != length(results$gene) ||
      !setequal(expression$gene, results$gene)) {
    stop(context$stem, ": expression and result gene sets differ")
  }

  context_rows <- lapply(candidate_genes, function(gene) {
    gene_index <- match(gene, expression$gene)
    result_index <- match(gene, results$gene)
    if (is.na(gene_index) || is.na(result_index)) {
      return(data.frame(
        context = context$stem,
        gene = gene,
        status = "gene_not_available",
        target_mean_check = NA,
        control_mean_check = NA,
        log2fc_check = NA,
        welch_t_check = NA,
        welch_df_check = NA,
        p_value_check = NA,
        check_method = "gene not available",
        passed = TRUE
      ))
    }

    target_values <- as.numeric(expression[gene_index, context$target])
    control_values <- as.numeric(expression[gene_index, context$control])
    result_row <- results[result_index, ]
    expected_delta <- mean(target_values) - mean(control_values)
    constant_zero_null <-
      expected_delta == 0 &&
      var(target_values) == 0 &&
      var(control_values) == 0
    if (constant_zero_null) {
      checked_t <- 0
      checked_df <- Inf
      checked_p <- 1
      check_method <- "R t.test was undefined for two constant-zero groups; exact null rule used"
    } else {
      test <- t.test(target_values, control_values, var.equal = FALSE)
      checked_t <- unname(test$statistic)
      checked_df <- unname(test$parameter)
      checked_p <- test$p.value
      check_method <- "R t.test"
    }
    delta_match <- same_number(expected_delta, result_row$log2FC)
    t_match <- same_number(checked_t, result_row$welch_t)
    df_match <- same_number(checked_df, result_row$welch_df)
    p_match <- same_number(checked_p, result_row$p_value)

    data.frame(
      context = context$stem,
      gene = gene,
      status = result_row$p_value_status,
      target_mean_check = mean(target_values),
      control_mean_check = mean(control_values),
      log2fc_check = expected_delta,
      welch_t_check = checked_t,
      welch_df_check = checked_df,
      p_value_check = checked_p,
      check_method = check_method,
      passed = delta_match && t_match && df_match && p_match
    )
  })
  context_rows <- do.call(rbind, context_rows)
  rows[[length(rows) + 1]] <- context_rows

  recalculated_fdr <- p.adjust(results$p_value, method = "BH")
  maximum_fdr_difference <- max(abs(recalculated_fdr - results$fdr_bh))
  summaries[[length(summaries) + 1]] <- data.frame(
    context = context$stem,
    tested_transcripts = nrow(results),
    candidate_checks = sum(context_rows$status != "gene_not_available"),
    candidate_failures = sum(!context_rows$passed),
    maximum_absolute_bh_difference = maximum_fdr_difference,
    bh_recalculation_passed = maximum_fdr_difference <= 1e-8,
    method = "independent R t.test and p.adjust checks"
  )
}

check_table <- do.call(rbind, rows)
summary_table <- do.call(rbind, summaries)
write.csv(
  check_table,
  file.path(results_dir, "external_welch_validation.csv"),
  row.names = FALSE
)
write.csv(
  summary_table,
  file.path(results_dir, "external_welch_validation_summary.csv"),
  row.names = FALSE
)

manifest <- list(
  generated_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  scripts = c(
    "code/R/run_external_welch_bh.R",
    "code/R/validate_external_welch_bh.R"
  ),
  checks = list(
    selected_candidate_rows_recomputed_with_t_test = sum(
      check_table$status != "gene_not_available"
    ),
    selected_candidate_failures = sum(!check_table$passed),
    bh_recalculated_for_all_contexts = all(summary_table$bh_recalculation_passed),
    maximum_absolute_bh_difference = max(
      summary_table$maximum_absolute_bh_difference
    )
  ),
  contexts = summary_table
)
if (requireNamespace("jsonlite", quietly = TRUE)) {
  jsonlite::write_json(
    manifest,
    file.path(results_dir, "external_welch_validation.json"),
    pretty = TRUE,
    auto_unbox = TRUE,
    na = "null"
  )
}

print(summary_table)
if (any(!check_table$passed) || any(!summary_table$bh_recalculation_passed)) {
  stop("External Welch/BH validation failed.")
}
