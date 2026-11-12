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
  library(jsonlite)
  library(Seurat)
})

REPO_ROOT <- v14_repo_root()
INPUT_RDS <- v14_path(
  "V14_MYELOID_RDS",
  file.path(
    REPO_ROOT, "results", "single_cell", "GSE274279_four_matrix",
    "GSE274279_four_matrix_myeloid.rds"
  )
)
OUTPUT_DIR <- v14_path(
  "V14_MYELOID_OUTPUT_DIR",
  dirname(INPUT_RDS)
)
CANDIDATES <- c("Cd68", "Clec7a", "Fcgr3", "Lyz1", "Lyz2", "Mpeg1", "Ptprc")

object <- readRDS(INPUT_RDS)
if (!inherits(object, "Seurat")) {
  stop("Input object is not a Seurat object: ", INPUT_RDS)
}
metadata <- object@meta.data
required <- c("age_months", "myeloid_state")
if (any(!required %in% colnames(metadata))) {
  stop("Required myeloid metadata are missing.")
}

data_layer <- GetAssayData(object, assay = "RNA", layer = "data")
rows <- list()
for (gene in CANDIDATES) {
  if (!gene %in% rownames(data_layer)) {
    next
  }
  for (age in sort(unique(metadata$age_months))) {
    selected <- metadata$age_months == age
    values <- as.numeric(data_layer[gene, selected])
    rows[[paste(gene, age, sep = "__")]] <- data.frame(
      gene = gene,
      cell_type = "Macrophages/Microglia",
      age_months = as.integer(age),
      cells = length(values),
      mean_log_expression = mean(values),
      percent_expressing = mean(values > 0) * 100,
      stringsAsFactors = FALSE
    )
  }
}
expression_table <- do.call(rbind, rows)
rownames(expression_table) <- NULL
write.csv(
  expression_table,
  file.path(OUTPUT_DIR, "GSE274279_myeloid_candidate_expression.csv"),
  row.names = FALSE
)

available <- intersect(CANDIDATES, rownames(data_layer))
missing <- setdiff(CANDIDATES, available)
write_json(
  list(
    source_object = v14_relative_path(INPUT_RDS, REPO_ROOT),
    nuclei = ncol(object),
    genes = nrow(object),
    candidates_requested = CANDIDATES,
    candidates_available = available,
    candidates_unavailable = missing,
    ages_months = sort(unique(as.integer(metadata$age_months))),
    note = sprintf(
      paste(
        "Six primary-intersection genes were measurable in the prepared",
        "%d-nucleus myeloid object. Clec7a was absent. Candidate genes",
        "Cd68, Fcgr3, Mpeg1, and Ptprc were also marker-panel genes;",
        "therefore these detections do not provide independent localization."
      ),
      ncol(object)
    )
  ),
  file.path(OUTPUT_DIR, "GSE274279_myeloid_candidate_expression_summary.json"),
  pretty = TRUE,
  auto_unbox = TRUE
)

cat(
  "Wrote six-gene myeloid expression table and availability summary to",
  OUTPUT_DIR,
  "\n"
)
print(expression_table)
