#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run with Rscript so repository paths can be resolved.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
source(file.path(dirname(script_file), "common_v14.R"))
REPO_ROOT <- v14_repo_root()
single_cell_output <- v14_path(
  "V14_SINGLE_CELL_OUTPUT_DIR",
  file.path(
    REPO_ROOT, "results", "single_cell", "GSE274279_four_matrix"
  )
)
results <- v14_path(
  "V14_MYELOID_AUDIT_OUTPUT_DIR",
  single_cell_output
)
dir.create(results, recursive = TRUE, showWarnings = FALSE)

metadata_path <- file.path(
  single_cell_output,
  "myeloid_four_matrix_cell_metadata.csv"
)
if (!file.exists(metadata_path)) {
  stop("Missing myeloid cell metadata: ", metadata_path)
}

metadata <- read.csv(metadata_path, row.names = 1, check.names = FALSE)
metadata$cell_id <- rownames(metadata)
required <- c("cell_id", "age_months", "myeloid_state")
if (any(!required %in% colnames(metadata))) {
  stop("Myeloid metadata lacks required columns: ",
       paste(setdiff(required, colnames(metadata)), collapse = ", "))
}
if (anyNA(metadata[, required])) {
  stop("Required myeloid metadata fields contain NA.")
}
if (anyDuplicated(metadata$cell_id)) {
  stop("Duplicate myeloid cell IDs are present.")
}
if (any(!metadata$age_months %in% c(3, 12, 24))) {
  stop("Unexpected age value in myeloid metadata.")
}

audit <- data.frame(
  cell_id = metadata$cell_id,
  age_months = metadata$age_months,
  myeloid_state = metadata$myeloid_state,
  included_in_245_total = TRUE,
  exclusion_reason = NA_character_,
  source_file = file.path(
    basename(single_cell_output),
    basename(metadata_path)
  ),
  stringsAsFactors = FALSE
)
audit <- audit[order(audit$age_months, audit$myeloid_state, audit$cell_id), ]
write.csv(
  audit,
  file.path(results, "myeloid_cell_id_audit.csv"),
  row.names = FALSE,
  na = "NA"
)

counts <- aggregate(
  cell_id ~ age_months + myeloid_state,
  data = audit,
  FUN = length
)
names(counts)[names(counts) == "cell_id"] <- "n_cells"
counts <- counts[order(counts$age_months, counts$myeloid_state), ]
counts$proportion_within_age <- ave(
  counts$n_cells,
  counts$age_months,
  FUN = function(x) x / sum(x)
)
write.csv(
  counts,
  file.path(results, "myeloid_count_audit.csv"),
  row.names = FALSE,
  na = "NA"
)

totals <- aggregate(
  cell_id ~ age_months,
  data = audit,
  FUN = length
)
names(totals)[names(totals) == "cell_id"] <- "n_cells"
EXPECTED_MYELOID_NUCLEI <- 245L
if (sum(totals$n_cells) != EXPECTED_MYELOID_NUCLEI) {
  stop(
    "Myeloid cell count audit does not total ",
    EXPECTED_MYELOID_NUCLEI,
    "."
  )
}
write.csv(
  totals,
  file.path(results, "myeloid_count_totals.csv"),
  row.names = FALSE
)

cat(
  "Myeloid cell-ID audit complete: ",
  sum(totals$n_cells),
  " unique nuclei.\n",
  sep = ""
)
print(totals, row.names = FALSE)
