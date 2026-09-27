#!/usr/bin/env Rscript
# Build the per-age input cache for the full CellChatDB.mouse run.
#
# All 2,019 ligand-receptor pairs of CellChatDB.mouse are evaluated. Database
# genes that are absent from the GSE274279 matrix (e.g. several H2/Ifna genes)
# are added as all-zero rows, because an undetected gene is equivalent to zero
# expression and this keeps every database interaction in the analysis.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript cellchat_full_cache.R <3|12|24>")
}
age <- args[[1]]

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
lib <- "/Users/jijiabiao/Desktop/耳聋课题/Rlibs"
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(CellChat)
})

data_dir <- "/tmp/ear_study_learning/GSE274279"
out_dir <- file.path(root, "cellchat_full", "cache")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

file_map <- list(
  "3" = c("3M_matrix.mtx.gz", "3M_barcodes.tsv.gz", "-0"),
  "12" = c("12M_matrix.mtx.gz", "12M_barcodes.tsv.gz", "-1"),
  "24" = c("24M_matrix.mtx.gz", "24M_barcodes.tsv.gz", "-2")
)
if (!age %in% names(file_map)) {
  stop("age must be one of 3, 12, 24")
}

features <- read.delim(
  gzfile(file.path(data_dir, "features.tsv.gz")),
  header = FALSE,
  stringsAsFactors = FALSE
)
gene_symbols <- features[[2]]

umap <- read.csv(
  file.path(root, "GSE274279_umap_coordinates.csv"),
  row.names = 1,
  check.names = FALSE
)

matrix_name <- file_map[[age]][1]
barcode_name <- file_map[[age]][2]
suffix <- file_map[[age]][3]
barcodes <- readLines(gzfile(file.path(data_dir, barcode_name)))
barcode_index <- paste0(barcodes, suffix)
age_meta <- umap[umap$age_months == as.numeric(age), , drop = FALSE]
keep <- barcode_index %in% rownames(age_meta)
barcodes <- barcodes[keep]
barcode_index <- barcode_index[keep]

counts <- Matrix::readMM(gzfile(file.path(data_dir, matrix_name)))
counts <- as(counts, "CsparseMatrix")
counts <- counts[, keep, drop = FALSE]
rownames(counts) <- gene_symbols
colnames(counts) <- barcode_index
counts <- counts[!duplicated(rownames(counts)), , drop = FALSE]

totals <- Matrix::colSums(counts)
totals[totals <= 0] <- 1
normalized <- t(t(counts) * (1e4 / totals))
normalized@x <- log1p(normalized@x)
rm(counts)
gc(verbose = FALSE)

db <- CellChatDB.mouse
complex_names <- rownames(db$complex)
cofactor_names <- rownames(db$cofactor)
complex_map <- setNames(
  lapply(seq_len(nrow(db$complex)), function(index) {
    subunits <- unlist(db$complex[index, , drop = TRUE], use.names = FALSE)
    unique(subunits[!is.na(subunits) & nzchar(subunits)])
  }),
  complex_names
)
expand_entity <- function(entity) {
  entity <- as.character(entity)
  if (entity %in% names(complex_map)) {
    return(complex_map[[entity]])
  }
  entity
}

# Only genuine gene symbols may be zero-filled. Complex and cofactor names must
# stay out of the expression matrix: computeExpr_LR treats any entity that is a
# row name as a single gene, which would replace the complex geometric mean by
# a constant zero.
entity_genes <- unique(c(
  as.character(db$interaction$ligand),
  as.character(db$interaction$receptor)
))
entity_genes <- unique(unlist(lapply(entity_genes, expand_entity), use.names = FALSE))
complex_genes <- unlist(db$complex, use.names = FALSE)
cofactor_genes <- unlist(db$cofactor, use.names = FALSE)
cofactor_genes <- setdiff(cofactor_genes, c(complex_names, cofactor_names))
db_genes <- unique(c(entity_genes, complex_genes, cofactor_genes))
db_genes <- db_genes[!is.na(db_genes) & nzchar(db_genes)]

present <- intersect(db_genes, rownames(normalized))
missing <- setdiff(db_genes, present)
mat <- normalized[present, , drop = FALSE]
if (length(missing) > 0) {
  zero_block <- Matrix::sparseMatrix(
    i = integer(0),
    j = integer(0),
    dims = c(length(missing), ncol(mat)),
    dimnames = list(missing, colnames(mat))
  )
  mat <- rbind(mat, zero_block)
}
mat <- mat[db_genes, , drop = FALSE]

meta <- data.frame(
  cell_type = age_meta[barcode_index, "cell_type"],
  row.names = barcode_index,
  stringsAsFactors = FALSE
)

cache <- list(
  age_months = as.numeric(age),
  mat = mat,
  meta = meta,
  lr = db$interaction,
  db_genes = db_genes,
  missing_genes = missing,
  n_cells_by_type = table(meta$cell_type)
)
saveRDS(
  cache,
  file.path(out_dir, paste0("cache_", age, "M.rds")),
  compress = "gzip"
)

cat(
  "Cached age", age, "M:",
  format(ncol(mat), big.mark = ","), "cells x",
  format(nrow(mat), big.mark = ","), "database genes;",
  "missing genes zero-filled:", length(missing), ";",
  "interactions:", nrow(db$interaction), "\n"
)
