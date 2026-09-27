#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript cellchat_full_by_age.R <3|12|24>")
}
age <- args[[1]]

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
lib <- "/Users/jijiabiao/Desktop/耳聋课题/Rlibs"
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(CellChat)
  library(future)
})

data_dir <- "/tmp/ear_study_learning/GSE274279"
out <- file.path(root, "cellchat_full")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(2026)
options(future.globals.maxSize = 8 * 1024^3)
future::plan("multisession", workers = 2)

file_map <- list(
  "3" = c("3M_matrix.mtx.gz", "3M_barcodes.tsv.gz", "-0"),
  "12" = c("12M_matrix.mtx.gz", "12M_barcodes.tsv.gz", "-1"),
  "24" = c("24M_matrix.mtx.gz", "24M_barcodes.tsv.gz", "-2")
)

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

meta <- data.frame(
  cell_type = age_meta[barcode_index, "cell_type"],
  row.names = barcode_index,
  stringsAsFactors = FALSE
)

complex_db <- CellChatDB.mouse$complex
complex_map <- setNames(
  lapply(seq_len(nrow(complex_db)), function(index) {
    subunits <- unlist(complex_db[index, , drop = TRUE], use.names = FALSE)
    unique(subunits[!is.na(subunits) & nzchar(subunits)])
  }),
  rownames(complex_db)
)
expand_entity <- function(entity) {
  entity <- as.character(entity)
  if (entity %in% names(complex_map)) {
    return(complex_map[[entity]])
  }
  entity
}

available_genes <- rownames(normalized)
interaction_db <- CellChatDB.mouse$interaction
keep_interaction <- vapply(
  seq_len(nrow(interaction_db)),
  function(index) {
    members <- c(
      expand_entity(interaction_db$ligand[index]),
      expand_entity(interaction_db$receptor[index])
    )
    all(members %in% available_genes)
  },
  logical(1)
)
interaction_db <- interaction_db[keep_interaction, , drop = FALSE]

cellchat <- createCellChat(
  object = normalized,
  meta = meta,
  group.by = "cell_type"
)
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat@LR$LRsig <- interaction_db
cat("Age", age, "M full CellChat interactions:", nrow(interaction_db), "\n")

cellchat <- computeCommunProb(cellchat, type = "triMean", nboot = 20)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

probabilities <- cellchat@net$prob
probability_index <- which(probabilities > 0, arr.ind = TRUE)
probability_frame <- data.frame(
  age_months = as.numeric(age),
  source = dimnames(probabilities)[[1]][probability_index[, 1]],
  target = dimnames(probabilities)[[2]][probability_index[, 2]],
  interaction_name = dimnames(probabilities)[[3]][probability_index[, 3]],
  prob = probabilities[probability_index],
  pval = cellchat@net$pval[probability_index],
  stringsAsFactors = FALSE
)

pathway_probabilities <- cellchat@netP$prob
if (length(dim(pathway_probabilities)) == 3 &&
    !is.null(dimnames(pathway_probabilities)[[3]])) {
  pathway_index <- which(pathway_probabilities > 0, arr.ind = TRUE)
  pathway_frame <- data.frame(
    age_months = as.numeric(age),
    source = dimnames(pathway_probabilities)[[1]][pathway_index[, 1]],
    target = dimnames(pathway_probabilities)[[2]][pathway_index[, 2]],
    pathway_name = dimnames(pathway_probabilities)[[3]][pathway_index[, 3]],
    prob = pathway_probabilities[pathway_index],
    pval = cellchat@netP$pval[pathway_index],
    stringsAsFactors = FALSE
  )
} else {
  pathway_frame <- data.frame()
}

saveRDS(cellchat, file.path(out, paste0("CellChat_full_GSE274279_", age, "M.rds")))
write.csv(
  probability_frame,
  file.path(out, paste0("CellChat_full_network_", age, "M.csv")),
  row.names = FALSE
)
write.csv(
  pathway_frame,
  file.path(out, paste0("CellChat_full_pathway_", age, "M.csv")),
  row.names = FALSE
)
future::plan("sequential")
cat("Full CellChat complete for", age, "M\n")
