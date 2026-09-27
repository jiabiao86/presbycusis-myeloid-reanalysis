#!/usr/bin/env Rscript
# Merge the chunked full CellChatDB.mouse run of one age into the final object.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript cellchat_full_merge.R <3|12|24> [nboot]")
}
age <- args[[1]]
nboot <- if (length(args) >= 2) as.integer(args[[2]]) else 100L

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
lib <- "/Users/jijiabiao/Desktop/耳聋课题/Rlibs"
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(CellChat)
  library(future)
})
future::plan("sequential")

chunk_dir <- file.path(root, "cellchat_full", "chunks")
final_dir <- file.path(root, "cellchat_full", "final")
dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)

cache <- readRDS(
  file.path(root, "cellchat_full", "cache", paste0("cache_", age, "M.rds"))
)
lr <- cache$lr
n_lr <- nrow(lr)

cellchat <- createCellChat(
  object = cache$mat,
  meta = cache$meta,
  group.by = "cell_type"
)
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
restored <- setdiff(rownames(cellchat@data), rownames(cellchat@data.signaling))
if (length(restored) > 0) {
  cellchat@data.signaling <- rbind(
    cellchat@data.signaling,
    cellchat@data[restored, , drop = FALSE]
  )
}
cellchat@LR$LRsig <- lr

groups <- levels(cellchat@idents)
n_group <- length(groups)
prob <- array(
  0,
  dim = c(n_group, n_group, n_lr),
  dimnames = list(groups, groups, rownames(lr))
)
pval <- array(
  1,
  dim = c(n_group, n_group, n_lr),
  dimnames = list(groups, groups, rownames(lr))
)

chunk_files <- list.files(
  chunk_dir,
  pattern = paste0("^chunk_", age, "M_"),
  full.names = TRUE
)
if (length(chunk_files) == 0) {
  stop("no chunk files for age ", age)
}

covered <- logical(n_lr)
elapsed_total <- 0
for (file in chunk_files) {
  chunk <- readRDS(file)
  if (chunk$nboot != nboot) {
    stop(file, " used nboot = ", chunk$nboot, ", expected ", nboot)
  }
  prob[, , chunk$start:chunk$end] <- chunk$prob
  pval[, , chunk$start:chunk$end] <- chunk$pval
  covered[chunk$start:chunk$end] <- TRUE
  elapsed_total <- elapsed_total + chunk$elapsed_sec
}
if (!all(covered)) {
  missing_ranges <- which(!covered)
  stop("interactions without chunk results: ", length(missing_ranges))
}

cellchat@net <- list(prob = prob, pval = pval)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat@options$parameter <- list(
  type.mean = "triMean",
  nboot = nboot,
  seed.use = 1L,
  database = "CellChatDB.mouse",
  n_interactions = n_lr,
  chunked = TRUE
)

saveRDS(
  cellchat,
  file.path(final_dir, paste0("CellChat_full_GSE274279_", age, "M.rds"))
)

network <- data.frame(
  age_months = as.numeric(age),
  source = rep(groups, times = n_group * n_lr),
  target = rep(groups, each = n_group, times = n_lr),
  interaction_name = rep(dimnames(prob)[[3]], each = n_group * n_group),
  prob = as.vector(prob),
  pval = as.vector(pval),
  stringsAsFactors = FALSE
)
network <- network[network$prob > 0, , drop = FALSE]
write.csv(
  network,
  file.path(final_dir, paste0("CellChat_full_network_", age, "M.csv")),
  row.names = FALSE
)
write.csv(
  network[network$pval < 0.05, , drop = FALSE],
  file.path(final_dir, paste0("CellChat_full_network_significant_", age, "M.csv")),
  row.names = FALSE
)

pathway_prob <- cellchat@netP$prob
pathway_names <- dimnames(pathway_prob)[[3]]
pathway_frame <- data.frame(
  age_months = as.numeric(age),
  source = rep(groups, times = n_group * length(pathway_names)),
  target = rep(groups, each = n_group, times = length(pathway_names)),
  pathway_name = rep(pathway_names, each = n_group * n_group),
  prob = as.vector(pathway_prob),
  stringsAsFactors = FALSE
)
pathway_frame <- pathway_frame[pathway_frame$prob > 0, , drop = FALSE]
write.csv(
  pathway_frame,
  file.path(final_dir, paste0("CellChat_full_pathway_", age, "M.csv")),
  row.names = FALSE
)

summary_frame <- data.frame(
  age_months = as.numeric(age),
  n_interactions_evaluated = n_lr,
  n_nonzero_edges = sum(prob > 0),
  n_significant_edges = sum(prob > 0 & pval < 0.05),
  n_significant_interactions = length(unique(
    network$interaction_name[network$pval < 0.05]
  )),
  n_significant_pathways = length(pathway_names),
  nboot = nboot,
  seed_use = 1L,
  chunk_elapsed_sec = elapsed_total
)
write.csv(
  summary_frame,
  file.path(final_dir, paste0("CellChat_full_summary_", age, "M.csv")),
  row.names = FALSE
)

print(summary_frame)
cat("Full CellChatDB.mouse run merged for", age, "M\n")
