#!/usr/bin/env Rscript
# Run one chunk of the full CellChatDB.mouse interaction set for one age.
#
# The permutation seed is fixed per (age, nboot), and CellChat draws the same
# permutation matrix for every chunk because it depends only on the number of
# cells and nboot. Chunk results are therefore identical to a single-pass run
# and can be merged into one network.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) {
  stop("Usage: Rscript cellchat_full_chunk.R <3|12|24> <start> <end> [nboot]")
}
age <- args[[1]]
start <- as.integer(args[[2]])
end <- as.integer(args[[3]])
nboot <- if (length(args) >= 4) as.integer(args[[4]]) else 100L

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
lib <- "/Users/jijiabiao/Desktop/耳聋课题/Rlibs"
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(CellChat)
  library(future)
})

# Sequential plan keeps the inner bootstrap loop on plain sapply; a multisession
# plan makes every coreceptor call dispatch a future job and dominates runtime.
future::plan("sequential")
options(future.globals.maxSize = 8 * 1024^3)

chunk_dir <- file.path(root, "cellchat_full", "chunks")
dir.create(chunk_dir, recursive = TRUE, showWarnings = FALSE)

cache <- readRDS(
  file.path(root, "cellchat_full", "cache", paste0("cache_", age, "M.rds"))
)
lr <- cache$lr
if (start < 1 || end > nrow(lr) || start > end) {
  stop("requested range outside 1..", nrow(lr))
}
chunk_lr <- lr[start:end, , drop = FALSE]

cellchat <- createCellChat(
  object = cache$mat,
  meta = cache$meta,
  group.by = "cell_type"
)
cellchat@DB <- CellChatDB.mouse
cellchat <- subsetData(cellchat)
# extractGene() drops database genes whose symbols are absent from CellChat's
# official mouse symbol table (here H2-BI and H2-Ea-ps). They are zero-filled
# rows in the cache and must stay in data.signaling: otherwise computeExpr_LR
# treats them as complex names and computeExpr_complex fails on the unknown row.
restored <- setdiff(rownames(cellchat@data), rownames(cellchat@data.signaling))
if (length(restored) > 0) {
  cellchat@data.signaling <- rbind(
    cellchat@data.signaling,
    cellchat@data[restored, , drop = FALSE]
  )
  cat("Restored zero-filled database genes:", paste(restored, collapse = ", "), "\n")
}
cellchat@LR$LRsig <- chunk_lr

cat(
  "Age", age, "M chunk", start, "-", end, ":", nrow(chunk_lr),
  "interactions, nboot =", nboot, "[", format(Sys.time()), "]\n"
)

t0 <- Sys.time()
cellchat <- computeCommunProb(
  cellchat,
  type = "triMean",
  nboot = nboot,
  seed.use = 1L
)
elapsed <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

result <- list(
  age_months = as.numeric(age),
  start = start,
  end = end,
  nboot = nboot,
  seed_use = 1L,
  prob = cellchat@net$prob,
  pval = cellchat@net$pval,
  lr = chunk_lr,
  elapsed_sec = elapsed
)
saveRDS(
  result,
  file.path(chunk_dir, paste0("chunk_", age, "M_", start, "-", end, ".rds"))
)
cat(
  "Chunk done in", round(elapsed, 1), "s;",
  "bootstrapped interactions:", sum(apply(result$prob, 3, sum) > 0),
  "[", format(Sys.time()), "]\n"
)
