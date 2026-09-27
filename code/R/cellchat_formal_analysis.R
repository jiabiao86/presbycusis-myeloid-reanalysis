#!/usr/bin/env Rscript

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
lib <- "/Users/jijiabiao/Desktop/耳聋课题/Rlibs"
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages({
  library(Matrix)
  library(CellChat)
  library(future)
})

data_dir <- "/tmp/ear_study_learning/GSE274279"
out <- file.path(root, "cellchat_formal")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

set.seed(2026)
options(future.globals.maxSize = 8 * 1024^3)
future::plan("multisession", workers = 8)

age_files <- list(
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

module_assignment <- read.csv(
  file.path(root, "wgcna_module_assignment.csv"),
  check.names = FALSE
)
severity_trend <- read.csv(
  file.path(root, "limma_severity_trend.csv"),
  check.names = FALSE
)
target_genes <- unique(c(
  module_assignment$gene[module_assignment$module == "M12"],
  severity_trend$gene[
    severity_trend$adj.P.Val < 0.05 &
      severity_trend$logFC > 0
  ]
))
target_genes <- target_genes[!is.na(target_genes)]

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

interaction_db <- CellChatDB.mouse$interaction
keep_interaction <- vapply(
  seq_len(nrow(interaction_db)),
  function(index) {
    members <- c(
      expand_entity(interaction_db$ligand[index]),
      expand_entity(interaction_db$receptor[index])
    )
    any(members %in% target_genes)
  },
  logical(1)
)
interaction_db <- interaction_db[keep_interaction, , drop = FALSE]
CellChatDB.focused <- CellChatDB.mouse
CellChatDB.focused$interaction <- interaction_db

cat(
  "CellChatDB interactions retained for severity-associated ligands/receptors:",
  nrow(interaction_db),
  "of",
  nrow(CellChatDB.mouse$interaction),
  "\n"
)

cellchat_list <- list()
communication_list <- list()
pathway_list <- list()
network_probability_list <- list()
pathway_probability_list <- list()

for (age in names(age_files)) {
  matrix_name <- age_files[[age]][1]
  barcode_name <- age_files[[age]][2]
  suffix <- age_files[[age]][3]
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

  cellchat <- createCellChat(
    object = normalized,
    meta = meta,
    group.by = "cell_type"
  )
  cellchat@DB <- CellChatDB.focused
  cellchat <- subsetData(cellchat)
  cellchat <- identifyOverExpressedGenes(cellchat)
  available_genes <- rownames(normalized)
  lr_available <- interaction_db[vapply(
    seq_len(nrow(interaction_db)),
    function(index) {
      members <- c(
        expand_entity(interaction_db$ligand[index]),
        expand_entity(interaction_db$receptor[index])
      )
      all(members %in% available_genes)
    },
    logical(1)
  ), , drop = FALSE]
  cellchat@LR$LRsig <- lr_available
  cat("Available CellChat interactions for", age, "M:", nrow(lr_available), "\n")
  cellchat <- computeCommunProb(cellchat, type = "triMean", nboot = 100)
  cellchat <- filterCommunication(cellchat, min.cells = 10)
  cellchat <- computeCommunProbPathway(cellchat)
  cellchat <- aggregateNet(cellchat)

  communication <- subsetCommunication(cellchat)
  communication$age_months <- as.numeric(age)
  pathway <- subsetCommunication(cellchat, slot.name = "netP")
  pathway$age_months <- as.numeric(age)

  network_probability <- cellchat@net$prob
  network_index <- which(network_probability > 0, arr.ind = TRUE)
  network_probability_frame <- data.frame(
    source = dimnames(network_probability)[[1]][network_index[, 1]],
    target = dimnames(network_probability)[[2]][network_index[, 2]],
    interaction_name = dimnames(network_probability)[[3]][network_index[, 3]],
    prob = network_probability[network_index],
    age_months = as.numeric(age),
    stringsAsFactors = FALSE
  )

  pathway_probability <- cellchat@netP$prob
  if (length(dim(pathway_probability)) == 3 &&
      !is.null(dimnames(pathway_probability)[[3]])) {
    pathway_index <- which(pathway_probability > 0, arr.ind = TRUE)
    pathway_probability_frame <- data.frame(
      source = dimnames(pathway_probability)[[1]][pathway_index[, 1]],
      target = dimnames(pathway_probability)[[2]][pathway_index[, 2]],
      pathway_name = dimnames(pathway_probability)[[3]][pathway_index[, 3]],
      prob = pathway_probability[pathway_index],
      age_months = as.numeric(age),
      stringsAsFactors = FALSE
    )
  } else {
    pathway_probability_frame <- data.frame(
      source = character(),
      target = character(),
      pathway_name = character(),
      prob = numeric(),
      age_months = numeric()
    )
  }

  cellchat_list[[age]] <- cellchat
  communication_list[[age]] <- communication
  pathway_list[[age]] <- pathway
  network_probability_list[[age]] <- network_probability_frame
  pathway_probability_list[[age]] <- pathway_probability_frame

  saveRDS(
    cellchat,
    file.path(out, paste0("CellChat_GSE274279_", age, "M.rds"))
  )
  write.csv(
    communication,
    file.path(out, paste0("CellChat_communication_", age, "M.csv")),
    row.names = FALSE
  )
  write.csv(
    pathway,
    file.path(out, paste0("CellChat_pathway_", age, "M.csv")),
    row.names = FALSE
  )
  write.csv(
    network_probability_frame,
    file.path(out, paste0("CellChat_network_all_probabilities_", age, "M.csv")),
    row.names = FALSE
  )
  write.csv(
    pathway_probability_frame,
    file.path(out, paste0("CellChat_pathway_all_probabilities_", age, "M.csv")),
    row.names = FALSE
  )
  cat("Completed CellChat for", age, "M:", nrow(communication), "communications\n")
}

all_communication <- do.call(rbind, communication_list)
all_pathway <- do.call(rbind, pathway_list)
all_network_probability <- do.call(rbind, network_probability_list)
all_pathway_probability <- do.call(rbind, pathway_probability_list)
write.csv(
  all_communication,
  file.path(out, "CellChat_all_communications.csv"),
  row.names = FALSE
)
write.csv(
  all_pathway,
  file.path(out, "CellChat_all_pathways.csv"),
  row.names = FALSE
)
write.csv(
  all_network_probability,
  file.path(out, "CellChat_network_all_probabilities.csv"),
  row.names = FALSE
)
write.csv(
  all_pathway_probability,
  file.path(out, "CellChat_pathway_all_probabilities.csv"),
  row.names = FALSE
)

communication_key <- c(
  "source",
  "target",
  "ligand",
  "receptor",
  "interaction_name",
  "pathway_name",
  "annotation"
)
communication_grouped <- aggregate(
  prob ~ age_months + source + target + ligand + receptor + interaction_name + pathway_name + annotation,
  data = all_communication,
  FUN = mean
)
comm_3 <- communication_grouped[communication_grouped$age_months == 3, ]
comm_24 <- communication_grouped[communication_grouped$age_months == 24, ]
if (nrow(comm_3) > 0 || nrow(comm_24) > 0) {
  comm_changes <- merge(
    comm_3,
    comm_24,
    by = communication_key,
    all = TRUE,
    suffixes = c("_3m", "_24m")
  )
  comm_changes$prob_3m[is.na(comm_changes$prob_3m)] <- 0
  comm_changes$prob_24m[is.na(comm_changes$prob_24m)] <- 0
  comm_changes$delta_24m_minus_3m <- comm_changes$prob_24m - comm_changes$prob_3m
  comm_changes$log2FC_24m_vs_3m <- log2(
    (comm_changes$prob_24m + 1e-6) / (comm_changes$prob_3m + 1e-6)
  )
  comm_changes$involving_myeloid <- comm_changes$source == "Macrophages/Microglia" |
    comm_changes$target == "Macrophages/Microglia"
  comm_changes <- comm_changes[
    order(comm_changes$delta_24m_minus_3m, decreasing = TRUE),
  ]
} else {
  comm_changes <- data.frame()
}
write.csv(
  comm_changes,
  file.path(out, "CellChat_communication_changes_24M_vs_3M.csv"),
  row.names = FALSE
)
write.csv(
  head(comm_changes[comm_changes$involving_myeloid, ], 100),
  file.path(out, "CellChat_top_myeloid_changes.csv"),
  row.names = FALSE
)

network_3 <- all_network_probability[all_network_probability$age_months == 3, ]
network_24 <- all_network_probability[all_network_probability$age_months == 24, ]
network_changes <- merge(
  network_3,
  network_24,
  by = c("source", "target", "interaction_name"),
  all = TRUE,
  suffixes = c("_3m", "_24m")
)
network_changes$prob_3m[is.na(network_changes$prob_3m)] <- 0
network_changes$prob_24m[is.na(network_changes$prob_24m)] <- 0
network_changes$delta_24m_minus_3m <- network_changes$prob_24m - network_changes$prob_3m
network_changes$involving_myeloid <- network_changes$source == "Macrophages/Microglia" |
  network_changes$target == "Macrophages/Microglia"
network_changes <- network_changes[
  order(network_changes$delta_24m_minus_3m, decreasing = TRUE),
]
write.csv(
  network_changes,
  file.path(out, "CellChat_network_probability_changes_24M_vs_3M.csv"),
  row.names = FALSE
)

pathway_probability_3 <- all_pathway_probability[all_pathway_probability$age_months == 3, ]
pathway_probability_24 <- all_pathway_probability[all_pathway_probability$age_months == 24, ]
pathway_probability_changes <- merge(
  pathway_probability_3,
  pathway_probability_24,
  by = c("source", "target", "pathway_name"),
  all = TRUE,
  suffixes = c("_3m", "_24m")
)
pathway_probability_changes$prob_3m[is.na(pathway_probability_changes$prob_3m)] <- 0
pathway_probability_changes$prob_24m[is.na(pathway_probability_changes$prob_24m)] <- 0
pathway_probability_changes$delta_24m_minus_3m <-
  pathway_probability_changes$prob_24m - pathway_probability_changes$prob_3m
pathway_probability_changes$involving_myeloid <-
  pathway_probability_changes$source == "Macrophages/Microglia" |
    pathway_probability_changes$target == "Macrophages/Microglia"
pathway_probability_changes <- pathway_probability_changes[
  order(pathway_probability_changes$delta_24m_minus_3m, decreasing = TRUE),
]
write.csv(
  pathway_probability_changes,
  file.path(out, "CellChat_pathway_probability_changes_24M_vs_3M.csv"),
  row.names = FALSE
)

pathway_grouped <- aggregate(
  prob ~ age_months + source + target + pathway_name,
  data = all_pathway,
  FUN = sum
)
pathway_3 <- pathway_grouped[pathway_grouped$age_months == 3, ]
pathway_24 <- pathway_grouped[pathway_grouped$age_months == 24, ]
if (nrow(pathway_3) > 0 || nrow(pathway_24) > 0) {
  pathway_changes <- merge(
    pathway_3,
    pathway_24,
    by = c("source", "target", "pathway_name"),
    all = TRUE,
    suffixes = c("_3m", "_24m")
  )
  pathway_changes$prob_3m[is.na(pathway_changes$prob_3m)] <- 0
  pathway_changes$prob_24m[is.na(pathway_changes$prob_24m)] <- 0
  pathway_changes$delta_24m_minus_3m <- pathway_changes$prob_24m - pathway_changes$prob_3m
  pathway_changes$involving_myeloid <- pathway_changes$source == "Macrophages/Microglia" |
    pathway_changes$target == "Macrophages/Microglia"
  pathway_changes <- pathway_changes[
    order(pathway_changes$delta_24m_minus_3m, decreasing = TRUE),
  ]
} else {
  pathway_changes <- data.frame()
}
write.csv(
  pathway_changes,
  file.path(out, "CellChat_pathway_changes_24M_vs_3M.csv"),
  row.names = FALSE
)

future::plan("sequential")
cat("Formal CellChat analysis complete:", out, "\n")
