#!/usr/bin/env Rscript
# Deep interpretation of the completed full CellChatDB.mouse run.
#
# Produces: cross-age significant interactions and pathways, cell-type
# communication strength, pathway-level decline, and the status of the
# interactions highlighted in the manuscript.

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
lib <- "/Users/jijiabiao/Desktop/耳聋课题/Rlibs"
.libPaths(c(lib, .libPaths()))

suppressPackageStartupMessages(library(CellChat))

final_dir <- file.path(root, "cellchat_full", "final")
out_dir <- file.path(root, "cellchat_full", "analysis")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ages <- c(3, 12, 24)
objects <- setNames(
  lapply(ages, function(age) {
    readRDS(file.path(final_dir, paste0("CellChat_full_GSE274279_", age, "M.rds")))
  }),
  as.character(ages)
)

edge_tables <- lapply(ages, function(age) {
  object <- objects[[as.character(age)]]
  prob <- object@net$prob
  pval <- object@net$pval
  index <- which(prob > 0, arr.ind = TRUE)
  data.frame(
    age_months = age,
    source = dimnames(prob)[[1]][index[, 1]],
    target = dimnames(prob)[[2]][index[, 2]],
    interaction = dimnames(prob)[[3]][index[, 3]],
    prob = prob[index],
    pval = pval[index],
    significant = pval[index] < 0.05,
    stringsAsFactors = FALSE
  )
})
edges <- do.call(rbind, edge_tables)
write.csv(edges, file.path(out_dir, "CellChat_full_edges_all_ages.csv"), row.names = FALSE)

interaction_presence <- do.call(rbind, lapply(ages, function(age) {
  object <- objects[[as.character(age)]]
  prob <- object@net$prob
  pval <- object@net$pval
  significant <- apply(prob > 0 & pval < 0.05, 3, any)
  detected <- apply(prob > 0, 3, any)
  data.frame(
    interaction = dimnames(prob)[[3]],
    age_months = age,
    detected = detected,
    significant = significant,
    max_prob = apply(prob, 3, max),
    stringsAsFactors = FALSE
  )
}))
presence_matrix <- reshape(
  interaction_presence[, c("interaction", "age_months", "significant")],
  idvar = "interaction",
  timevar = "age_months",
  direction = "wide"
)
names(presence_matrix) <- sub("significant\\.", "sig_", names(presence_matrix))
for (column in c("sig_3", "sig_12", "sig_24")) {
  if (!column %in% names(presence_matrix)) presence_matrix[[column]] <- FALSE
}
presence_matrix$pattern <- paste0(
  ifelse(presence_matrix$sig_3, "3", "-"),
  ifelse(presence_matrix$sig_12, "12", "-"),
  ifelse(presence_matrix$sig_24, "24", "-")
)
presence_matrix$n_ages <- rowSums(presence_matrix[, c("sig_3", "sig_12", "sig_24")])
write.csv(
  presence_matrix,
  file.path(out_dir, "CellChat_full_interaction_presence.csv"),
  row.names = FALSE
)

pathway_tables <- lapply(ages, function(age) {
  object <- objects[[as.character(age)]]
  prob <- object@netP$prob
  index <- which(prob > 0, arr.ind = TRUE)
  data.frame(
    age_months = age,
    source = dimnames(prob)[[1]][index[, 1]],
    target = dimnames(prob)[[2]][index[, 2]],
    pathway = dimnames(prob)[[3]][index[, 3]],
    prob = prob[index],
    stringsAsFactors = FALSE
  )
})
pathways <- do.call(rbind, pathway_tables)
write.csv(pathways, file.path(out_dir, "CellChat_full_pathway_edges_all_ages.csv"), row.names = FALSE)

pathway_strength <- do.call(rbind, lapply(ages, function(age) {
  frame <- pathways[pathways$age_months == age, , drop = FALSE]
  total <- aggregate(prob ~ pathway, frame, sum)
  edges_n <- aggregate(prob ~ pathway, frame, length)
  data.frame(
    age_months = age,
    pathway = total$pathway,
    total_probability = total$prob,
    n_edges = edges_n$prob,
    stringsAsFactors = FALSE
  )
}))
write.csv(
  pathway_strength,
  file.path(out_dir, "CellChat_full_pathway_strength.csv"),
  row.names = FALSE
)

celltype_strength <- do.call(rbind, lapply(ages, function(age) {
  frame <- edges[edges$age_months == age & edges$significant, , drop = FALSE]
  outgoing <- aggregate(prob ~ source, frame, sum)
  incoming <- aggregate(prob ~ target, frame, sum)
  cell_types <- sort(unique(c(outgoing$source, incoming$target)))
  data.frame(
    age_months = age,
    cell_type = cell_types,
    outgoing_probability = outgoing$prob[match(cell_types, outgoing$source)],
    incoming_probability = incoming$prob[match(cell_types, incoming$target)],
    stringsAsFactors = FALSE
  )
}))
celltype_strength[is.na(celltype_strength)] <- 0
write.csv(
  celltype_strength,
  file.path(out_dir, "CellChat_full_celltype_strength.csv"),
  row.names = FALSE
)

focus_interactions <- c(
  "APP_CD74", "PTPRC_MRC1", "PTN_PTPRZ1", "PTN_ALK", "PTN_SDC2",
  "SEMA3A_NRP1_PLXNA2", "SEMA3A_NRP1_PLXNA1", "NRG3_ERBB4", "NRG1_ERBB4",
  "NRXN1_NLGN1", "NRXN3_NLGN1", "NRXN1_NLGN2",
  "LAMA2_ITGA6_ITGB1", "LAMA4_ITGA6_ITGB1", "LAMB1_ITGA6_ITGB1", "COL4A3_ITGA3_ITGB1",
  "FN1_ITGAV_ITGB1", "SPP1_ITGAV_ITGB1", "SPP1_ITGA9_ITGB1",
  "MPZ_MPZ", "MPZ_MPZL1", "CADM1_CADM1", "CDH2_CDH2", "NCAM1_NCAM1",
  "GAS6_MERTK", "PROS1_AXL", "GAS6_AXL", "VEGFA_VEGFR1", "IGF1_IGF1R",
  "TGFB2_TGFBR1_TGFBR2", "BMP5_ACVR1_ACVR2A", "FGF1_FGFR1", "FGF1_FGFR2",
  "EFNA5_EPHA3", "EFNA5_EPHA4", "EFNB2_EPHA4", "NTF3_NTRK2", "NTF3_NTRK3",
  "SPP1_ITGAV_ITGB5", "AGRN_DAG1", "LAMA5_SV2A", "GAS6_AXL"
)
focus <- edges[edges$interaction %in% focus_interactions, , drop = FALSE]
focus <- focus[order(focus$interaction, focus$age_months, -focus$prob), , drop = FALSE]
write.csv(
  focus,
  file.path(out_dir, "CellChat_full_focus_interactions.csv"),
  row.names = FALSE
)

age_summary <- do.call(rbind, lapply(ages, function(age) {
  frame <- edges[edges$age_months == age, , drop = FALSE]
  sig <- frame[frame$significant, , drop = FALSE]
  data.frame(
    age_months = age,
    nonzero_edges = nrow(frame),
    significant_edges = nrow(sig),
    significant_interactions = length(unique(sig$interaction)),
    significant_pathways = length(unique(
      pathways$pathway[pathways$age_months == age]
    )),
    total_significant_probability = sum(sig$prob),
    mean_probability_per_significant_edge = mean(sig$prob),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  age_summary,
  file.path(out_dir, "CellChat_full_age_summary.csv"),
  row.names = FALSE
)

report_path <- file.path(out_dir, "CellChat_full_interpretation.md")
con <- file(report_path, "w")
writeLines("# Full CellChatDB.mouse interpretation (GSE274279, 2,019 interactions)", con)
writeLines("", con)
writeLines("## Communication strength by age", con)
writeLines("", con)
writeLines("| age | nonzero edges | significant edges | significant interactions | significant pathways | total significant probability |", con)
writeLines("| --- | --- | --- | --- | --- | --- |", con)
for (index in seq_len(nrow(age_summary))) {
  row <- age_summary[index, ]
  writeLines(
    paste0(
      "| ", row$age_months, "M | ", row$nonzero_edges, " | ", row$significant_edges,
      " | ", row$significant_interactions, " | ", row$significant_pathways,
      " | ", signif(row$total_significant_probability, 4), " |"
    ),
    con
  )
}
writeLines("", con)

writeLines("## Persistence of significant interactions", con)
writeLines("", con)
writeLines(paste0("- significant in all three ages: ", sum(presence_matrix$n_ages == 3)), con)
writeLines(paste0("- significant in 3M only: ", sum(presence_matrix$pattern == "3---")), con)
writeLines(paste0("- significant in 24M only: ", sum(presence_matrix$pattern == "--24")), con)
writeLines(paste0("- significant in 3M and 12M: ", sum(presence_matrix$pattern == "312---")), con)
writeLines("", con)
core <- presence_matrix[presence_matrix$n_ages == 3, "interaction"]
writeLines(paste0("Core interactions significant at every age (", length(core), "):"), con)
writeLines(paste0("  ", paste(sort(core), collapse = ", ")), con)
writeLines("", con)

writeLines("## Pathways with the largest total probability per age", con)
writeLines("", con)
for (age in ages) {
  frame <- pathway_strength[pathway_strength$age_months == age, , drop = FALSE]
  frame <- frame[order(-frame$total_probability), , drop = FALSE]
  writeLines(paste0("### ", age, "M"), con)
  for (index in seq_len(min(8, nrow(frame)))) {
    row <- frame[index, ]
    writeLines(
      paste0("- ", row$pathway, ": ", signif(row$total_probability, 4),
             " (", row$n_edges, " edges)"),
      con
    )
  }
  writeLines("", con)
}

writeLines("## Communication strength by cell type (significant edges)", con)
writeLines("", con)
for (age in ages) {
  frame <- celltype_strength[celltype_strength$age_months == age, , drop = FALSE]
  frame <- frame[order(-frame$outgoing_probability), , drop = FALSE]
  writeLines(paste0("### ", age, "M"), con)
  for (index in seq_len(nrow(frame))) {
    row <- frame[index, ]
    writeLines(
      paste0("- ", row$cell_type, ": outgoing ", signif(row$outgoing_probability, 4),
             " | incoming ", signif(row$incoming_probability, 4)),
      con
    )
  }
  writeLines("", con)
}

writeLines("## Focus interactions (manuscript candidates)", con)
writeLines("", con)
for (name in sort(unique(focus$interaction))) {
  frame <- focus[focus$interaction == name, , drop = FALSE]
  frame <- frame[order(frame$age_months), , drop = FALSE]
  entries <- vapply(seq_len(nrow(frame)), function(index) {
    row <- frame[index, ]
    paste0(
      row$age_months, "M ", ifelse(row$significant, "sig", "ns"),
      " (prob ", signif(row$prob, 3), ", p ", signif(row$pval, 3), ")"
    )
  }, character(1))
  writeLines(paste0("- ", name, ": ", paste(entries, collapse = "; ")), con)
}
close(con)

cat("interpretation written to", report_path, "\n")
print(age_summary)
cat("\ninteractions significant at all ages:", length(core), "\n")
cat("young-only (3M):", sum(presence_matrix$pattern == "3---"), "\n")
cat("late-only (24M):", sum(presence_matrix$pattern == "--24"), "\n")
