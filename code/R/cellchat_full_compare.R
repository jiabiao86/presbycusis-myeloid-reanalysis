#!/usr/bin/env Rscript
# Compare the full CellChatDB.mouse run with the controlled formal run
# (severity-associated M12/limma-positive subset) and report what the full
# database adds.

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
full_dir <- file.path(root, "cellchat_full", "final")
formal_dir <- file.path(root, "cellchat_formal")
ages <- c(3, 12, 24)

read_age <- function(suffix) {
  frames <- lapply(ages, function(age) {
    file <- file.path(full_dir, paste0("CellChat_full_", suffix, "_", age, "M.csv"))
    if (!file.exists(file)) {
      return(NULL)
    }
    read.csv(file, check.names = FALSE)
  })
  frames <- Filter(Negate(is.null), frames)
  if (length(frames) == 0) {
    return(NULL)
  }
  do.call(rbind, frames)
}

network_all <- read_age("network")
network_sig <- read_age("network_significant")
pathway <- read_age("pathway")
summary_frame <- read_age("summary")

if (is.null(network_all) || is.null(summary_frame)) {
  stop("full-run outputs are incomplete; run the merge step first")
}

write.csv(
  network_sig,
  file.path(full_dir, "CellChat_full_significant_all_ages.csv"),
  row.names = FALSE
)
write.csv(
  pathway,
  file.path(full_dir, "CellChat_full_pathway_all_ages.csv"),
  row.names = FALSE
)

controlled <- read.csv(
  file.path(formal_dir, "CellChat_all_communications.csv"),
  check.names = FALSE
)
controlled$run <- "controlled_subset"
if (!is.null(network_sig) && nrow(network_sig) > 0) {
  full_sig <- network_sig
  full_sig$run <- "full_database"
  combined <- rbind(
    controlled[, c("age_months", "source", "target", "interaction_name", "prob", "pval", "run")],
    full_sig[, c("age_months", "source", "target", "interaction_name", "prob", "pval", "run")]
  )
} else {
  combined <- controlled[, c("age_months", "source", "target", "interaction_name", "prob", "pval", "run")]
}
combined <- combined[order(combined$run, combined$age_months, -combined$prob), , drop = FALSE]
write.csv(
  combined,
  file.path(full_dir, "CellChat_full_vs_controlled_significant.csv"),
  row.names = FALSE
)

top_interactions <- do.call(rbind, lapply(ages, function(age) {
  frame <- network_all[network_all$age_months == age, , drop = FALSE]
  frame <- frame[order(-frame$prob), , drop = FALSE]
  head(frame, 15)
}))
write.csv(
  top_interactions,
  file.path(full_dir, "CellChat_full_top_interactions_by_age.csv"),
  row.names = FALSE
)

report_path <- file.path(full_dir, "CellChat_full_report.md")
con <- file(report_path, "w")
writeLines("# Full CellChatDB.mouse run (all 2,019 interactions)", con)
writeLines("", con)
writeLines("## Per-age summary", con)
writeLines("", con)
writeLines(
  paste(
    "| age | interactions evaluated | nonzero edges | significant edges | significant interactions | significant pathways |",
    "| --- | --- | --- | --- | --- | --- |",
    sep = "\n"
  ),
  con
)
for (index in seq_len(nrow(summary_frame))) {
  row <- summary_frame[index, ]
  writeLines(
    paste0(
      "| ", row$age_months, "M | ", row$n_interactions_evaluated,
      " | ", row$n_nonzero_edges,
      " | ", row$n_significant_edges,
      " | ", row$n_significant_interactions,
      " | ", row$n_significant_pathways, " |"
    ),
    con
  )
}
writeLines("", con)

if (!is.null(network_sig) && nrow(network_sig) > 0) {
  writeLines("## Significant interactions (permutation p < 0.05)", con)
  writeLines("", con)
  for (index in seq_len(nrow(network_sig))) {
    row <- network_sig[index, ]
    writeLines(
      paste0(
        "- ", row$age_months, "M: ", row$source, " -> ", row$target, " | ",
        row$interaction_name, " | prob ", signif(row$prob, 4),
        " | p ", signif(row$pval, 4)
      ),
      con
    )
  }
} else {
  writeLines("## Significant interactions", con)
  writeLines("", con)
  writeLines(
    "No interaction reached the permutation threshold in the full database run.",
    con
  )
}
writeLines("", con)
writeLines("## Controlled subset run for comparison", con)
writeLines("", con)
for (index in seq_len(nrow(controlled))) {
  row <- controlled[index, ]
  writeLines(
    paste0(
      "- ", row$age_months, "M: ", row$source, " -> ", row$target, " | ",
      row$interaction_name, " | prob ", signif(row$prob, 4),
      " | p ", signif(row$pval, 4)
    ),
    con
  )
}
close(con)

cat("per-age summary\n")
print(summary_frame)
cat("\nsignificant interactions in the full run:", nrow(network_sig), "\n")
if (nrow(network_sig) > 0) print(network_sig)
cat("\nreport written to", report_path, "\n")
