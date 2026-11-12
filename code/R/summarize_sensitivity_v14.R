#!/usr/bin/env Rscript

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(script_arg) == 0) {
  stop("Run with Rscript so repository paths can be resolved.")
}
script_file <- normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
source(file.path(dirname(script_file), "common_v14.R"))
configure_v14_library()

REPO_ROOT <- v14_repo_root()
PRIMARY <- v14_primary_results_dir(REPO_ROOT)
RESULTS <- file.path(REPO_ROOT, "results")
OUT <- file.path(RESULTS, "sensitivity_summary")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

scenarios <- data.frame(
  scenario = c(
    "primary",
    "exclude_GSM1201568",
    "exclude_GSM1201602",
    "exclude_both_unmatched"
  ),
  label = c(
    "Primary (41 cochlea / 39 IC)",
    "Exclude GSM1201568 (40 / 39)",
    "Exclude GSM1201602 (40 / 39)",
    "Exclude both unmatched (39 / 39)"
  ),
  directory = c(
    PRIMARY,
    file.path(RESULTS, "sensitivity_exclude_GSM1201568"),
    file.path(RESULTS, "sensitivity_exclude_GSM1201602"),
    file.path(RESULTS, "sensitivity_exclude_both_unmatched")
  ),
  stringsAsFactors = FALSE
)

read_shared <- function(path) {
  data <- read.csv(
    file.path(path, "shared_genes.csv"),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  data[data$contrast == "SP_vs_YC", , drop = FALSE]
}

shared <- lapply(scenarios$directory, read_shared)
names(shared) <- scenarios$scenario
universe <- sort(unique(unlist(lapply(shared, `[[`, "gene"))))

membership <- data.frame(gene = universe, stringsAsFactors = FALSE)
for (scenario in names(shared)) {
  membership[[scenario]] <- universe %in% shared[[scenario]]$gene
}

loo_frequency <- read.csv(
  file.path(RESULTS, "leave_one_animal_out", "leave_one_animal_out_gene_frequency.csv"),
  check.names = FALSE,
  stringsAsFactors = FALSE
)
membership$leave_one_animal_out_frequency <- loo_frequency$frequency[
  match(membership$gene, loo_frequency$gene)
]

summary <- do.call(rbind, lapply(seq_len(nrow(scenarios)), function(index) {
  scenario <- scenarios$scenario[index]
  data <- shared[[scenario]]
  primary_genes <- shared$primary$gene
  current_genes <- data$gene
  data.frame(
    scenario = scenario,
    label = scenarios$label[index],
    cochlea_samples = if (scenario %in% c("exclude_both_unmatched")) 39 else if (scenario == "primary") 41 else 40,
    inferior_colliculus_samples = 39,
    shared_gene_count = nrow(data),
    shared_genes = paste(current_genes, collapse = ";"),
    lost_vs_primary = paste(setdiff(primary_genes, current_genes), collapse = ";"),
    gained_vs_primary = paste(setdiff(current_genes, primary_genes), collapse = ";"),
    stringsAsFactors = FALSE
  )
}))

loo_summary <- jsonlite::fromJSON(
  file.path(RESULTS, "leave_one_animal_out", "leave_one_animal_out_summary.json")
)
summary <- rbind(
  summary,
  data.frame(
    scenario = "leave_one_animal_out",
    label = "Leave one animal out (41 models)",
    cochlea_samples = NA_integer_,
    inferior_colliculus_samples = NA_integer_,
    shared_gene_count = paste0(
      loo_summary$shared_gene_count$min, "-",
      loo_summary$shared_gene_count$max,
      "; median ", loo_summary$shared_gene_count$median
    ),
    shared_genes = paste(membership$gene, collapse = ";"),
    lost_vs_primary = "",
    gained_vs_primary = "",
    stringsAsFactors = FALSE
  )
)

write.csv(summary, file.path(OUT, "sensitivity_scenarios.csv"), row.names = FALSE)
write.csv(membership, file.path(OUT, "sensitivity_gene_membership.csv"), row.names = FALSE)
write.csv(
  loo_summary$shared_gene_count$distribution,
  file.path(OUT, "leave_one_animal_out_distribution.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()),
  file.path(OUT, "sessionInfo.txt")
)

cat("Sensitivity summary written to", OUT, "\n")
