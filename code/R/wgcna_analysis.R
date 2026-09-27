suppressPackageStartupMessages({
  library(WGCNA)
  library(matrixStats)
})

allowWGCNAThreads(2)

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else "/Users/jijiabiao/Desktop/耳聋课题"
out <- file.path(root, "2026_supplement_analysis")

read_expression <- function(path) {
  frame <- read.csv(path, check.names = FALSE, row.names = 1)
  frame <- as.matrix(frame)
  storage.mode(frame) <- "double"
  frame
}

cochlea <- read_expression(file.path(out, "expression_GSE49543_log2.csv"))
central <- read_expression(file.path(out, "expression_GSE49522_log2.csv"))

common_genes <- intersect(rownames(cochlea), rownames(central))
cochlea_common <- cochlea[common_genes, , drop = FALSE]
central_common <- central[common_genes, , drop = FALSE]

gene_variance <- rowVars(cochlea_common)
variance_rank <- order(gene_variance, decreasing = TRUE)
keep <- variance_rank[seq_len(min(4000, length(variance_rank)))]
cochlea_hvg <- cochlea_common[keep, , drop = FALSE]
central_hvg <- central_common[keep, , drop = FALSE]

expr1 <- t(cochlea_hvg)
expr2 <- t(central_hvg)
gsg <- goodSamplesGenes(expr1, verbose = 0)
if (!gsg$allOK) {
  expr1 <- expr1[gsg$goodSamples, gsg$goodGenes]
  expr2 <- expr2[gsg$goodSamples, gsg$goodGenes]
}

powers <- c(1:10, seq(from = 12, to = 20, by = 2))
sft <- pickSoftThreshold(expr1, powerVector = powers, verbose = 0)
write.csv(
  sft$fitIndices,
  file.path(out, "wgcna_soft_power.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)
power <- sft$powerEstimate
if (!is.finite(power)) {
  power <- 6
}

network <- blockwiseModules(
  expr1,
  power = power,
  maxBlockSize = 5000,
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 30,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 3
)

module_labels <- data.frame(
  gene = colnames(expr1),
  module = paste0("M", network$colors)
)
write.csv(
  module_labels,
  file.path(out, "wgcna_module_assignment.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

module_eigengenes <- orderMEs(network$MEs)
severity <- c(YC = 0, MA = 1, MP = 2, SP = 3)[sub("-.*$", "", rownames(expr1))]
trait_table <- data.frame(severity = severity)
module_trait <- cor(module_eigengenes, trait_table$severity, use = "pairwise.complete.obs")
module_trait <- data.frame(
  module = rownames(module_trait),
  severity_correlation = module_trait[, 1]
)
write.csv(
  module_trait,
  file.path(out, "wgcna_module_trait_correlation.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

preservation_path <- file.path(out, "wgcna_module_preservation.csv")
preservation_result <- tryCatch(
  {
    multi_expr <- list(
      cochlea = list(data = expr1),
      central = list(data = expr2)
    )
    multi_color <- list(
      cochlea = network$colors
    )
    module_preservation <- modulePreservation(
      multi_expr,
      multi_color,
      referenceNetworks = 1,
      testNetworks = 2,
      nPermutations = 50,
      verbose = 2,
      indent = 0,
      parallel = FALSE
    )
    saveRDS(
      module_preservation,
      file.path(out, "wgcna_module_preservation.rds")
    )
    observed <- as.data.frame(module_preservation$preservation$observed[[1]][[2]])
    zsummary <- as.data.frame(module_preservation$preservation$Z[[1]][[2]])
    log_p_adjusted <- as.data.frame(module_preservation$preservation$log.pBonf[[1]][[2]])
    observed$module <- rownames(observed)
    zsummary$module <- rownames(zsummary)
    log_p_adjusted$module <- rownames(log_p_adjusted)
    merged <- merge(
      zsummary[, c("module", "Zsummary.pres")],
      observed[, c("module", "medianRank.pres")],
      by = "module",
      all = TRUE
    )
    merged <- merge(
      merged,
      log_p_adjusted[, c("module", "log.p.Bonfsummary.pres")],
      by = "module",
      all = TRUE
    )
    data.frame(
      module = merged$module,
      observed_rank = merged$medianRank.pres,
      observed_p = 10^merged$log.p.Bonfsummary.pres,
      zsummary = merged$Zsummary.pres,
      zsummary_p = 10^merged$log.p.Bonfsummary.pres
    )
  },
  error = function(error) {
    data.frame(
      module = "ERROR",
      observed_rank = NA_real_,
      observed_p = NA_real_,
      zsummary = NA_real_,
      zsummary_p = NA_real_,
      error = conditionMessage(error)
    )
  }
)
write.csv(
  preservation_result,
  preservation_path,
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("WGCNA complete; selected power =", power, "\n")
