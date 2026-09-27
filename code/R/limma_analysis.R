suppressPackageStartupMessages({
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args) >= 1) args[[1]] else "/Users/jijiabiao/Desktop/耳聋课题"
out <- file.path(root, "2026_supplement_analysis")
expression_path <- file.path(root, "2022耳聋耳鸣", "public", "exp.txt")

expression <- read.delim(
  expression_path,
  check.names = FALSE,
  row.names = 1,
  stringsAsFactors = FALSE
)
expression <- as.matrix(expression)
storage.mode(expression) <- "double"
expression <- log2(expression + 1)

sample_names <- colnames(expression)
groups <- factor(sub("-.*$", "", sample_names), levels = c("YC", "MA", "MP", "SP"))
design <- model.matrix(~ 0 + groups)
colnames(design) <- levels(groups)
fit <- lmFit(expression, design)

contrast_matrix <- makeContrasts(
  MA_vs_YC = MA - YC,
  MP_vs_MA = MP - MA,
  MP_vs_YC = MP - YC,
  SP_vs_MA = SP - MA,
  SP_vs_MP = SP - MP,
  SP_vs_YC = SP - YC,
  levels = design
)
fit_contrasts <- eBayes(contrasts.fit(fit, contrast_matrix))

for (comparison in colnames(contrast_matrix)) {
  table <- topTable(
    fit_contrasts,
    coef = comparison,
    number = Inf,
    sort.by = "P"
  )
  table$gene <- rownames(table)
  table <- table[, c("gene", setdiff(colnames(table), "gene"))]
  write.csv(
    table,
    file.path(out, paste0("limma_", comparison, ".csv")),
    row.names = FALSE,
    fileEncoding = "UTF-8"
  )
}

severity <- as.numeric(groups) - 1
severity_design <- model.matrix(~ severity)
severity_fit <- eBayes(lmFit(expression, severity_design))
severity_table <- topTable(
  severity_fit,
  coef = "severity",
  number = Inf,
  sort.by = "P"
)
severity_table$gene <- rownames(severity_table)
severity_table <- severity_table[, c("gene", setdiff(colnames(severity_table), "gene"))]
write.csv(
  severity_table,
  file.path(out, "limma_severity_trend.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

summary <- lapply(colnames(contrast_matrix), function(comparison) {
  table <- topTable(fit_contrasts, coef = comparison, number = Inf, sort.by = "P")
  data.frame(
    comparison = comparison,
    FDR_0.05_log2FC_1_up = sum(table$adj.P.Val <= 0.05 & table$logFC >= 1),
    FDR_0.05_log2FC_1_down = sum(table$adj.P.Val <= 0.05 & table$logFC <= -1),
    FDR_0.05_log2FC_0.585_up = sum(table$adj.P.Val <= 0.05 & table$logFC >= 0.585),
    FDR_0.05_log2FC_0.585_down = sum(table$adj.P.Val <= 0.05 & table$logFC <= -0.585)
  )
})
summary <- do.call(rbind, summary)
write.csv(
  summary,
  file.path(out, "limma_summary.csv"),
  row.names = FALSE,
  fileEncoding = "UTF-8"
)

cat("limma analysis complete\n")
