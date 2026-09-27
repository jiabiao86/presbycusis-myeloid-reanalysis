#!/usr/bin/env Rscript

root <- "/Users/jijiabiao/Desktop/耳聋课题"
lib <- file.path(root, "Rlibs")
source_dir <- "/private/tmp/cellchat_build/sqjin-CellChat-e4f6862"

dir.create(lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(lib, .libPaths()))
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", lib = lib, type = "source")
}

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", lib = lib, type = "source")
}

remotes::install_local(
  source_dir,
  lib = lib,
  dependencies = TRUE,
  upgrade = "never",
  force = TRUE,
  type = "source",
  INSTALL_opts = "--no-multiarch"
)

cat("CellChat installed at", find.package("CellChat", lib.loc = lib), "\n")
