#!/usr/bin/env Rscript

root <- "/Users/jijiabiao/Desktop/耳聋课题/2026_supplement_analysis"
database_path <- file.path(root, "data_CellChatDB_mouse.rda")
output_dir <- file.path(root, "cellchat_enhancement")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

environment <- new.env()
load(database_path, envir = environment)
database <- get(ls(environment)[1], envir = environment)

write.csv(
  database$interaction,
  file.path(output_dir, "cellchatdb_mouse_interactions.csv"),
  row.names = FALSE
)

complex <- database$complex
complex$complex_name <- rownames(complex)
write.csv(
  complex,
  file.path(output_dir, "cellchatdb_mouse_complexes.csv"),
  row.names = FALSE
)

cat("Exported", nrow(database$interaction), "interactions and", nrow(complex), "complexes\n")
