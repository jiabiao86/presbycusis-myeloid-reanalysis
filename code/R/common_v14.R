v14_script_file <- function() {
  script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
  if (length(script_arg) == 0) {
    stop("Run this script with Rscript so the repository root can be resolved.")
  }
  normalizePath(sub("^--file=", "", script_arg[[1]]), mustWork = TRUE)
}

v14_repo_root <- function() {
  configured <- trimws(Sys.getenv("V14_REPO_ROOT", unset = ""))
  if (nzchar(configured)) {
    return(normalizePath(configured, mustWork = TRUE))
  }
  normalizePath(
    dirname(dirname(dirname(v14_script_file()))),
    mustWork = TRUE
  )
}

v14_path <- function(env_var, default) {
  configured <- trimws(Sys.getenv(env_var, unset = ""))
  if (nzchar(configured)) {
    return(normalizePath(configured, mustWork = FALSE))
  }
  normalizePath(default, mustWork = FALSE)
}

v14_relative_path <- function(path, root) {
  normalized_path <- normalizePath(path, mustWork = FALSE)
  normalized_root <- normalizePath(root, mustWork = FALSE)
  prefix <- paste0(normalized_root, .Platform$file.sep)
  if (startsWith(normalized_path, prefix)) {
    return(substring(normalized_path, nchar(prefix) + 1L))
  }
  normalized_path
}

v14_prepared_dir <- function(repo_root) {
  v14_path(
    "V14_PREPARED_DIR",
    file.path(repo_root, "inputs", "prepared")
  )
}

v14_primary_results_dir <- function(repo_root) {
  v14_path(
    "V14_PRIMARY_RESULTS_DIR",
    file.path(repo_root, "results", "primary")
  )
}

v14_sensitivity_results_dir <- function(repo_root) {
  v14_path(
    "V14_SENSITIVITY_RESULTS_DIR",
    file.path(repo_root, "results", "sensitivity_exclude_GSM1201568")
  )
}

v14_single_cell_dir <- function(repo_root) {
  v14_path(
    "V14_SINGLE_CELL_DIR",
    file.path(repo_root, "inputs", "single_cell", "GSE274279_10x")
  )
}

v14_derived_dir <- function(repo_root) {
  v14_path(
    "V14_DERIVED_DIR",
    file.path(repo_root, "results", "derived")
  )
}

v14_figure_dir <- function(repo_root) {
  v14_path(
    "V14_FIGURE_DIR",
    file.path(repo_root, "figures", "source")
  )
}

configure_v14_library <- function() {
  configured <- trimws(Sys.getenv("V14_R_LIB", unset = ""))
  if (nzchar(configured)) {
    .libPaths(c(normalizePath(configured, mustWork = TRUE), .libPaths()))
  }
}
