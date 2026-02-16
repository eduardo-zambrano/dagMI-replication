# helpers.R — Shared utilities for dagMI replication scripts
# Writes LaTeX tabular blocks for \input{} from the paper

#' Format a number for LaTeX table display
#' @param x numeric value
#' @param digits number of decimal places
#' @return formatted string
format_num <- function(x, digits = 2) {
  if (is.na(x) || is.null(x)) return("---")
  formatC(x, format = "f", digits = digits)
}

#' Format a rate (rejection rate, size, power) for display
#' @param x numeric value between 0 and 1
#' @param digits number of decimal places
#' @return formatted string
format_rate <- function(x, digits = 2) {
  if (is.na(x) || is.null(x)) return("---")
  formatC(x, format = "f", digits = digits)
}

#' Format time in seconds for display
#' @param x numeric value (seconds)
#' @param digits number of decimal places
#' @return formatted string
format_time <- function(x, digits = 1) {
  if (is.na(x) || is.null(x)) return("---")
  formatC(x, format = "f", digits = digits)
}

#' Write a LaTeX tabular block to a file
#'
#' Writes only the \\begin{tabular}...\\end{tabular} block.
#' The paper wraps this with \\begin{table}, \\caption, \\label, etc.
#'
#' @param lines character vector of lines inside the tabular environment
#' @param file output file path
#' @param align column alignment string (e.g., "lcccc")
write_tex_table <- function(lines, file, align) {
  out <- c(
    paste0("\\begin{tabular}{", align, "}"),
    "\\toprule",
    lines,
    "\\bottomrule",
    "\\end{tabular}"
  )
  writeLines(out, file)
  cat(sprintf("  -> Wrote %s\n", file))
}

#' Escape special LaTeX characters in a string
#' @param x character string
#' @return escaped string
escape_tex <- function(x) {
  x <- gsub("_", "\\_", x, fixed = TRUE)
  x <- gsub("%", "\\%", x, fixed = TRUE)
  x <- gsub("&", "\\&", x, fixed = TRUE)
  x
}

#' Get the output directory path (relative to repo root)
#' @return path string
output_dir <- function() {
  file.path(find_repo_root(), "output")
}

#' Find the repository root by looking for run_all.R
#' @return path to repo root
find_repo_root <- function() {
  # Try current directory first
  if (file.exists("run_all.R")) return(".")
  # Try parent
  if (file.exists("../run_all.R")) return("..")
  # Try grandparent
  if (file.exists("../../run_all.R")) return("../..")
  # Fallback
  stop("Cannot find repo root. Run scripts from the dagMI-replication directory.")
}

#' Source helpers if not already loaded (for standalone script execution)
#' Call this at the top of each experiment script
ensure_helpers <- function() {
  if (!exists("write_tex_table", mode = "function")) {
    # Find and source helpers.R
    root <- find_repo_root()
    source(file.path(root, "R", "helpers.R"))
  }
}

#' Check if QUICK_MODE is set
#' @return logical
is_quick_mode <- function() {
  exists("QUICK_MODE", envir = .GlobalEnv) && isTRUE(get("QUICK_MODE", envir = .GlobalEnv))
}

#' Get simulation parameters depending on mode
#' @param full_n_sims number of simulations in full mode
#' @param full_B number of bootstrap replicates in full mode
#' @return list with n_sims and B
sim_params <- function(full_n_sims = 200, full_B = 500) {
  if (is_quick_mode()) {
    list(n_sims = 5, B = 30)
  } else {
    list(n_sims = full_n_sims, B = full_B)
  }
}

cat("helpers.R loaded\n")
