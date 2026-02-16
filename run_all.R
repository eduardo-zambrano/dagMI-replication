#!/usr/bin/env Rscript
# run_all.R — Master script for dagMI replication
#
# Usage:
#   Rscript run_all.R              # Full run (publication quality, ~15-20 hours)
#   Rscript run_all.R --quick      # Quick validation (~25 min)
#   Rscript run_all.R --section 7  # Only Section 7
#   Rscript run_all.R --section 8  # Only Section 8
#   Rscript run_all.R --section 9  # Only Section 9
#   Rscript run_all.R --experiment 3  # Only Experiment 3

args <- commandArgs(trailingOnly = TRUE)

# Parse arguments
QUICK_MODE <- "--quick" %in% args
section_filter <- NULL
experiment_filter <- NULL

if ("--section" %in% args) {
  idx <- which(args == "--section")
  if (idx < length(args)) section_filter <- as.integer(args[idx + 1])
}
if ("--experiment" %in% args) {
  idx <- which(args == "--experiment")
  if (idx < length(args)) experiment_filter <- as.integer(args[idx + 1])
}

# Set QUICK_MODE as global variable for scripts to read
assign("QUICK_MODE", QUICK_MODE, envir = .GlobalEnv)

cat("========================================\n")
cat("dagMI Replication Package\n")
cat("========================================\n")
cat(sprintf("Mode: %s\n", ifelse(QUICK_MODE, "QUICK (validation)", "FULL (publication)")))
if (!is.null(section_filter)) cat(sprintf("Section filter: %d\n", section_filter))
if (!is.null(experiment_filter)) cat(sprintf("Experiment filter: %d\n", experiment_filter))
cat(sprintf("Started: %s\n", Sys.time()))
cat("========================================\n\n")

# Ensure output directory exists
if (!dir.exists("output")) dir.create("output")

# Source helpers into global environment so scripts can use them
source("R/helpers.R")

# Define all scripts
scripts <- list(
  list(section = 7, experiment = NA,
       file = "section7/benchmark_computation.R",
       name = "Table 1: Computational benchmarks"),
  list(section = 8, experiment = 1,
       file = "section8/experiment01_size.R",
       name = "Table 2: Empirical size"),
  list(section = 8, experiment = 2,
       file = "section8/experiment02_power.R",
       name = "Table 3: Power MI vs CI"),
  list(section = 8, experiment = 3,
       file = "section8/experiment03_comparison.R",
       name = "Table 4: Method comparison"),
  list(section = 8, experiment = 4,
       file = "section8/experiment04_bandwidth.R",
       name = "Table 6: Bandwidth sensitivity"),
  list(section = 8, experiment = 5,
       file = "section8/experiment05_ordering.R",
       name = "Table 7: Ordering strategies"),
  list(section = 8, experiment = 6,
       file = "section8/experiment06_power_beta.R",
       name = "Table 8: Power vs beta"),
  list(section = 8, experiment = 7,
       file = "section8/experiment07_large_dag.R",
       name = "Table 9: Large DAG n=8"),
  list(section = 8, experiment = 8,
       file = "section8/experiment08_robustness.R",
       name = "Table 10: Robustness"),
  list(section = 9, experiment = NA,
       file = "section9/iv_returns.R",
       name = "Section 9.1: IV returns"),
  list(section = 9, experiment = NA,
       file = "section9/gene_network.R",
       name = "Table 11: Gene network")
)

# Apply filters
if (!is.null(section_filter)) {
  scripts <- scripts[sapply(scripts, function(s) s$section == section_filter)]
}
if (!is.null(experiment_filter)) {
  scripts <- scripts[sapply(scripts, function(s) !is.na(s$experiment) && s$experiment == experiment_filter)]
}

if (length(scripts) == 0) {
  cat("No scripts match the specified filters.\n")
  quit("no", status = 1)
}

cat(sprintf("Running %d script(s):\n", length(scripts)))
for (s in scripts) cat(sprintf("  - %s\n", s$name))
cat("\n")

# Run each script
timings <- data.frame(
  script = character(0),
  time_sec = numeric(0),
  status = character(0)
)

for (s in scripts) {
  cat(sprintf("\n--- %s ---\n", s$name))
  t0 <- proc.time()[3]

  status <- tryCatch({
    source(s$file, local = FALSE)
    "OK"
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    sprintf("FAILED: %s", e$message)
  })

  elapsed <- proc.time()[3] - t0
  timings <- rbind(timings, data.frame(
    script = s$name,
    time_sec = elapsed,
    status = status
  ))
  cat(sprintf("  Time: %.1f seconds [%s]\n", elapsed, status))
}

# Summary
cat("\n========================================\n")
cat("SUMMARY\n")
cat("========================================\n")
cat(sprintf("Finished: %s\n", Sys.time()))
cat(sprintf("Total time: %.1f seconds (%.1f minutes)\n",
    sum(timings$time_sec), sum(timings$time_sec) / 60))
cat("\n")
print(timings[, c("script", "time_sec", "status")])

# Check output files
cat("\n--- Output files ---\n")
output_files <- list.files("output", pattern = "\\.tex$", full.names = TRUE)
for (f in output_files) {
  cat(sprintf("  %s (%d bytes)\n", f, file.info(f)$size))
}

n_failed <- sum(timings$status != "OK")
if (n_failed > 0) {
  cat(sprintf("\nWARNING: %d script(s) failed.\n", n_failed))
  quit("no", status = 1)
} else {
  cat("\nAll scripts completed successfully.\n")
  cat("To update paper tables: cp output/*.tex ../DAG_Testing/tables/\n")
}
