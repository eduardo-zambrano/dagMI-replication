#!/usr/bin/env Rscript
# install.R — Install dagMI and dependencies for replication
#
# Usage: Rscript install.R

cat("=== Installing dagMI replication dependencies ===\n\n")

# Core dependency: dagMI package
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes", repos = "https://cloud.r-project.org")
}

cat("Installing dagMI from GitHub...\n")
remotes::install_github("eduardo-zambrano/dagMI@v0.2.0", upgrade = "never")

# Check dagMI loads
if (requireNamespace("dagMI", quietly = TRUE)) {
  cat(sprintf("  dagMI %s installed successfully.\n\n",
              as.character(packageVersion("dagMI"))))
} else {
  stop("Failed to install dagMI. See errors above.")
}

# Required: wooldridge (Card 1995 data for Section 9.1)
cat("Installing required packages...\n")
for (pkg in c("wooldridge")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %s %s: installed\n", pkg, as.character(packageVersion(pkg))))
  } else {
    cat(sprintf("  %s: FAILED to install\n", pkg))
  }
}

# Optional dependencies for Experiment 3 (method comparison)
cat("\nChecking optional packages for Experiment 3 (method comparison)...\n")
for (pkg in c("pcalg", "CAM")) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %s: installed\n", pkg))
  } else {
    cat(sprintf("  %s: NOT installed (Experiment 3 will skip this method)\n", pkg))
    cat(sprintf("    To install: install.packages('%s')\n", pkg))
  }
}

cat("\n=== Installation complete ===\n")
