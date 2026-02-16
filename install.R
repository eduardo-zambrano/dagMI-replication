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
remotes::install_github("eduardo-zambrano/dagMI", upgrade = "never")

# Check dagMI loads
if (requireNamespace("dagMI", quietly = TRUE)) {
  cat(sprintf("dagMI %s installed successfully.\n\n",
              as.character(packageVersion("dagMI"))))
} else {
  stop("Failed to install dagMI. See errors above.")
}

# Optional dependencies for Experiment 3 (method comparison)
optional_pkgs <- c("pcalg", "bnlearn")
cat("Checking optional packages for Experiment 3 (method comparison)...\n")
for (pkg in optional_pkgs) {
  if (requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("  %s: installed\n", pkg))
  } else {
    cat(sprintf("  %s: NOT installed (Experiment 3 will skip %s columns)\n", pkg, pkg))
    cat(sprintf("    To install: install.packages('%s')\n", pkg))
  }
}

# Optional: AER for IV data
if (requireNamespace("AER", quietly = TRUE)) {
  cat("  AER: installed (for Section 9.1 IV data)\n")
} else {
  cat("  AER: NOT installed (Section 9.1 will use bundled data)\n")
}

cat("\n=== Installation complete ===\n")
