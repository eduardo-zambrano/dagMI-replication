#!/usr/bin/env Rscript
# Section 9.1: Returns to Schooling with IV (Card 1995)
# Output: output/section9_iv_results.tex (inline values for paper narrative)

library(dagMI)
if (!exists("write_tex_table", mode = "function")) source("R/helpers.R")

cat("=== Section 9.1: IV Returns to Schooling ===\n")

set.seed(2026)
params <- sim_params(full_n_sims = 1, full_B = 500)
B <- params$B

# --- Load or simulate data ---
# Try wooldridge::card (Card 1995) — classic IV: nearc4 -> educ -> lwage
# Fall back to simulated data if unavailable
data_loaded <- FALSE

if (requireNamespace("wooldridge", quietly = TRUE)) {
  tryCatch({
    data("card", package = "wooldridge", envir = environment())
    d <- get("card", envir = environment())
    # IV structure: nearc4 (proximity to 4-year college) -> educ -> lwage
    Z <- d$nearc4           # instrument: grew up near 4-year college
    D <- d$educ             # years of education
    Y <- d$lwage            # log hourly wage
    # Remove rows with NA
    complete <- complete.cases(Z, D, Y)
    Z <- scale(Z[complete])[, 1]
    D <- scale(D[complete])[, 1]
    Y <- scale(Y[complete])[, 1]
    data_mat <- cbind(Z = Z, D = D, Y = Y)
    data_loaded <- TRUE
    cat(sprintf("  Loaded Card (1995) data from wooldridge package (N=%d)\n", sum(complete)))
  }, error = function(e) {
    cat(sprintf("  wooldridge data load failed: %s\n", e$message))
  })
}

if (!data_loaded) {
  # Simulate IV-valid data: Z -> D -> Y (Z _||_ Y | D)
  cat("  Using simulated IV data\n")
  N <- 1000
  Z <- rnorm(N)                              # instrument
  D <- 0.3 * Z + rnorm(N)                    # schooling
  Y <- 0.5 * D + rnorm(N)                    # log wage
  data_mat <- cbind(Z = Z, D = D, Y = Y)
}

# --- Define IV DAG: Z -> D -> Y ---
A_iv <- matrix(c(0,1,0, 0,0,1, 0,0,0), 3, 3, byrow = TRUE)
g_iv <- dag(A_iv, nodes = c("Z", "D", "Y"))

# --- Run MI test ---
cat("  Running MI test...\n")
res <- mi_test(data_mat, g_iv, B = B, ordering = "first", verbose = FALSE, seed = 1)

cat(sprintf("  T_hat = %.2f, p-value = %.2f, decision = %s\n",
            res$statistic, res$p_value, res$decision))

# --- Write results as LaTeX commands ---
out_lines <- c(
  sprintf("\\newcommand{\\ivThat}{%.2f}", res$statistic),
  sprintf("\\newcommand{\\ivPvalue}{%.2f}", res$p_value),
  sprintf("\\newcommand{\\ivDecision}{%s}",
          ifelse(res$decision == "reject", "reject", "does not reject"))
)
writeLines(out_lines, file.path(output_dir(), "section9_iv_results.tex"))
cat(sprintf("  -> Wrote %s\n", file.path(output_dir(), "section9_iv_results.tex")))
cat("Section 9.1 complete.\n")
