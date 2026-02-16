#!/usr/bin/env Rscript
# Experiment 6: Power as a function of effect size beta
# Tests Chain -> Fork (Markov equivalent) at N=1000 for varying beta
# Output: output/table08_power_beta.tex (Paper Table 8)

library(dagMI)
if (!exists("write_tex_table", mode = "function")) source("R/helpers.R")

cat("=== Experiment 6: Power vs. effect size ===\n")

set.seed(2026)
params <- sim_params(full_n_sims = 50, full_B = 100)
n_sims <- params$n_sims
B <- params$B
N <- 1000
betas <- c(0.05, 0.1, 0.2, 0.3, 0.5, 0.7)

# True DAG: Chain X1 -> X2 -> X3
# Tested DAG: Fork X1 <- X2 -> X3 (Markov equivalent)
A_fork <- matrix(c(0, 0, 0,
                    1, 0, 1,
                    0, 0, 0), nrow = 3, byrow = TRUE)
g_fork <- dag(A_fork, nodes = c("X1", "X2", "X3"))

results <- data.frame(
  beta = numeric(0),
  rejection_rate = numeric(0),
  mean_T = numeric(0),
  sd_T = numeric(0)
)

for (beta in betas) {
  cat(sprintf("  beta = %.2f: ", beta))
  rejections <- 0
  T_vals <- numeric(n_sims)

  for (sim in 1:n_sims) {
    X1 <- rnorm(N)
    X2 <- beta * X1 + rnorm(N, sd = sqrt(1 - beta^2))
    X3 <- beta * X2 + rnorm(N, sd = sqrt(1 - beta^2))
    data_mat <- cbind(X1, X2, X3)
    colnames(data_mat) <- c("X1", "X2", "X3")

    res <- mi_test(data_mat, g_fork, B = B, ordering = "first",
                   verbose = FALSE, seed = sim)
    T_vals[sim] <- res$statistic
    if (res$decision == "reject") rejections <- rejections + 1
  }

  results <- rbind(results, data.frame(
    beta = beta,
    rejection_rate = rejections / n_sims,
    mean_T = mean(T_vals),
    sd_T = sd(T_vals)
  ))
  cat(sprintf("power = %.3f, mean T = %.3f\n", rejections / n_sims, mean(T_vals)))
}

# Write LaTeX table
lines <- c(
  "$\\beta$ & Rej.\\ rate & Mean $\\hat{T}$ & SD $\\hat{T}$ \\\\",
  "\\midrule"
)
for (i in 1:nrow(results)) {
  lines <- c(lines, sprintf("%s & %s & %s & %s \\\\",
    format_num(results$beta[i], 2),
    format_rate(results$rejection_rate[i], 2),
    format_num(results$mean_T[i], 2),
    format_num(results$sd_T[i], 2)))
}

write_tex_table(lines, file.path(output_dir(), "table08_power_beta.tex"), "lccc")
cat("Experiment 6 complete.\n")
