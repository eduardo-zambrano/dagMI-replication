#!/usr/bin/env Rscript
# Experiment 8: Robustness — misspecification and moment boundary cases
# Output: output/table10_robustness.tex (Paper Table 10)

library(dagMI)
if (!exists("write_tex_table", mode = "function")) source("R/helpers.R")

cat("=== Experiment 8: Robustness checks ===\n")

set.seed(2026)
params <- sim_params(full_n_sims = 50, full_B = 100)
n_sims <- params$n_sims
B <- params$B
N <- 1000
beta <- 0.5

# DAG: Chain X1 -> X2 -> X3
A_chain <- matrix(c(0, 1, 0,
                     0, 0, 1,
                     0, 0, 0), nrow = 3, byrow = TRUE)
g_chain <- dag(A_chain, nodes = c("X1", "X2", "X3"))

run_dgp <- function(dgp_name, generate_fn) {
  cat(sprintf("  %s: ", dgp_name))
  rejections <- 0
  T_vals <- numeric(n_sims)
  for (sim in 1:n_sims) {
    data_mat <- generate_fn(N)
    res <- tryCatch(
      mi_test(data_mat, g_chain, B = B, ordering = "first",
              verbose = FALSE, seed = sim),
      error = function(e) list(statistic = NA, decision = "error")
    )
    T_vals[sim] <- res$statistic
    if (!is.na(res$decision) && res$decision == "reject")
      rejections <- rejections + 1
  }
  n_valid <- sum(!is.na(T_vals))
  cat(sprintf("rej = %.3f, mean T = %.3f\n",
              rejections / n_valid, mean(T_vals, na.rm = TRUE)))
  data.frame(
    dgp = dgp_name,
    rejection_rate = rejections / n_valid,
    mean_T = mean(T_vals, na.rm = TRUE),
    sd_T = sd(T_vals, na.rm = TRUE)
  )
}

# Gaussian baseline
res_gauss <- run_dgp("Gaussian (baseline)", function(N) {
  X1 <- rnorm(N)
  X2 <- beta * X1 + rnorm(N, sd = sqrt(1 - beta^2))
  X3 <- beta * X2 + rnorm(N, sd = sqrt(1 - beta^2))
  m <- cbind(X1, X2, X3); colnames(m) <- c("X1","X2","X3"); m
})

# t_10 errors
res_t10 <- run_dgp("$t_{10}$ errors", function(N) {
  X1 <- rt(N, df = 10) / sqrt(10 / 8)
  X2 <- beta * X1 + rt(N, df = 10) / sqrt(10 / 8) * sqrt(1 - beta^2)
  X3 <- beta * X2 + rt(N, df = 10) / sqrt(10 / 8) * sqrt(1 - beta^2)
  m <- cbind(X1, X2, X3); colnames(m) <- c("X1","X2","X3"); m
})

# t_3 errors (moment violation)
res_t3 <- run_dgp("$t_3$ errors", function(N) {
  X1 <- rt(N, df = 3)
  X2 <- beta * X1 + rt(N, df = 3) * sqrt(1 - beta^2) / sqrt(3)
  X3 <- beta * X2 + rt(N, df = 3) * sqrt(1 - beta^2) / sqrt(3)
  m <- cbind(X1, X2, X3); colnames(m) <- c("X1","X2","X3"); m
})

# Nonlinear tanh DGP
res_tanh <- run_dgp("Nonlinear $\\tanh$ DGP", function(N) {
  X1 <- rnorm(N)
  X2 <- tanh(2 * beta * X1) + rnorm(N, sd = 0.3)
  X3 <- tanh(2 * beta * X2) + rnorm(N, sd = 0.3)
  m <- cbind(X1, X2, X3); colnames(m) <- c("X1","X2","X3"); m
})

results <- rbind(res_gauss, res_t10, res_t3, res_tanh)

# Tail classification
tails <- c("Normal", "Moderate", "Heavy", "Normal")

# Write LaTeX table
lines <- c(
  "DGP & Rej.\\ rate & Mean $\\hat{T}$ & SD $\\hat{T}$ & Tails \\\\",
  "\\midrule"
)
for (i in 1:nrow(results)) {
  lines <- c(lines, sprintf("%s & %s & %s & %s & %s \\\\",
    results$dgp[i],
    format_rate(results$rejection_rate[i], 2),
    format_num(results$mean_T[i], 2),
    format_num(results$sd_T[i], 2),
    tails[i]))
}

write_tex_table(lines, file.path(output_dir(), "table10_robustness.tex"), "lcccc")
cat("Experiment 8 complete.\n")
