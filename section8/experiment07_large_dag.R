#!/usr/bin/env Rscript
# Experiment 7: Larger DAG (n=8) — size and power
# Output: output/table09_large_dag.tex (Paper Table 9)

library(dagMI)
if (!exists("write_tex_table", mode = "function")) source("R/helpers.R")

cat("=== Experiment 7: Large DAG (n=8) ===\n")

set.seed(2026)
params <- sim_params(full_n_sims = 50, full_B = 100)
n_sims <- params$n_sims
B <- params$B
sample_sizes <- c(500, 1000)
beta <- 0.5

# 8-node chain: X1 -> X2 -> ... -> X8
n_nodes <- 8
A_chain8 <- matrix(0, nrow = n_nodes, ncol = n_nodes)
for (i in 1:(n_nodes - 1)) A_chain8[i, i + 1] <- 1
node_names <- paste0("X", 1:n_nodes)
g_chain8 <- dag(A_chain8, nodes = node_names)

# 8-node fork: X1 -> X2, X1 -> X3, ..., X1 -> X8
A_fork8 <- matrix(0, nrow = n_nodes, ncol = n_nodes)
for (i in 2:n_nodes) A_fork8[1, i] <- 1
g_fork8 <- dag(A_fork8, nodes = node_names)

generate_chain8 <- function(N, beta) {
  X <- matrix(0, nrow = N, ncol = n_nodes)
  X[, 1] <- rnorm(N)
  for (j in 2:n_nodes) {
    X[, j] <- beta * X[, j - 1] + rnorm(N, sd = sqrt(1 - beta^2))
  }
  colnames(X) <- node_names
  X
}

results <- data.frame(
  test = character(0),
  N = integer(0),
  rejection_rate = numeric(0),
  mean_T = numeric(0),
  mean_time_sec = numeric(0)
)

for (N in sample_sizes) {
  cat(sprintf("\n  N = %d\n", N))

  # Size: test correct DAG (chain)
  rejections_size <- 0
  T_size <- numeric(n_sims)
  times_size <- numeric(n_sims)
  for (sim in 1:n_sims) {
    data_mat <- generate_chain8(N, beta)
    t0 <- proc.time()[3]
    res <- mi_test(data_mat, g_chain8, B = B, ordering = "first",
                   verbose = FALSE, seed = sim)
    times_size[sim] <- proc.time()[3] - t0
    T_size[sim] <- res$statistic
    if (res$decision == "reject") rejections_size <- rejections_size + 1
  }
  cat(sprintf("    Size: %.3f, time: %.1fs\n",
              rejections_size / n_sims, mean(times_size)))
  results <- rbind(results, data.frame(
    test = "Size (chain$_8$, correct)",
    N = N,
    rejection_rate = rejections_size / n_sims,
    mean_T = mean(T_size),
    mean_time_sec = mean(times_size)
  ))

  # Power: test wrong DAG (fork) with chain data
  rejections_power <- 0
  T_power <- numeric(n_sims)
  times_power <- numeric(n_sims)
  for (sim in 1:n_sims) {
    data_mat <- generate_chain8(N, beta)
    t0 <- proc.time()[3]
    res <- mi_test(data_mat, g_fork8, B = B, ordering = "first",
                   verbose = FALSE, seed = sim)
    times_power[sim] <- proc.time()[3] - t0
    T_power[sim] <- res$statistic
    if (res$decision == "reject") rejections_power <- rejections_power + 1
  }
  cat(sprintf("    Power: %.3f, time: %.1fs\n",
              rejections_power / n_sims, mean(times_power)))
  results <- rbind(results, data.frame(
    test = "Power (fork$_8$ vs chain$_8$)",
    N = N,
    rejection_rate = rejections_power / n_sims,
    mean_T = mean(T_power),
    mean_time_sec = mean(times_power)
  ))
}

# Write LaTeX table
lines <- c(
  "& & & & Mean time \\\\",
  "Test & $N$ & Rej.\\ rate & Mean $\\hat{T}$ & (sec/sim) \\\\",
  "\\midrule"
)
for (i in 1:nrow(results)) {
  lines <- c(lines, sprintf("%s & %d & %s & %s & %s \\\\",
    results$test[i],
    results$N[i],
    format_rate(results$rejection_rate[i], 2),
    format_num(results$mean_T[i], 2),
    format_time(results$mean_time_sec[i], 1)))
}

write_tex_table(lines, file.path(output_dir(), "table09_large_dag.tex"), "lcccc")
cat("Experiment 7 complete.\n")
