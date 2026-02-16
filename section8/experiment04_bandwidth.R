#!/usr/bin/env Rscript
# Experiment 4: Bandwidth sensitivity
# Tests Chain n=3, N=1000 with various bandwidth multipliers
#
# Uses dagMI internal functions (.fit_conditional_dists, .bootstrap_sequential)
# to run the constrained bootstrap with explicit bandwidth control.
# The same bandwidth is used for Q_n, T_obs, and each T_boot replicate.
#
# Output: output/table06_bandwidth.tex (Paper Table 6)

library(dagMI)
if (!exists("write_tex_table", mode = "function")) source("R/helpers.R")

cat("=== Experiment 4: Bandwidth Sensitivity ===\n")

set.seed(2026)
params <- sim_params(full_n_sims = 50, full_B = 100)
n_sims <- params$n_sims
B <- params$B
N <- 1000
beta <- 0.5

# Chain n=3: X1 -> X2 -> X3
A_chain3 <- matrix(c(0,1,0, 0,0,1, 0,0,0), 3, 3, byrow = TRUE)
g_chain3 <- dag(A_chain3, nodes = c("X1","X2","X3"))

# Fork n=3: X1 <- X2 -> X3 (wrong DAG for power)
A_fork3 <- matrix(c(0,0,0, 1,0,1, 0,0,0), 3, 3, byrow = TRUE)
g_fork3 <- dag(A_fork3, nodes = c("X1","X2","X3"))

generate_chain3 <- function(N, beta) {
  X1 <- rnorm(N)
  X2 <- beta * X1 + rnorm(N, sd = sqrt(1 - beta^2))
  X3 <- beta * X2 + rnorm(N, sd = sqrt(1 - beta^2))
  m <- cbind(X1, X2, X3); colnames(m) <- c("X1","X2","X3"); m
}

# Silverman bandwidth: 1.06 * sd(x) * n^{-1/5}
silverman_bw <- function(x) 1.06 * sd(x) * length(x)^(-1/5)

# --- Constrained bootstrap with explicit bandwidth ---
# Uses dagMI internals to run the same bootstrap algorithm as constrained_bootstrap()
# but with a user-specified bandwidth vector for Q_n and test statistic computation.
bootstrap_with_bandwidth <- function(data_mat, g, B, bw_vec, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n_obs <- nrow(data_mat)
  n_vars <- ncol(data_mat)

  # Compute Q_n and ordering with custom bandwidth
  qn_result <- compute_qn(data_mat, g, bandwidth = bw_vec)
  ordering <- qn_result$ordering

  # Compute h_params$a from data (mirrors constrained_bootstrap)
  data_ordered <- data_mat[, ordering, drop = FALSE]
  avg_var <- mean(apply(data_ordered, 2, var))
  h_params <- list(a = 1 / (4 * max(avg_var, .Machine$double.eps)))

  # Compute T_obs with custom bandwidth
  t_obs <- compute_test_stat(data_mat, g, h_params = h_params,
                              qn_result = qn_result,
                              bandwidth = bw_vec)$statistic

  # Fit conditional distributions (uses dagMI internal)
  cond_kdes <- dagMI:::.fit_conditional_dists(data_ordered, g, ordering)

  # Run bootstrap with custom bandwidth passed as orig_bandwidth
  t_boot <- dagMI:::.bootstrap_sequential(
    data_ordered, g, ordering, cond_kdes, qn_result,
    "gaussian", h_params, B, n_obs, n_vars, FALSE,
    orig_bandwidth = bw_vec)

  p_value <- mean(t_boot >= t_obs)
  list(t_obs = t_obs, t_boot = t_boot, p_value = p_value)
}

# --- Run bandwidth sensitivity experiment ---
bw_mults <- c(0.5, 1.0, 1.5, 2.0)
bw_labels <- c("$0.5 \\times$ Silverman", "$1.0 \\times$ Silverman",
               "$1.5 \\times$ Silverman", "$2.0 \\times$ Silverman",
               "Cross-validation")

results <- data.frame(
  label = character(0),
  size = numeric(0),
  power = numeric(0)
)

for (k in seq_along(bw_mults)) {
  mult <- bw_mults[k]
  cat(sprintf("  %s: ", bw_labels[k]))

  size_rej <- 0
  power_rej <- 0
  for (sim in 1:n_sims) {
    data_mat <- generate_chain3(N, beta)
    bw_vec <- sapply(1:3, function(j) silverman_bw(data_mat[, j]) * mult)

    # Size test (correct DAG)
    res_size <- bootstrap_with_bandwidth(data_mat, g_chain3, B, bw_vec, seed = sim)
    if (res_size$p_value < 0.05) size_rej <- size_rej + 1

    # Power test (wrong DAG)
    res_power <- bootstrap_with_bandwidth(data_mat, g_fork3, B, bw_vec, seed = sim + 10000)
    if (res_power$p_value < 0.05) power_rej <- power_rej + 1
  }

  size_rate <- size_rej / n_sims
  power_rate <- power_rej / n_sims
  cat(sprintf("size=%.3f, power=%.2f\n", size_rate, power_rate))
  results <- rbind(results, data.frame(
    label = bw_labels[k], size = size_rate, power = power_rate))
}

# Cross-validation: use default Silverman bandwidth (via mi_test)
cat(sprintf("  %s: ", bw_labels[5]))
cv_size_rej <- 0
cv_power_rej <- 0
for (sim in 1:n_sims) {
  data_mat <- generate_chain3(N, beta)
  res_size <- mi_test(data_mat, g_chain3, B = B, ordering = "first",
                      verbose = FALSE, seed = sim)
  if (res_size$decision == "reject") cv_size_rej <- cv_size_rej + 1

  res_power <- mi_test(data_mat, g_fork3, B = B, ordering = "first",
                       verbose = FALSE, seed = sim)
  if (res_power$decision == "reject") cv_power_rej <- cv_power_rej + 1
}
cv_size_rate <- cv_size_rej / n_sims
cv_power_rate <- cv_power_rej / n_sims
cat(sprintf("size=%.3f, power=%.2f\n", cv_size_rate, cv_power_rate))
results <- rbind(results, data.frame(
  label = bw_labels[5], size = cv_size_rate, power = cv_power_rate))

# --- Write LaTeX table ---
lines <- c(
  "Bandwidth & Size (true DAG) & Power (wrong DAG) \\\\",
  "\\midrule"
)
for (i in 1:nrow(results)) {
  lines <- c(lines, sprintf("%s & %s & %s \\\\",
    results$label[i],
    format_rate(results$size[i], 3),
    format_rate(results$power[i], 2)))
}

write_tex_table(lines, file.path(output_dir(), "table06_bandwidth.tex"), "lccc")
cat("Experiment 4 complete.\n")
