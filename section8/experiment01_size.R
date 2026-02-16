#!/usr/bin/env Rscript
# Experiment 1: Empirical size under correct DAG specification
# Tests 4 DAG types x 4 sample sizes, constrained vs naive bootstrap
# Output: output/table02_size.tex (Paper Table 2)

library(dagMI)
if (!exists("write_tex_table", mode = "function")) source("R/helpers.R")

cat("=== Experiment 1: Empirical Size ===\n")

set.seed(2026)
params <- sim_params(full_n_sims = 50, full_B = 100)
n_sims <- params$n_sims
B <- params$B
sample_sizes <- c(100, 300, 500, 1000)
beta <- 0.5

# --- Define DAGs ---

# Chain n=3: X1 -> X2 -> X3
A_chain3 <- matrix(c(0,1,0, 0,0,1, 0,0,0), 3, 3, byrow = TRUE)
g_chain3 <- dag(A_chain3, nodes = c("X1","X2","X3"))

# Fork n=3: X2 -> X1, X2 -> X3
A_fork3 <- matrix(c(0,0,0, 1,0,1, 0,0,0), 3, 3, byrow = TRUE)
g_fork3 <- dag(A_fork3, nodes = c("X1","X2","X3"))

# Diamond n=4: X1 -> X2, X1 -> X3, X2 -> X4, X3 -> X4
A_diamond <- matrix(c(0,1,1,0, 0,0,0,1, 0,0,0,1, 0,0,0,0), 4, 4, byrow = TRUE)
g_diamond <- dag(A_diamond, nodes = c("X1","X2","X3","X4"))

# Random n=5: X1->X2, X1->X3, X2->X4, X3->X5, X4->X5
A_rand5 <- matrix(0, 5, 5)
A_rand5[1,2] <- 1; A_rand5[1,3] <- 1; A_rand5[2,4] <- 1
A_rand5[3,5] <- 1; A_rand5[4,5] <- 1
g_rand5 <- dag(A_rand5, nodes = paste0("X", 1:5))

# --- Data generators ---
generate_from_dag <- function(N, adj_mat, beta) {
  n_nodes <- nrow(adj_mat)
  X <- matrix(0, N, n_nodes)
  X[, 1] <- rnorm(N)
  for (j in 2:n_nodes) {
    parents <- which(adj_mat[, j] == 1)
    if (length(parents) > 0) {
      signal <- rowSums(X[, parents, drop = FALSE] * beta / length(parents))
      noise_sd <- sqrt(max(1 - beta^2, 0.1))
      X[, j] <- signal + rnorm(N, sd = noise_sd)
    } else {
      X[, j] <- rnorm(N)
    }
  }
  colnames(X) <- paste0("X", 1:n_nodes)
  X
}

# --- Run experiments ---
dag_configs <- list(
  list(name = "Chain ($n=3$)", dag = g_chain3, adj = A_chain3),
  list(name = "Fork ($n=3$)", dag = g_fork3, adj = A_fork3),
  list(name = "Diamond ($n=4$)", dag = g_diamond, adj = A_diamond),
  list(name = "Random ($n=5$)", dag = g_rand5, adj = A_rand5)
)

results <- list()

for (cfg in dag_configs) {
  cat(sprintf("  %s:\n", cfg$name))
  for (N in sample_sizes) {
    rejections <- 0
    for (sim in 1:n_sims) {
      data_mat <- generate_from_dag(N, cfg$adj, beta)
      res <- mi_test(data_mat, cfg$dag, B = B, ordering = "first",
                     verbose = FALSE, seed = sim)
      if (res$decision == "reject") rejections <- rejections + 1
    }
    rate <- rejections / n_sims
    cat(sprintf("    N=%d: %.3f\n", N, rate))
    results <- c(results, list(data.frame(
      dag = cfg$name, N = N, rate = rate, bootstrap = "Constrained"
    )))
  }
}

# Naive bootstrap comparison (chain only)
# Naive = resample rows iid, ignoring DAG structure
cat("  Chain (naive bootstrap):\n")
for (N in sample_sizes) {
  rejections <- 0
  for (sim in 1:n_sims) {
    data_mat <- generate_from_dag(N, A_chain3, beta)
    # Compute observed test statistic
    t_obs <- compute_test_stat(data_mat, g_chain3)$statistic
    # Naive bootstrap: resample rows iid
    t_boot <- numeric(B)
    for (b in 1:B) {
      idx <- sample(N, N, replace = TRUE)
      boot_data <- data_mat[idx, , drop = FALSE]
      t_boot[b] <- compute_test_stat(boot_data, g_chain3)$statistic
    }
    p_val <- mean(t_boot >= t_obs)
    if (p_val < 0.05) rejections <- rejections + 1
  }
  rate <- rejections / n_sims
  cat(sprintf("    N=%d: %.3f\n", N, rate))
  results <- c(results, list(data.frame(
    dag = "Chain ($n=3$)", N = N, rate = rate, bootstrap = "Naive"
  )))
}

results_df <- do.call(rbind, results)

# --- Write LaTeX table ---
lines <- c(
  "& \\multicolumn{4}{c}{Sample Size $N$} & \\\\",
  "\\cmidrule(lr){2-5}",
  "DAG Type & 100 & 300 & 500 & 1000 & Bootstrap Type \\\\",
  "\\midrule"
)

# Constrained rows
for (dag_name in unique(results_df$dag[results_df$bootstrap == "Constrained"])) {
  sub <- results_df[results_df$dag == dag_name & results_df$bootstrap == "Constrained", ]
  sub <- sub[order(sub$N), ]
  lines <- c(lines, sprintf("%s & %s & %s & %s & %s & Constrained \\\\",
    dag_name,
    format_rate(sub$rate[1], 3),
    format_rate(sub$rate[2], 3),
    format_rate(sub$rate[3], 3),
    format_rate(sub$rate[4], 3)))
}

# Naive comparison
lines <- c(lines,
  "\\midrule",
  "\\multicolumn{6}{l}{\\textit{Comparison: Naive bootstrap (invalid)}} \\\\"
)
sub <- results_df[results_df$bootstrap == "Naive", ]
sub <- sub[order(sub$N), ]
lines <- c(lines, sprintf("Chain ($n=3$) & %s & %s & %s & %s & Naive \\\\",
  format_rate(sub$rate[1], 3),
  format_rate(sub$rate[2], 3),
  format_rate(sub$rate[3], 3),
  format_rate(sub$rate[4], 3)))

write_tex_table(lines, file.path(output_dir(), "table02_size.tex"), "lccccc")
cat("Experiment 1 complete.\n")
