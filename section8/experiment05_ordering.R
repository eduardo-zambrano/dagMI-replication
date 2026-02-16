#!/usr/bin/env Rscript
# Experiment 5: Ordering selection strategies
# Diamond DAG n=4, N=1000 with 5 ordering strategies
# Output: output/table07_ordering.tex (Paper Table 7)

library(dagMI)
if (!exists("write_tex_table", mode = "function")) source("R/helpers.R")

cat("=== Experiment 5: Ordering Strategies ===\n")

set.seed(2026)
params <- sim_params(full_n_sims = 50, full_B = 100)
n_sims <- params$n_sims
B <- params$B
N <- 1000
beta <- 0.5

# Diamond n=4: X1 -> X2, X1 -> X3, X2 -> X4, X3 -> X4
A_diamond <- matrix(c(0,1,1,0, 0,0,0,1, 0,0,0,1, 0,0,0,0), 4, 4, byrow = TRUE)
g_diamond <- dag(A_diamond, nodes = c("X1","X2","X3","X4"))

# Alternative DAG for power: Hub X2 -> X1, X2 -> X3, X2 -> X4
A_hub <- matrix(c(0,0,0,0, 1,0,1,1, 0,0,0,0, 0,0,0,0), 4, 4, byrow = TRUE)
g_hub <- dag(A_hub, nodes = c("X1","X2","X3","X4"))

generate_diamond <- function(N, beta) {
  X1 <- rnorm(N)
  X2 <- beta * X1 + rnorm(N, sd = sqrt(1 - beta^2))
  X3 <- beta * X1 + rnorm(N, sd = sqrt(1 - beta^2))
  X4 <- beta/2 * (X2 + X3) + rnorm(N, sd = sqrt(1 - beta^2))
  m <- cbind(X1, X2, X3, X4)
  colnames(m) <- c("X1","X2","X3","X4")
  m
}

# Get all topological orderings for each DAG
all_orders_diamond <- topological_orders(g_diamond)
all_orders_hub <- topological_orders(g_hub)
cat(sprintf("  Diamond has %d topological orderings\n", length(all_orders_diamond)))
cat(sprintf("  Hub has %d topological orderings\n", length(all_orders_hub)))

# Helper: run MI test with a specific ordering (as a permutation vector)
run_test_ordering <- function(data_mat, g, ordering_perm, B, seed) {
  mi_test(data_mat, g, B = B, ordering = ordering_perm,
          verbose = FALSE, seed = seed)
}

# Helper: run MI test for all orderings and collect p-values
run_all_orderings <- function(data_mat, g, all_orders, B, seed) {
  p_vals <- numeric(length(all_orders))
  t_vals <- numeric(length(all_orders))
  for (k in seq_along(all_orders)) {
    res <- mi_test(data_mat, g, B = B, ordering = all_orders[[k]],
                   verbose = FALSE, seed = seed)
    p_vals[k] <- res$p_value
    t_vals[k] <- res$statistic
  }
  list(p_values = p_vals, statistics = t_vals)
}

# --- Strategy definitions ---
strategies <- list(
  list(name = "Single (arbitrary) ordering"),
  list(name = "Data-driven selection (Eq.\\ \\ref{eq:data-driven-ordering})"),
  list(name = "Min over all orderings"),
  list(name = "Bonferroni (all orderings)"),
  list(name = "Fisher combination")
)

# --- Run strategies 1 & 2 (single ordering each) ---
for (s_idx in 1:2) {
  strat <- strategies[[s_idx]]
  cat(sprintf("  %s: ", strat$name))
  ord_type <- if (s_idx == 1) "first" else "optimal"
  size_rej <- 0
  power_rej <- 0
  for (sim in 1:n_sims) {
    set.seed(2026 + sim)
    data_mat <- generate_diamond(N, beta)
    res_size <- mi_test(data_mat, g_diamond, B = B, ordering = ord_type,
                        verbose = FALSE, seed = sim)
    res_power <- mi_test(data_mat, g_hub, B = B, ordering = ord_type,
                         verbose = FALSE, seed = sim)
    if (res_size$decision == "reject") size_rej <- size_rej + 1
    if (res_power$decision == "reject") power_rej <- power_rej + 1
  }
  size_rate <- size_rej / n_sims
  power_rate <- power_rej / n_sims
  cat(sprintf("size=%.3f, power=%.2f\n", size_rate, power_rate))
  strategies[[s_idx]]$size <- size_rate
  strategies[[s_idx]]$power <- power_rate
}

# --- Run all-orderings once, reuse for strategies 3-5 ---
cat("  Computing all-orderings p-values (shared by Min/Bonferroni/Fisher)...\n")
rej_min_size <- 0; rej_min_power <- 0
rej_bonf_size <- 0; rej_bonf_power <- 0
rej_fisher_size <- 0; rej_fisher_power <- 0

for (sim in 1:n_sims) {
  set.seed(2026 + sim)
  data_mat <- generate_diamond(N, beta)

  # Compute p-values for all orderings (expensive — done ONCE per sim)
  res_all_size <- run_all_orderings(data_mat, g_diamond, all_orders_diamond, B, sim)
  res_all_power <- run_all_orderings(data_mat, g_hub, all_orders_hub, B, sim)

  # Strategy 3: Min over all orderings
  if (min(res_all_size$p_values) < 0.05) rej_min_size <- rej_min_size + 1
  if (min(res_all_power$p_values) < 0.05) rej_min_power <- rej_min_power + 1

  # Strategy 4: Bonferroni
  bonf_size <- bonferroni_orderings(res_all_size$p_values, alpha = 0.05)
  bonf_power <- bonferroni_orderings(res_all_power$p_values, alpha = 0.05)
  if (bonf_size$global_reject) rej_bonf_size <- rej_bonf_size + 1
  if (bonf_power$global_reject) rej_bonf_power <- rej_bonf_power + 1

  # Strategy 5: Fisher combination
  fisher_size <- fisher_combine(res_all_size$p_values)
  fisher_power <- fisher_combine(res_all_power$p_values)
  if (fisher_size$combined_p < 0.05) rej_fisher_size <- rej_fisher_size + 1
  if (fisher_power$combined_p < 0.05) rej_fisher_power <- rej_fisher_power + 1

  if (sim %% 10 == 0) cat(sprintf("    sim %d/%d\n", sim, n_sims))
}

strategies[[3]]$size <- rej_min_size / n_sims
strategies[[3]]$power <- rej_min_power / n_sims
cat(sprintf("  %s: size=%.3f, power=%.2f\n", strategies[[3]]$name,
    strategies[[3]]$size, strategies[[3]]$power))

strategies[[4]]$size <- rej_bonf_size / n_sims
strategies[[4]]$power <- rej_bonf_power / n_sims
cat(sprintf("  %s: size=%.3f, power=%.2f\n", strategies[[4]]$name,
    strategies[[4]]$size, strategies[[4]]$power))

strategies[[5]]$size <- rej_fisher_size / n_sims
strategies[[5]]$power <- rej_fisher_power / n_sims
cat(sprintf("  %s: size=%.3f, power=%.2f\n", strategies[[5]]$name,
    strategies[[5]]$size, strategies[[5]]$power))

# Collect results
results <- data.frame(
  strategy = sapply(strategies, `[[`, "name"),
  size = sapply(strategies, `[[`, "size"),
  power = sapply(strategies, `[[`, "power")
)

# --- Write LaTeX table ---
lines <- c(
  "Strategy & Size & Power \\\\",
  "\\midrule"
)
for (i in 1:nrow(results)) {
  lines <- c(lines, sprintf("%s & %s & %s \\\\",
    results$strategy[i],
    format_rate(results$size[i], 3),
    format_rate(results$power[i], 2)))
}

write_tex_table(lines, file.path(output_dir(), "table07_ordering.tex"), "lcc")
cat("Experiment 5 complete.\n")
