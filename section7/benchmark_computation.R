#!/usr/bin/env Rscript
# Section 7.3: Computational benchmarks
# Times compute_qn() and full bootstrap test for various n and N
# Output: output/table01_computation.tex (Paper Table 1)

library(dagMI)
if (!exists("write_tex_table", mode = "function")) source("R/helpers.R")

cat("=== Section 7: Computational Benchmarks ===\n")

set.seed(2026)

# In quick mode: fewer/smaller configs, shorter timeout
if (is_quick_mode()) {
  node_sizes <- c(3, 5)
  sample_sizes <- c(500, 1000)
  B_bench <- 30
  timeout_sec <- 120  # 2 min timeout
} else {
  node_sizes <- c(3, 5, 8, 10, 15, 20)
  sample_sizes <- c(500, 1000, 5000)
  B_bench <- 100
  timeout_sec <- 3600  # 1 hour timeout per config
}

beta <- 0.5

# Build a chain DAG of size n
make_chain_dag <- function(n) {
  A <- matrix(0, n, n)
  for (i in 1:(n - 1)) A[i, i + 1] <- 1
  nodes <- paste0("X", 1:n)
  dag(A, nodes = nodes)
}

# Generate data from chain DAG
generate_chain <- function(N, n, beta) {
  X <- matrix(0, N, n)
  X[, 1] <- rnorm(N)
  for (j in 2:n) {
    X[, j] <- beta * X[, j - 1] + rnorm(N, sd = sqrt(1 - beta^2))
  }
  colnames(X) <- paste0("X", 1:n)
  X
}

# Results matrices
qn_times <- matrix(NA, length(node_sizes), length(sample_sizes))
rownames(qn_times) <- node_sizes
colnames(qn_times) <- sample_sizes

boot_times <- matrix(NA, length(node_sizes), length(sample_sizes))
rownames(boot_times) <- node_sizes
colnames(boot_times) <- sample_sizes

for (i in seq_along(node_sizes)) {
  n <- node_sizes[i]
  g <- make_chain_dag(n)

  for (j in seq_along(sample_sizes)) {
    N <- sample_sizes[j]
    cat(sprintf("  n=%d, N=%d: ", n, N))

    data_mat <- generate_chain(N, n, beta)

    # Time Q_n computation only
    t0 <- proc.time()[3]
    tryCatch({
      qn_res <- compute_qn(data_mat, g)
      qn_time <- proc.time()[3] - t0
      if (qn_time > timeout_sec) {
        cat("Q_n timeout\n")
        qn_times[i, j] <- NA
        boot_times[i, j] <- NA
        next
      }
      qn_times[i, j] <- qn_time
      cat(sprintf("Q_n=%.1fs, ", qn_time))
    }, error = function(e) {
      cat(sprintf("Q_n error: %s\n", e$message))
      qn_times[i, j] <<- NA
    })

    # Time full bootstrap test (with polling-based timeout)
    t0 <- proc.time()[3]
    tryCatch({
      job <- parallel::mcparallel(
        mi_test(data_mat, g, B = B_bench, ordering = "first",
                verbose = FALSE, seed = 1))
      result <- NULL
      while ((proc.time()[3] - t0) < timeout_sec) {
        res_list <- parallel::mccollect(job, wait = FALSE)
        if (!is.null(res_list) && length(res_list) > 0) {
          result <- res_list[[1]]
          break
        }
        Sys.sleep(5)
      }
      boot_time <- proc.time()[3] - t0
      if (is.null(result) || inherits(result, "try-error")) {
        try(tools::pskill(job$pid, signal = 9L), silent = TRUE)
        parallel::mccollect(job, wait = FALSE)  # clean up
        cat("bootstrap timeout\n")
        boot_times[i, j] <- NA
        next
      }
      boot_times[i, j] <- boot_time
      cat(sprintf("full=%.1fs\n", boot_time))
    }, error = function(e) {
      cat(sprintf("bootstrap error: %s\n", e$message))
      boot_times[i, j] <<- NA
    })
  }
}

# --- Format time for display ---
fmt_time <- function(t) {
  if (is.na(t)) return("---")
  if (t < 1) return(format_num(t, 2))
  if (t < 100) return(format_num(t, 1))
  return(as.character(round(t)))
}

# --- Write LaTeX table ---
lines <- c(
  paste0("& \\multicolumn{", length(sample_sizes), "}{c}{$\\widehat{Q}_n^{\\mathcal{G}}$ only}",
         " & \\multicolumn{", length(sample_sizes), "}{c}{Full test with bootstrap} \\\\"),
  paste0("\\cmidrule(lr){2-", 1 + length(sample_sizes), "}",
         " \\cmidrule(lr){", 2 + length(sample_sizes), "-",
         1 + 2 * length(sample_sizes), "}"),
  paste0("$n$ & ",
         paste(sprintf("$N=%s$", sample_sizes), collapse = " & "),
         " & ",
         paste(sprintf("$N=%s$", sample_sizes), collapse = " & "),
         " \\\\"),
  "\\midrule"
)

for (i in seq_along(node_sizes)) {
  n <- node_sizes[i]
  qn_vals <- sapply(1:length(sample_sizes), function(j) fmt_time(qn_times[i, j]))
  boot_vals <- sapply(1:length(sample_sizes), function(j) fmt_time(boot_times[i, j]))
  lines <- c(lines, sprintf("%d & %s & %s \\\\",
    n,
    paste(qn_vals, collapse = " & "),
    paste(boot_vals, collapse = " & ")))
}

ncols <- 1 + 2 * length(sample_sizes)
align <- paste0("l", paste(rep("c", ncols - 1), collapse = ""))
write_tex_table(lines, file.path(output_dir(), "table01_computation.tex"), align)
cat("Section 7 benchmarks complete.\n")
