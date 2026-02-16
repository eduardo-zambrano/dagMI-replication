# dagMI Replication Package

Replication scripts for all data-driven tables in:

> Zambrano, E. (2026). "Testing Causal DAG Structures via Multilinear Inequalities."

## Quick Start

```bash
# 1. Install dependencies
Rscript install.R

# 2. Run all experiments (quick validation, ~25 min)
Rscript run_all.R --quick

# 3. Run all experiments (full publication quality, ~15-20 hours)
Rscript run_all.R

# 4. Copy generated tables to paper
cp output/*.tex ../DAG_Testing/tables/
```

## What This Reproduces

| Script | Output | Paper Table |
|--------|--------|-------------|
| `section7/benchmark_computation.R` | `table01_computation.tex` | Table 1 (Sec 7.3) |
| `section8/experiment01_size.R` | `table02_size.tex` | Table 2 (Sec 8.1) |
| `section8/experiment02_power.R` | `table03_power.tex` | Table 3 (Sec 8.2) |
| `section8/experiment03_comparison.R` | `table04_comparison.tex` | Table 4 (Sec 8.3) |
| `section8/experiment04_bandwidth.R` | `table06_bandwidth.tex` | Table 6 (Sec 8.4) |
| `section8/experiment05_ordering.R` | `table07_ordering.tex` | Table 7 (Sec 8.5) |
| `section8/experiment06_power_beta.R` | `table08_power_beta.tex` | Table 8 (Sec 8.6) |
| `section8/experiment07_large_dag.R` | `table09_large_dag.tex` | Table 9 (Sec 8.7) |
| `section8/experiment08_robustness.R` | `table10_robustness.tex` | Table 10 (Sec 8.8) |
| `section9/iv_returns.R` | `section9_iv_results.tex` | Section 9.1 (inline) |
| `section9/gene_network.R` | `table11_gene_network.tex` | Table 11 (Sec 9.2) |

Tables 5 and 12 are qualitative/non-data and remain hardcoded in the paper.

## Usage Options

```bash
Rscript run_all.R                  # Full run (publication quality)
Rscript run_all.R --quick          # Quick validation (~25 min)
Rscript run_all.R --section 8     # Only Section 8 experiments
Rscript run_all.R --experiment 3  # Only Experiment 3
```

## Dependencies

**Required:** R (>= 4.0.0), dagMI package

**Optional (for Experiment 3 method comparison):** pcalg, bnlearn

Run `Rscript install.R` to install all dependencies.

## Output Format

Each script produces a `.tex` file containing only the `\begin{tabular}...\end{tabular}` block. The paper imports these via `\input{tables/tableXX.tex}` and wraps them with `\begin{table}`, `\caption{}`, `\label{}`.

## License

MIT
