# Data

## `gtex_jak_stat.csv` -- GTEx v8 Gene Expression (Section 9.2)

RNA-seq expression data for four genes in the JAK-STAT-interferon
signaling pathway, extracted from the GTEx v8 bulk tissue gene TPM
matrix (Whole Blood samples).

| Column | Gene | Role in pathway |
|--------|------|-----------------|
| JAK2 | Janus kinase 2 | Upstream kinase |
| STAT1 | Signal transducer and activator of transcription 1 | Transcription factor (hub) |
| IRF1 | Interferon regulatory factor 1 | STAT1 target |
| GBP1 | Guanylate-binding protein 1 | STAT1 target |

- **Samples:** 755 Whole Blood donors
- **Units:** log2(TPM + 1)
- **Source:** GTEx Portal (https://gtexportal.org/home/datasets), file
  `GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.gct.gz`

**Reference:**

> The GTEx Consortium (2020). "The GTEx Consortium Atlas of Genetic
> Regulatory Effects across Human Tissues." *Science*, 369(6509),
> 1318--1330. https://doi.org/10.1126/science.aaz1776

## Card (1995) -- Section 9.1

The IV analysis data are accessed programmatically via the `wooldridge`
R package (`data("card")`). No data file is stored in this directory.

- **Sample:** N = 3,010 men from the National Longitudinal Survey of
  Young Men (NLSYM)
- **Variables used:** `nearc4` (binary indicator for proximity to a
  four-year college), `educ` (years of completed schooling), `lwage`
  (log hourly wage)

**Reference:**

> Card, D. (1995). "Using Geographic Variation in College Proximity to
> Estimate the Return to Schooling." In *Aspects of Labour Market
> Behaviour: Essays in Honour of John Vanderkamp*, ed. L. N. Christofides,
> E. K. Grant, and R. Swidinsky, 201--222. University of Toronto Press.
