# Data

## DREAM5 Gene Expression Data

The gene network analysis (Section 9.2) uses expression data from the DREAM5 Network Inference Challenge.

**Source:** Marbach et al. (2012), "Wisdom of crowds for robust gene network inference." *Nature Methods*, 9(8), 796-804.

**Subset used:** 3 genes from E. coli (lexA, recA, umuD) involved in the SOS DNA repair response.

**To obtain the data:**
1. Create an account at https://www.synapse.org
2. Download DREAM5 challenge data (E. coli expression matrix)
3. Extract columns for lexA, recA, umuD
4. Save as `dream5_subset.csv` in this directory

If this file is not present, the scripts will use simulated data with similar properties.

## IV Data (Angrist & Krueger 1991)

The IV analysis (Section 9.1) uses data from the AER R package (`SchoolingReturns` dataset). Install with:

```r
install.packages("AER")
```

If AER is not available, simulated IV-valid data is used.
