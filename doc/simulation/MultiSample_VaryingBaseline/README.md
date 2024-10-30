# Simulated Data sets

## With NO Batch Effect

In this category, `batch.rmEffect = TRUE` for all simulations

### No Heterogeneity Between Samples in Cells

-   **simulated data output path:** `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_noBatch_StromalCell/sims/`

-   **scSTM with `content = NULL`**: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_noBatch_StromalCell/scSTM_noContent_Prevalence_TimeandResponse`
- **adjRandIndex:**
`"res/adjRandIndex_MultiSample_VaryingBaseline_noBatch_StromalCell.csv"`

### Exist Heterogeneity Between Samples in Cancer Cell

This category uses `cancer_cell = 2`

-   **simulated data output path:** `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_noBatch_CancerCell/sims/`

-   **scSTM with `content = ~Sample`**: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_noBatch_CancerCell/scSTM_Content_Prevalence_TimeandResponse`

-   **scSTM with `content = NULL`**: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_noBatch_CancerCell/scSTM_noContent_Prevalence_TimeandResponse`

-   **adjRandIndex**: `res/adjRandIndex_MultiSample_VaryingBaseline_noBatch_CancerCell.csv`

## With Batch Effect

In this category, `batch.rmEffect = FALSE` and `batch.facLoc <- 0.5` for all simulations

### No Heterogeneity Between Samples in Cells

**simulated data output path:**`/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_StromalCell/sims/`
**scSTM with `content = ~Sample`:** `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_StromalCell/scSTM_noContent_Prevalence_TimeandResponse`
**scSTM with `content = NULL`:** 
`/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_StromalCell/scSTM_Content_Prevalence_TimeandResponse`

### Exist Heterogeneity Between Samples in Cancer Cell

This category uses `cancer_cell = 2`

**simulated data output path:** `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_CancerCell/sims/`
**scSTM with `content = ~Sample`:** `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_CancerCell/scSTM_noContent_Prevalence_TimeandResponse`
**scSTM with `content = NULL`:** 
`/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_CancerCell/scSTM_Content_Prevalence_TimeandResponse`

