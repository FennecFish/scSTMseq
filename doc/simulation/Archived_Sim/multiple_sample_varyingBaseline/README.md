# MultiSample_VaryingBaseline

## No Batch Effect

### Cancer Cell

#### nSample = 6:

1.  `nCell = 5`

2.  `nCell = 8`

    -   sims: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample6_nCellType8_noBatch_CancerCell/`
    -   scSTM with L1: 49118387
    -   scSTM with Pooled: 49110398
    -   scSTM with linear Regression

#### nSample 12:

1.  `nCell = 8`

    -   sims: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample12_nCellType8_noBatch_CancerCell/`
    -   scSTM with pooled: 49158906

#### 

### Stromal Cell

#### nSample = 6:

1.  `nCell = 5`

    -   sims: `"/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample5_nCellType5_noBatch_StromalCell/sims"`

    -   scSTM_Pooled: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample5_nCellType5_noBatch_StromalCell/scSTM_Pooled_noContent_Prevalence_TimeandResponse`

    -   scSTM_L1: 49110074, `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample5_nCellType5_noBatch_StromalCell/scSTM_L1_noContent_Prevalence_TimeandResponse`

    -   scSTM_linearRegression: 49160904

#### nSample 12:

## Batch Effect

### Cancer Cell

#### nSample = 6:

1.  `nCell = 5`

    -   sims: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample6_nCellType5_Batch_CancerCell/`
    -   scSTM

2.  `nCell = 8`

    -   sims: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample6_nCellType8_Batch_CancerCell/sims/`

    -   scSTM with Pooled Gamma: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample6_nCellType8_Batch_CancerCell/scSTM_Pooled_Content_Prevalence_TimeandResponse`

    -   scSTM with L1 Gamma: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample6_nCellType8_Batch_CancerCell/scSTM_L1_Content_Prevalence_TimeandResponse`

    -   scSTM with Linear Regression 49160942

#### nSample 12:

1.  `nCell = 5`
2.  `nCell = 8`
    -   sims: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample12_nCellType8_Batch_CancerCell/sims/`
    -   scSTM with Pooled Gamma

### Stromal Cell

#### nSample = 6:

#### nSample 12:

# Simulated Data sets

## With NO Batch Effect

In this category, `batch.rmEffect = TRUE` for all simulations

### No Heterogeneity Between Samples in Cells

-   **simulated data output path:** `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_noBatch_StromalCell/sims/`

-   **scSTM with `content = NULL`**: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_noBatch_StromalCell/scSTM_noContent_Prevalence_TimeandResponse`

-   **adjRandIndex:** `"res/adjRandIndex_MultiSample_VaryingBaseline_noBatch_StromalCell.csv"`

### Exist Heterogeneity Between Samples in Cancer Cell

This category uses `cancer_cell = 2`

-   **simulated data output path:** `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_noBatch_CancerCell/sims/`

-   **scSTM with `content = ~Sample`**: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_noBatch_CancerCell/scSTM_Content_Prevalence_TimeandResponse`

-   **scSTM with `content = NULL`**: `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_noBatch_CancerCell/scSTM_noContent_Prevalence_TimeandResponse`

-   **adjRandIndex**: `res/adjRandIndex_MultiSample_VaryingBaseline_noBatch_CancerCell.csv`

## With Batch Effect

In this category, `batch.rmEffect = FALSE` and `batch.facLoc <- 0.5` for all simulations

### No Heterogeneity Between Samples in Cells

**simulated data output path:**`/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_StromalCell/sims/` **scSTM with `content = ~Sample`:** `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_StromalCell/scSTM_noContent_Prevalence_TimeandResponse` **scSTM with `content = NULL`:** `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_StromalCell/scSTM_Content_Prevalence_TimeandResponse`

### Exist Heterogeneity Between Samples in Cancer Cell

This category uses `cancer_cell = 2`

**simulated data output path:** `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_CancerCell/sims/` **scSTM with `content = ~Sample`:** `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_CancerCell/scSTM_noContent_Prevalence_TimeandResponse` **scSTM with `content = NULL`:** `/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_CancerCell/scSTM_Content_Prevalence_TimeandResponse`
