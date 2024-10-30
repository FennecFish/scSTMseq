#!/bin/bash

#SBATCH --job-name=scSTM_LR_S6_C5_noBatch_Stromal
#SBATCH --output=scSTM_LR_S6_C5_noBatch_Stromal%A_%a.out
#SBATCH --error=scSTM_LR_S6_C5_noBatch_Stromal%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=40G
#SBATCH --array=1-40
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample6_nCellType5_noBatch_StromalCell/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_LinearRegression_S6_Content_SampleInteraction_Prevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







