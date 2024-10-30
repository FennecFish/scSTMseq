#!/bin/bash

#SBATCH --job-name=scSTM_Pooled_S3_C5_noBatch_Stromal
#SBATCH --output=scSTM_Pooled_S3_C5_noBatch_Stromal%A_%a.out
#SBATCH --error=scSTM_Pooled_S3_C5_noBatch_Stromal%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=100G
#SBATCH --array=1-49
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SingleResponse/nSample3_nCellType5_noBatch_StromalCell/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*Null*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_Pooled_MultiSample.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







