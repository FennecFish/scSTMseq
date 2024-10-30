#!/bin/bash

#SBATCH --job-name=VB_noBatch_Stromal
#SBATCH --output=VB_noBatch_Stromal_%A_%a.out
#SBATCH --error=VB_noBatch_Stromal_%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=30G
#SBATCH --array=1-40
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_noBatch_StromalCell/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_noBatch_Stromal.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







