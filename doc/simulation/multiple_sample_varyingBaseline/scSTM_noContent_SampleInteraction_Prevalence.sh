#!/bin/bash

#SBATCH --job-name=scSTM_L1_noBatch_Stromal
#SBATCH --output=scSTM_L1_noBatch_Stromal%a.out
#SBATCH --error=scSTM_L1_noBatch_Stromal%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=40G
#SBATCH --array=1-40
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_varyingBaseline/MultiSample_VaryingBaseline_noBatch_StromalCell/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_noContent_SampleInteraction_Prevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







