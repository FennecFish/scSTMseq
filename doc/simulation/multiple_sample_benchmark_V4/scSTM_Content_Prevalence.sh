#!/bin/bash

#SBATCH --job-name=V4_n12_multi_C_P
#SBATCH --output=V4_n12_multi_C_P_%A_%a.out
#SBATCH --error=V4_n12_multi_C_P_%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=35G
#SBATCH --array=1-60
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V4/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*nsample12.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_Content_Prevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







