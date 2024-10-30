#!/bin/bash

#SBATCH --job-name=final_n36_multi_C_noP
#SBATCH --output=final_n36_multi_C_noP_%A_%a.out
#SBATCH --error=final_n36_multi_C_noP_%A_%a.err
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=300G
#SBATCH --array=1
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*nsample36.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_Content_SampleandTime_noPrevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







