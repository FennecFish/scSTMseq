#!/bin/bash

#SBATCH --job-name=final_n36_multi_C_P
#SBATCH --output=final_n36_multi_C_P_%A_%a.out
#SBATCH --error=final_n36_multi_C_P_%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=2-
#SBATCH --mem=250G
#SBATCH --array=2-160
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*nsample36.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_Content_SampleandTime_Prevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







