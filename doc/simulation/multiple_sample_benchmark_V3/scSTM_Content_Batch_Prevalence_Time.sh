#!/bin/bash

#SBATCH --job-name=V3_n6_batch
#SBATCH --output=V3_n6_batch_%A_%a.out
#SBATCH --error=V3_n6_batch_%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=2-
#SBATCH --mem=50G
#SBATCH --array=1-90
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims_*_nsample6.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_Content_Batch_Prevalence_Time.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







