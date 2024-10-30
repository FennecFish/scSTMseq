#!/bin/bash

#SBATCH --job-name=n24_seurat
#SBATCH --output=n24_seurat_%A_%a.out
#SBATCH --error=n24_seurat_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=8G
#SBATCH --array=2-160
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_final/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*nsample24.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript seurat.R "${FILES[$INDEX]}"






