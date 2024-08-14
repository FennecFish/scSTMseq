#!/bin/bash

#SBATCH --job-name=null_scSTM
#SBATCH --output=null_scSTM_%A_%a.out
#SBATCH --error=null_scSTM_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-606
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/fig1/V3_single/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript null_scSTM_single_filterGenes_noContent_run.R "${FILES[$INDEX]}"







