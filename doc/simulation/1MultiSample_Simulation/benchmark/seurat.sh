#!/bin/bash

#SBATCH --job-name=seurat
#SBATCH --output=seurat_%A_%a.out
#SBATCH --error=seurat_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --array=11-988
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"
FILES=($(find "$DIR" -type f -path "*_noBatch_*/1000sims/*.rds"))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript seurat.R "$PARENT_DIR" "$FILE"







