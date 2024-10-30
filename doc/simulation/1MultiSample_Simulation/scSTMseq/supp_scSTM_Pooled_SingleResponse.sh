#!/bin/bash

#SBATCH --job-name=supp_scSTM_Pooled_Time
#SBATCH --output=supp_scSTM_Pooled_Time%A_%a.out
#SBATCH --error=supp_scSTM_Pooled_Time%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=60G
#SBATCH --array=1-3
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/nSample20_nCellType15_noBatch_StromalCell/sims"
FILES=("sims_1729998879_HighVar1.rds" "sims_1729998867_HighVar1.rds" "sims_1729998879_HighVar0.6.rds")

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Construct the full file path
FILE_PATH="$DIR/$FILE"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE_PATH")")

Rscript scSTM_Pooled_SingleResponse.R "$PARENT_DIR" "$FILE_PATH" "$SLURM_NTASKS"

