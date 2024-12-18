#!/bin/bash

#SBATCH --job-name=Manova1000
#SBATCH --output=Manova%A_%a.out
#SBATCH --error=Manova%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem=5G
#SBATCH --array=1-988
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/"
STM_DIR="1000scSTM_Pooled_noContent_Prevalence_Time"
FILES=($(find "$DIR" -type f -path "*/$STM_DIR/*Null*.rds"))

# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript ManovaRM_SingleResponse.R "$PARENT_DIR" "$STM_DIR" "$FILE" 
# Rscript ManovaRM_SingleResponse.R







