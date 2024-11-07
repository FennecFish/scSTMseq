#!/bin/bash

#SBATCH --job-name=supp_scSTM_Pooled_Time
#SBATCH --output=supp_scSTM_Pooled_Time%A_%a.out
#SBATCH --error=supp_scSTM_Pooled_Time%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=60G
#SBATCH --array=11-55
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

module load r/4.3.1

SIMS_DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/nSample20_nCellType5_noBatch_StromalCell/Manualsims"
SCSTM_DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/nSample20_nCellType5_noBatch_StromalCell/scSTM_Manualsims_Pooled_noContent_Prevalence_Time"

# Find files that have 'sims_' in their names
SIMS_FILES=($(find "$SIMS_DIR" -maxdepth 1 -type f -name "sims*.rds"))

# Filter out files that already have a corresponding 'scSTM_' file
FILES=()
for SIM_FILE in "${SIMS_FILES[@]}"; do
    BASE_NAME=$(basename "$SIM_FILE")
    SCSTM_FILE="${SCSTM_DIR}/${BASE_NAME/sims_/scSTM_}"
    if [ ! -f "$SCSTM_FILE" ]; then
        FILES+=("$SIM_FILE")
    fi
done

# Show all files in the FILES array
# for FILE in "${FILES[@]}"; do
#     echo "$FILE"
# done

# Display the count of files in the FILES array
# echo "There are ${#FILES[@]} files in total."


# Calculate the index for the SLURM array
INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

# Get the file for the current array task
FILE="${FILES[$INDEX]}"

# # Construct the full file path
# FILE_PATH="$SIMS_DIR/$FILE"

# Extract the parent directory of the current file
PARENT_DIR=$(dirname "$(dirname "$FILE")")

Rscript scSTM_noSampleVariation.R "$PARENT_DIR" "$FILE" "$SLURM_NTASKS"

