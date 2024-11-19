#!/bin/bash

#SBATCH --job-name=supp_
#SBATCH --output=supp_scSTM_Pooled_Time%A_%a.out
#SBATCH --error=supp_scSTM_Pooled_Time%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=70G
#SBATCH --array=1-4
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

module load r/4.3.1

S_DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/*/seurat"
P_DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/1MultiSample/SingleResponse/*/propeller"

# Find files that have 'sims_' in their names
S_FILES=($(find "$S_DIR" -maxdepth 1 -type f -name "*/sims/*.rds"))

# Filter out files that already have a corresponding 'scSTM_' file
FILES=()
for SIM_FILE in $(find $S_FILES -type f -name "*.rds"); do
    # Extract the base name of the file
    BASE_NAME=$(basename "$SIM_FILE")
    
    # Identify the corresponding scSTM directory path by substituting 'sims' with 'scSTM_Pooled_noContent_Prevalence_Time'
    DIRECTORY=$(dirname "$SIM_FILE" | sed "s|sims|scSTM_Pooled_noContent_Prevalence_Time|")

    # Construct the corresponding scSTM file name
    SCSTM_FILE="${DIRECTORY}/${BASE_NAME/sims_/scSTM_}"
    
    # Check if the corresponding scSTM file does not exist
    if [ ! -f "$SCSTM_FILE" ]; then
        FILES+=("$SIM_FILE")
    fi
done

# # Show all files in the FILES array
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

Rscript scSTM_Pooled_SingleResponse.R "$PARENT_DIR" "$FILE" "$SLURM_NTASKS"


