#!/bin/bash

#SBATCH --job-name=scSTM_Pooled_S6_C5_noBatch_Stromal
#SBATCH --output=scSTM_Pooled_S6_C5_noBatch_Stromal%A_%a.out
#SBATCH --error=scSTM_Pooled_S6_C5_noBatch_Stromal%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=70G
#SBATCH --array=1-5
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

SIMS_DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma/nSample4_nCellType5_noBatch_StromalCell/sims"
SCSTM_DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma/nSample4_nCellType5_noBatch_StromalCell/scSTM_Pooled_noContent_Prevalence_TimeandResponse"

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
INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_Pooled_S6.R "${FILES[$INDEX]}" "$SLURM_NTASKS"

# Show all files in the FILES array
# for FILE in "${FILES[@]}"; do
#     echo "$FILE"
# done

# Display the count of files in the FILES array
# echo "There are ${#FILES[@]} files in total."





