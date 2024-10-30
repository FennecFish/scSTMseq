#!/bin/bash

#SBATCH --job-name=supp_V3_n12_time
#SBATCH --output=supp_V3_n12_time_%A_%a.out
#SBATCH --error=supp_V3_n12_time_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --time=1-
#SBATCH --mem=60G
#SBATCH --array=1-6
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

module load r/4.3.1
SIMS_DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims"
SCSTM_DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/scSTM_Content_TimeNoInteraction_Prevalence_Time"
# Find files that have 'sims_' in their names
SIMS_FILES=($(find "$SIMS_DIR" -maxdepth 1 -type f -name "sims*nsample12.rds"))

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

# Show all files in the FILES array
# for FILE in "${FILES[@]}"; do
#     echo "$FILE"
# done

# Display the count of files in the FILES array
# echo "There are ${#FILES[@]} files in total."

Rscript scSTM_Content_Time_Prevalence_Time.R "${FILES[$INDEX]}" "$SLURM_NTASKS"
