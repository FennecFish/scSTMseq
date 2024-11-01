#!/bin/bash
#SBATCH --job-name=MeanTimeCalculation_linearMixed
#SBATCH --output=MeanTimeCalculation_linearMixed%j.out
#SBATCH --error=MeanTimeCalculation_linearMixed%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --mem=1G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

# Job ID
ARRAY_JOB_ID=52498463

# output file
OUTPUT_FILE="mean_time_${ARRAY_JOB_ID}.txt"
echo -e "Array Job ID\tComputing Time" > "$OUTPUT_FILE"

# calculate the total time and count of completed tasks
total_time=0
count=0

# elapsed times for COMPLETED tasks in the array job
for time in $(sacct -j "$ARRAY_JOB_ID" --state=COMPLETED -o Elapsed -n | awk '{print $1}'); do
    # Convert HH:MM:SS to seconds
    IFS=: read -r hours minutes seconds <<< "$time"
    elapsed_seconds=$((10#$hours * 3600 + 10#$minutes * 60 + 10#$seconds))
    
    # Add to total time and increment count
    total_time=$((total_time + elapsed_seconds))
    count=$((count + 1))
done

# Calculate the mean time if there are any completed tasks
if (( count > 0 )); then
    mean_seconds=$((total_time / count))
    # Convert mean time back to HH:MM:SS format
    mean_time_formatted=$(printf "%02d:%02d:%02d" $((mean_seconds / 3600)) $(( (mean_seconds % 3600) / 60)) $((mean_seconds % 60)))
    echo -e "${ARRAY_JOB_ID}\t${mean_time_formatted}" >> "$OUTPUT_FILE"
    echo "Mean time for COMPLETED tasks written to $OUTPUT_FILE"
else
    echo -e "${ARRAY_JOB_ID}\tNo completed tasks found" >> "$OUTPUT_FILE"
    echo "No completed tasks found for array job $ARRAY_JOB_ID."
fi

