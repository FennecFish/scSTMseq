#!/bin/bash

#SBATCH --job-name=V3_fastTopics
#SBATCH --output=V3_fastTopics_%A_%a.out
#SBATCH --error=V3_fastTopics_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --array=1-270
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript fastTopics.R "${FILES[$INDEX]}"





