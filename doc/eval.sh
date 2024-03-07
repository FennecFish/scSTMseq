#!/bin/bash

#SBATCH --job-name=eval
#SBATCH --output=eval_%A_%a.out
#SBATCH --error=eval_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=15G
#SBATCH --array=1-2
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="../data"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

echo ${FILES[$INDEX]}
Rscript eval.R "${FILES[$INDEX]}"


