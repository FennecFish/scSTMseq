#!/bin/bash

#SBATCH --job-name=V3_n12_monocle3
#SBATCH --output=V3_n12_monocle3_%A_%a.out
#SBATCH --error=V3_n12_monocle3_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=5G
#SBATCH --array=1-90
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.2

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/multi_sample_benchmark_V3/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*nsample12.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript monocle3.R "${FILES[$INDEX]}"







