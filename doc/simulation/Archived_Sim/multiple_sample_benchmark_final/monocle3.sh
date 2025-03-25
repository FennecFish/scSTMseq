#!/bin/bash

#SBATCH --job-name=n48_monocle3
#SBATCH --output=n48_monocle3_%A_%a.out
#SBATCH --error=n48_monocle3_%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2:00:00
#SBATCH --mem-per-cpu=12G
#SBATCH --array=2-160
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.2

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_noBatch_CancerCell/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*nsample48.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript monocle3.R "${FILES[$INDEX]}"







