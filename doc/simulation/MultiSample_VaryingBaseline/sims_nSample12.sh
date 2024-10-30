#!/bin/bash

#SBATCH --job-name=sim_S12
#SBATCH --output=sim_S12%A_%a.out
#SBATCH --error=sim_S12%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem-per-cpu=30G
#SBATCH --array=1-10
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

INDEX=$(($SLURM_ARRAY_TASK_ID - 1))

Rscript sims_nSample12.R $INDEX






