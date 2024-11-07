#!/bin/bash

#SBATCH --job-name=adjRandIndex
#SBATCH --output=adjRandIndex.out
#SBATCH --error=adjRandIndex.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=2-
#SBATCH --mem=50G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu

module load r/4.3.1
Rscript adjustedRandIndex.R 

