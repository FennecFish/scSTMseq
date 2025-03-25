#!/bin/bash

#SBATCH --job-name=adjRandIndex
#SBATCH --output=adjRandIndex.out
#SBATCH --error=adjRandIndex.err
#SBATCH -n 1
#SBATCH --time=8:00:00
#SBATCH --mem=90g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript eval.R
