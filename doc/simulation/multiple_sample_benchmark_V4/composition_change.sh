#!/bin/bash

#SBATCH --job-name=V4_cc
#SBATCH --output=V4_cc.out
#SBATCH --error=V4_cc.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=15g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript composition_change.R
