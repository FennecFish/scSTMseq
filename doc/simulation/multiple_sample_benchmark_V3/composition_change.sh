#!/bin/bash

#SBATCH --job-name=V3_cc
#SBATCH --output=V3_cc.out
#SBATCH --error=V3_cc.err
#SBATCH -n 1
#SBATCH --time=2-
#SBATCH --mem=15g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript composition_change.R
