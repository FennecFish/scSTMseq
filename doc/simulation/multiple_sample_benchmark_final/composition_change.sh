#!/bin/bash

#SBATCH --job-name=final_cc
#SBATCH --output=final_cc.out
#SBATCH --error=final_cc.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript composition_change.R
