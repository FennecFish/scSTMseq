#!/bin/bash

#SBATCH --job-name=psi_LR
#SBATCH --output=psi_LR.out
#SBATCH --error=psi_LR.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript Psi_LinearRegression_eval.R
