#!/bin/bash

#SBATCH --job-name=theta_LR
#SBATCH --output=theta_LR.out
#SBATCH --error=theta_LR.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=30g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript Theta_LinearRegression_eval.R
