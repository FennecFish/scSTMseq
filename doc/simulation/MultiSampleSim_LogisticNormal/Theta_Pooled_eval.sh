#!/bin/bash

#SBATCH --job-name=theta_Pooled
#SBATCH --output=theta_Pooled.out
#SBATCH --error=theta_Pooled.err
#SBATCH -n 1
#SBATCH --time=5:00:00
#SBATCH --mem=30g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript Theta_Pooled_eval.R
