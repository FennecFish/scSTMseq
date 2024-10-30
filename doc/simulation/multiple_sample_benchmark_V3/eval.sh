#!/bin/bash

#SBATCH --job-name=eval_multi
#SBATCH --output=eval_multi.out
#SBATCH --error=eval_multi.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=30g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript eval.R
