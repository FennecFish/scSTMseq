#!/bin/bash

#SBATCH --job-name=V4_eval_multi
#SBATCH --output=V4_eval_multi.out
#SBATCH --error=V4_eval_multi.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=30g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript eval.R
