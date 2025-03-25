#!/bin/bash

#SBATCH --job-name=final_eval_multi
#SBATCH --output=final_eval_multi.out
#SBATCH --error=final_eval_multi.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=45g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript eval.R
