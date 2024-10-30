#!/bin/bash

#SBATCH --job-name=eval_de
#SBATCH --output=eval_de.out
#SBATCH --error=eval_de.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=13g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript cellType_Topic_Match.R
