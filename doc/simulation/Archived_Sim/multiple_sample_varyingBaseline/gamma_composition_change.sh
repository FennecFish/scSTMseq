#!/bin/bash

#SBATCH --job-name=gamma_cc
#SBATCH --output=gamma_cc.out
#SBATCH --error=gamma_cc.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=20g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript gamma_composition_change.R
