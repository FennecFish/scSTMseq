#!/bin/bash

#SBATCH --job-name=Manova_Theta_SingleResponse_eval
#SBATCH --output=Manova_Theta_SingleResponse_eval%A_%a.out
#SBATCH --error=Manova_Theta_SingleResponse_eval%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=6G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1
Rscript Manova_eval.R





