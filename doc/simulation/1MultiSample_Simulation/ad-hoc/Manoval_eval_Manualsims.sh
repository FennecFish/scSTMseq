#!/bin/bash

#SBATCH --job-name=Manova1000
#SBATCH --output=Manova1000%A_%a.out
#SBATCH --error=Manova1000%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=10G
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1
Rscript Manoval_eval_Manualsims.R





