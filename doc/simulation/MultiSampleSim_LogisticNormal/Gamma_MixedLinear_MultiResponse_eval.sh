#!/bin/bash

#SBATCH --job-name=Gamma_MixedLinear
#SBATCH --output=Gamma_MixedLinear.out
#SBATCH --error=Gamma_MixedLinear.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=80g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript Gamma_MixedLinear_MultiResponse_eval.R
