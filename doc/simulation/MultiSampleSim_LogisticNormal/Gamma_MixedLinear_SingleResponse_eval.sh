#!/bin/bash

#SBATCH --job-name=Gamma_SingleResponse
#SBATCH --output=Gamma_SingleResponse.out
#SBATCH --error=Gamma_SingleResponse.err
#SBATCH -n 1
#SBATCH --time=1-
#SBATCH --mem=30g
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1
Rscript Gamma_MixedLinear_SingleResponse_eval.R
