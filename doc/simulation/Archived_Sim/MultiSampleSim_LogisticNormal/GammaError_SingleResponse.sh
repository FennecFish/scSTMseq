#!/bin/bash

#SBATCH --job-name=GammaError_SingleResponse
#SBATCH --output=GammaError_SingleResponse%A_%a.out
#SBATCH --error=GammaError_SingleResponse%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=30G
#SBATCH --array=1
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

# DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SingleResponse/nSample3_nCellType5_noBatch_StromalCell/scSTM_LinearMixed_noSample_noContent_Prevalence_TimeandResponse"
# FILES=($(find "$DIR" -maxdepth 1 -type f -name "scSTM*.rds"))

# INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

# Rscript GammaError_SingleResponse.R "${FILES[$INDEX]}" 
Rscript GammaError_SingleResponse.R 







