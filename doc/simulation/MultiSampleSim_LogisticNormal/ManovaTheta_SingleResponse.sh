#!/bin/bash

#SBATCH --job-name=Manova_Theta_SingleResponse_eval
#SBATCH --output=Manova_Theta_SingleResponse_eval%A_%a.out
#SBATCH --error=Manova_Theta_SingleResponse_eval%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=8:00:00
#SBATCH --mem=30G
#SBATCH --array=1
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

# DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_LogisticNormal/NormalGamma_SingleResponse/nSample20_nCellType5_noBatch_CancerCell/scSTM_LinearMixed_noContent_Prevalence_Time"
# FILES=($(find "$DIR" -maxdepth 1 -type f -name "scSTM*.rds"))
# 
# INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index
# 
# Rscript Manova_Theta_SingleResponse_eval.R "${FILES[$INDEX]}"
Rscript Manova_Theta_SingleResponse_eval.R







