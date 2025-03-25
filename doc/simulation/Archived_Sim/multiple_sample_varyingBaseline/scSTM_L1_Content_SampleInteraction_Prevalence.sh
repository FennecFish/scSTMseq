#!/bin/bash

#SBATCH --job-name=scSTM_L1_Content_Prevalence
#SBATCH --output=scSTM_L1_Content_Prevalence%A_%a.out
#SBATCH --error=scSTM_L1_Content_Prevalence%A_%a.err
#SBATCH --ntasks=10
#SBATCH --cpus-per-task=1
#SBATCH --time=1-
#SBATCH --mem=200G
#SBATCH --array=1
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu


module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_nSample12_nCellType5_Batch_CancerCell/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript scSTM_L1_Content_SampleInteraction_Prevalence.R "${FILES[$INDEX]}" "$SLURM_NTASKS"







