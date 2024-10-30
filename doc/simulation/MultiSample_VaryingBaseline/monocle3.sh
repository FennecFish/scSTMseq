#!/bin/bash

#SBATCH --job-name=monocle3_Cancer_Batch
#SBATCH --output=monocle3_Cancer_Batch%A_%a.out
#SBATCH --error=monocle3_Cancer_Batch%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=5:00:00
#SBATCH --mem-per-cpu=10G
#SBATCH --array=1-40
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.2

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline/MultiSample_VaryingBaseline_nSample6_nCellType8_noBatch_CancerCell/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript monocle3.R "${FILES[$INDEX]}"







