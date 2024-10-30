#!/bin/bash

#SBATCH --job-name=seurat_Cancer_Batch
#SBATCH --output=seurat_Cancer_Batch%A_%a.out
#SBATCH --error=seurat_Cancer_Batch%A_%a.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=3:00:00
#SBATCH --mem-per-cpu=3G
#SBATCH --array=1-44
#SBATCH --mail-type=all
#SBATCH --mail-user=euphyw@live.unc.edu



module load r/4.3.1

DIR="/work/users/e/u/euphyw/scLDAseq/data/simulation/MultiSample_VaryingBaseline_Batch0.5_CancerCell/sims"
FILES=($(find "$DIR" -maxdepth 1 -type f -name "sims*.rds"))

INDEX=$(($SLURM_ARRAY_TASK_ID - 1)) # Calculate array index

Rscript seurat.R "${FILES[$INDEX]}"





