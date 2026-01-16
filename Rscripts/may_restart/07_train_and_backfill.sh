#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=05:00:00
#SBATCH --job-name=test_jan7_2026
#SBATCH --array=114,116,117,118,160,221,246,247,304,423,543,546,547,553,565,57,62,98%18
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla 07_train_and_backfill.R ${SLURM_ARRAY_TASK_ID}
