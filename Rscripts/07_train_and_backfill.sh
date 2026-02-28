#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --job-name=2020_run_backfill
#SBATCH --array=1-674%10
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla 07_train_and_backfill.R ${SLURM_ARRAY_TASK_ID}

# these settings worked for the vast majority of subbasins
# I ran:
# grep -L "writing .* layers to .*/subbasin_[0-9]\+" Y2020_*.log | \
# xargs -r grep -L "^done$"
# to identify subbasins that were OOM killed or timed out and re-ran with 07_train_and_backfill_larger.sh