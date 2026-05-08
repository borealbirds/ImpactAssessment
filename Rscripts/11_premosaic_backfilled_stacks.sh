#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=24:00:00
#SBATCH --job-name=premosaic
#SBATCH --array=1-19
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla 11_premosaic_backfilled_stacks.R

# some BCRs needed more memory:
# so we used --array=1,2,3,5,6,7,8,10,11,12,13,14,15,16,17 --mem=256G 11_premosaic_backfilled_stacks.sh

# some needed more memory yet:
# so we used --array=3,6,11,12,16,17 --mem=512G 11_premosaic_backfilled_stacks.sh
