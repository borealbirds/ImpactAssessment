#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=128G
#SBATCH --time=06:00:00
#SBATCH --job-name=2020_run_backfill_larger
#SBATCH --array=19,24,27,30,48,53,54,55,56,57,58,61,62,63,64,77,82,83,91,92,98,100,102,103,104,105,106,107,108,122,179,184,190,203,220,228,235,239,260,272,273,276,278,296,314,317,318,319,324,339,347,352,360,376,377,379,386,397,405,427,438,461,467,472,486,488,512,516,530,538,540,543,546,547,553,558,559,560,565,572,586,588,596,599,603,608,621,629,630,637,638,642,651,662,666,672,673
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla 07_train_and_backfill.R ${SLURM_ARRAY_TASK_ID}
