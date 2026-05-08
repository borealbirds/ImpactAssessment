#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=4:00:00
#SBATCH --job-name=combine_obs_bf
#SBATCH --array=1-2
#SBATCH --mail-user=mannfred@ualberta.ca

# Phase 3: fill obs columns into bf-only density tables after 12A completes.
# Run after 12A_observed.sh finishes:
#   OBS_JOB_ID=$(sbatch --parsable 12A_observed.sh)
#   sbatch --dependency=afterok:$OBS_JOB_ID 12D_combine.sh
# Adjust --array to match the number of species.

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla 12D_combine.R
