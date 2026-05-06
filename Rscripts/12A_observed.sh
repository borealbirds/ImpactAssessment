#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=2:00:00
#SBATCH --job-name=obs_predict
#SBATCH --array=1-2
#SBATCH --mail-user=mannfred@ualberta.ca

# Phase 1: canonical observed-landscape predictions.
# Run this BEFORE 12B_repredict_birds.sh (coalition runs).
# Adjust --array to match the number of species.

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla 12A_observed.R
