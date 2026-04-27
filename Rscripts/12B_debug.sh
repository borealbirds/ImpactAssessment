#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=02:00:00
#SBATCH --job-name=12B_debug
#SBATCH --output=/home/mannfred/scratch/impact_assessment/logs/12B_debug_%j.out
#SBATCH --error=/home/mannfred/scratch/impact_assessment/logs/12B_debug_%j.err

module load StdEnv/2023 gcc/12.3 gdal/3.9.1 udunits/2.2.28 r/4.4.0

ulimit -s unlimited

Rscript --vanilla 12B_debug.R