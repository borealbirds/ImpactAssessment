#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --job-name=change_maps
#SBATCH --mail-user=mannfred@ualberta.ca
# submit: sbatch 18_map_population_change.sh   (every species in 14B's tables; SPECIES=CAWA,OVEN sbatch ... for a subset)
# year:   YEARS=1995 sbatch 18_map_population_change.sh   (one year; default 2020)

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla 18_map_population_change.R
