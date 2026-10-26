#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=pop_samples
#SBATCH --mail-user=mannfred@ualberta.ca
# submit: sbatch 19A_population_samples.sh   (1995-2020; YEARS=2015,2020 sbatch ... for others)

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla 19A_population_samples.R
