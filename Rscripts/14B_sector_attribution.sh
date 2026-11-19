#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --job-name=shapley_tables
#SBATCH --mail-user=mannfred@ualberta.ca
# submit: bash 14B_sector_attribution.sh [12H job id]   (14B + 14C, then 18; SPECIES=CAWA,OVEN limits 18's maps; YEARS=1995,2000 for other years)

if [ -z "${SLURM_JOB_ID:-}" ]; then
  set -euo pipefail
  for y in $(echo "${YEARS:-2020}" | tr ',' ' '); do
    export YEARS=$y
    J=$(sbatch --parsable --export=ALL ${1:+--dependency=afterok:$1} "$0")
    M=$(sbatch --parsable --export=ALL --dependency=afterok:$J 18_map_population_change.sh)
    echo "REPORT year=$y${SPECIES:+ SPECIES=$SPECIES} 14B=$J${1:+ (after $1)} 18=$M"
  done
  exit 0
fi

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla 14B_sector_attribution.R && Rscript --vanilla 14C_results_tables.R
