#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=8:00:00
#SBATCH --job-name=obs_predict
#SBATCH --mail-user=mannfred@ualberta.ca
# submit: bash 12A_observed.sh SPECIES [SPECIES ...]   (e.g. OSFL GRSP)

if [ -z "${SLURM_JOB_ID:-}" ]; then
  set -euo pipefail
  [ $# -ge 1 ] || { echo "usage: bash 12A_observed.sh SPECIES [SPECIES ...]"; exit 1; }
  export SPECIES=$(IFS=,; echo "$*")
  J=$(sbatch --parsable --export=ALL --array=1-$# "$0")
  echo "OBS species=$SPECIES 12A=$J (1-$#)"
  exit 0
fi

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla 12A_observed.R
