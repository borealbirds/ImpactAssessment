#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --job-name=covmosaic
#SBATCH --mail-user=mannfred@ualberta.ca
# submit: bash 06_build_covariate_stacks.sh YEAR [YEAR ...]   (e.g. 1995 2000)

if [ -z "${SLURM_JOB_ID:-}" ]; then
  set -euo pipefail
  [ $# -ge 1 ] || { echo "usage: bash 06_build_covariate_stacks.sh YEAR [YEAR ...]"; exit 1; }
  export YEARS=$(IFS=,; echo "$*")
  J=$(sbatch --parsable --export=ALL --array=1-$# "$0")
  echo "COVMOSAIC years=$YEARS 06=$J (1-$#)"
  exit 0
fi

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

export TMPDIR=${SLURM_TMPDIR:-$TMPDIR}
export YEARS=$(echo "$YEARS" | cut -d, -f"${SLURM_ARRAY_TASK_ID}")
Rscript --vanilla 06_build_covariate_stacks.R
