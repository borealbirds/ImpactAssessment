#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=48G
#SBATCH --time=03:00:00
#SBATCH --job-name=coal_bcr
#SBATCH --mail-user=mannfred@ualberta.ca
# submit: bash 12D_repredict_all_coalitions.sh SPECIES [SPECIES ...]   (e.g. OSFL GRSP)
# smoke:  TEST_BCR=can60 TEST_N_BOOT=2 bash 12D_repredict_all_coalitions.sh SPECIES [SPECIES ...]

if [ -z "${SLURM_JOB_ID:-}" ]; then
  set -euo pipefail
  [ $# -ge 1 ] || { echo "usage: bash 12D_repredict_all_coalitions.sh SPECIES [SPECIES ...]"; exit 1; }
  BOOT=/home/mannfred/projects/def-ecknight/NationalModels/output/06_bootstraps
  export SPECIES=$(IFS=,; echo "$*")
  N=0
  for sp in "$@"; do
    [ -d "$BOOT/$sp" ] || { echo "no bootstrap models for $sp in $BOOT"; exit 1; }
    for b in $(ls "$BOOT/$sp" | grep -E 'can.*\.Rdata$' | sed -E "s/^${sp}_(.*)\.Rdata$/\1/"); do
      if [ -z "${TEST_BCR:-}" ] || [[ ",$TEST_BCR," == *",$b,"* ]]; then N=$((N + 1)); fi
    done
  done
  [ $N -ge 1 ] || { echo "no can* models for SPECIES=$SPECIES${TEST_BCR:+ in TEST_BCR=$TEST_BCR}"; exit 1; }
  J=$(sbatch --parsable --export=ALL --array=1-$N "$0")
  H=$(sbatch --parsable --export=ALL --dependency=afterany:$J 12H_merge_bcr_tables.sh)
  echo "COAL species=$SPECIES${TEST_BCR:+ TEST_BCR=$TEST_BCR}${TEST_N_BOOT:+ TEST_N_BOOT=$TEST_N_BOOT}${DROP_Q99:+ DROP_Q99=$DROP_Q99} 12D=$J (1-$N) 12H=$H"
  exit 0
fi

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla 12D_repredict_all_coalitions.R
