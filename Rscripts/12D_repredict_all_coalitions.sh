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
  BIG_BCRS="can11"
  MEM_BIG=64G
  export SPECIES=$(IFS=,; echo "$*")

  complement() {
    echo "$2" | awk -v n="$1" 'BEGIN { p = 1; out = "" }
      { if ($1 > p) out = out (out == "" ? "" : ",") ($1 - 1 > p ? p "-" ($1 - 1) : p); p = $1 + 1 }
      END { if (p <= n) out = out (out == "" ? "" : ",") (n > p ? p "-" n : p); print out }'
  }

  N=0
  BIG=""
  for sp in "$@"; do
    [ -d "$BOOT/$sp" ] || { echo "no bootstrap models for $sp in $BOOT"; exit 1; }
    for b in $(ls "$BOOT/$sp" | grep -E 'can.*\.Rdata$' | sed -E "s/^${sp}_(.*)\.Rdata$/\1/"); do
      if [ -z "${TEST_BCR:-}" ] || [[ ",$TEST_BCR," == *",$b,"* ]]; then
        N=$((N + 1))
        if [[ " $BIG_BCRS " == *" $b "* ]]; then BIG="$BIG $N"; fi
      fi
    done
  done
  [ $N -ge 1 ] || { echo "no can* models for SPECIES=$SPECIES${TEST_BCR:+ in TEST_BCR=$TEST_BCR}"; exit 1; }
  if [ -n "$BIG" ]; then SMALL=$(complement $N "$(printf '%s\n' $BIG)"); else SMALL="1-$N"; fi
  BIG=$(echo $BIG | tr ' ' ',')

  JS=""
  JB=""
  if [ -n "$SMALL" ]; then JS=$(sbatch --parsable --export=ALL --array=$SMALL "$0"); fi
  if [ -n "$BIG" ]; then JB=$(sbatch --parsable --export=ALL --mem=$MEM_BIG --array=$BIG "$0"); fi
  H=$(sbatch --parsable --export=ALL --dependency=afterany${JS:+:$JS}${JB:+:$JB} 12H_merge_bcr_tables.sh)
  echo "COAL species=$SPECIES${TEST_BCR:+ TEST_BCR=$TEST_BCR}${TEST_N_BOOT:+ TEST_N_BOOT=$TEST_N_BOOT}${DROP_Q99:+ DROP_Q99=$DROP_Q99}${DENS_MIN:+ DENS_MIN=$DENS_MIN} 12D=${JS:-none} (48G: ${SMALL:-none}) + ${JB:-none} ($MEM_BIG: ${BIG:-none}) 12H=$H"
  exit 0
fi

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla 12D_repredict_all_coalitions.R
