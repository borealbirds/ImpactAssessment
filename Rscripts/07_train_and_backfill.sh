#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --job-name=backfill
#SBATCH --array=1-674%30
#SBATCH --mail-user=mannfred@ualberta.ca
# submit: bash 07_train_and_backfill.sh YEAR [YEAR ...]   (e.g. 2015 2020)
# pre-2020 years with the 2020 footprint masks: HF_MASK_YEAR=2020 bash 07_train_and_backfill.sh 1995 2000

if [ -z "${SLURM_JOB_ID:-}" ]; then
  set -euo pipefail
  [ $# -ge 1 ] || { echo "usage: bash 07_train_and_backfill.sh YEAR [YEAR ...]"; exit 1; }
  N_SUB=674
  BIG_SUBS="57 61 62 98 107"
  MEM_SMALL=24G
  MEM_BIG=72G
  BCRS="can10 can11 can12 can13 can14 can3 can40 can41 can42 can5 can60 can61 can70 can71 can72 can80 can81 can82 can9"
  BIG_BCRS="can11 can12 can13 can61"
  MEM_BIG11=128G
  TIME_BIG11=08:00:00

  complement() {
    echo "$2" | awk -v n="$1" 'BEGIN { p = 1; out = "" }
      { if ($1 > p) out = out (out == "" ? "" : ",") ($1 - 1 > p ? p "-" ($1 - 1) : p); p = $1 + 1 }
      END { if (p <= n) out = out (out == "" ? "" : ",") (n > p ? p "-" n : p); print out }'
  }

  N_BCR=0
  BIG11=""
  for b in $BCRS; do
    N_BCR=$((N_BCR + 1))
    if [[ " $BIG_BCRS " == *" $b "* ]]; then BIG11="$BIG11 $N_BCR"; fi
  done
  SMALL=$(complement $N_SUB "$(printf '%s\n' $BIG_SUBS | sort -n)")
  SMALL11=$(complement $N_BCR "$(printf '%s\n' $BIG11)")
  BIG=$(echo $BIG_SUBS | tr ' ' ',')
  BIG11=$(echo $BIG11 | tr ' ' ',')

  for y in "$@"; do
    export YEARS=$y
    J07S=$(sbatch --parsable --export=ALL --mem=$MEM_SMALL --array=${SMALL}%30 "$0")
    J07B=$(sbatch --parsable --export=ALL --mem=$MEM_BIG --array=$BIG "$0")
    DEP=afterany:${J07S%%;*}:${J07B%%;*}
    J11S=$(sbatch --parsable --export=ALL --dependency=$DEP --array=$SMALL11 11_premosaic_backfilled_stacks.sh)
    J11B=$(sbatch --parsable --export=ALL --dependency=$DEP --mem=$MEM_BIG11 --time=$TIME_BIG11 --array=$BIG11 11_premosaic_backfilled_stacks.sh)
    echo "BACKFILL year=$y masks=${HF_MASK_YEAR:-$y} 07=$J07S ($MEM_SMALL: $SMALL) + $J07B ($MEM_BIG: $BIG) 11=$J11S ($SMALL11) + $J11B ($MEM_BIG11: $BIG11), after both 07s"
  done
  exit 0
fi

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla 07_train_and_backfill.R ${SLURM_ARRAY_TASK_ID}
