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

if [ -z "${SLURM_JOB_ID:-}" ]; then
  set -euo pipefail
  [ $# -ge 1 ] || { echo "usage: bash 07_train_and_backfill.sh YEAR [YEAR ...]"; exit 1; }
  N_SUB=674
  N_BCR=19
  BIG_SUBS="57 61 62 98 107"
  MEM_SMALL=24G
  MEM_BIG=72G
  export YEARS=$(IFS=,; echo "$*")

  in_blocks() { for b in $(seq 0 $(($1 - 1))); do for s in $3; do echo $((b * $2 + s)); done; done | sort -n; }
  complement() {
    echo "$2" | awk -v n="$1" 'BEGIN { p = 1; out = "" }
      { if ($1 > p) out = out (out == "" ? "" : ",") ($1 - 1 > p ? p "-" ($1 - 1) : p); p = $1 + 1 }
      END { if (p <= n) out = out (out == "" ? "" : ",") (n > p ? p "-" n : p); print out }'
  }

  BIG=$(in_blocks $# $N_SUB "$BIG_SUBS")
  SMALL=$(complement $((N_SUB * $#)) "$BIG")
  BIG=$(echo $BIG | tr ' ' ',')

  J07S=$(sbatch --parsable --export=ALL --mem=$MEM_SMALL --array=${SMALL}%30 "$0")
  J07B=$(sbatch --parsable --export=ALL --mem=$MEM_BIG --array=$BIG "$0")
  DEP=afterany:${J07S%%;*}:${J07B%%;*}
  J11=$(sbatch --parsable --export=ALL --dependency=$DEP --array=1-$((N_BCR * $#)) 11_premosaic_backfilled_stacks.sh)
  echo "BACKFILL years=$YEARS 07=$J07S ($MEM_SMALL: $SMALL) + $J07B ($MEM_BIG: $BIG) 11=$J11 (1-$((N_BCR * $#))), after both 07s"
  exit 0
fi

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla 07_train_and_backfill.R ${SLURM_ARRAY_TASK_ID}
