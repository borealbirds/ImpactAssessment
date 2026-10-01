#!/bin/bash
# usage: bash 07_submit_backfill_years.sh YEAR [YEAR ...]   (e.g. 2015 2020)
set -euo pipefail
[ $# -ge 1 ] || { echo "usage: bash 07_submit_backfill_years.sh YEAR [YEAR ...]"; exit 1; }
N_SUB=674
N_BCR=19
BIG_SUBS="57 61 62 98 107"
MEM_SMALL=24G
MEM_BIG=72G
LONG_BCRS="6"
TIME_LONG=24:00:00
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
LONG=$(in_blocks $# $N_BCR "$LONG_BCRS")
SHORT=$(complement $((N_BCR * $#)) "$LONG")
LONG=$(echo $LONG | tr ' ' ',')

J07S=$(sbatch --parsable --export=ALL --mem=$MEM_SMALL --array=${SMALL}%30 07_train_and_backfill.sh)
J07B=$(sbatch --parsable --export=ALL --mem=$MEM_BIG --array=$BIG 07_train_and_backfill.sh)
DEP=afterany:${J07S%%;*}:${J07B%%;*}
J11S=$(sbatch --parsable --export=ALL --dependency=$DEP --array=$SHORT 11_premosaic_backfilled_stacks.sh)
J11L=$(sbatch --parsable --export=ALL --dependency=$DEP --time=$TIME_LONG --array=$LONG 11_premosaic_backfilled_stacks.sh)
echo "BACKFILL years=$YEARS 07=$J07S ($MEM_SMALL: $SMALL) + $J07B ($MEM_BIG: $BIG) 11=$J11S ($SHORT) + $J11L ($TIME_LONG: $LONG), after both 07s"
