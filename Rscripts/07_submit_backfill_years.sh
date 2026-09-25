#!/bin/bash
# usage: bash 07_submit_backfill_years.sh YEAR [YEAR ...]   (e.g. 2015 2020)
set -euo pipefail
[ $# -ge 1 ] || { echo "usage: bash 07_submit_backfill_years.sh YEAR [YEAR ...]"; exit 1; }
N_SUB=674
N_BCR=19
export YEARS=$(IFS=,; echo "$*")
J07=$(sbatch --parsable --export=ALL --array=1-$((N_SUB * $#))%30 07_train_and_backfill.sh)
J11=$(sbatch --parsable --export=ALL --dependency=afterany:${J07%%;*} --array=1-$((N_BCR * $#)) 11_premosaic_backfilled_stacks.sh)
echo "BACKFILL years=$YEARS 07=$J07 ($((N_SUB * $#)) tasks) 11=$J11 ($((N_BCR * $#)) tasks, after 07)"
