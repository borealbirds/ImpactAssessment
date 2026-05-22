#!/bin/bash
# 12B_submit_tiered.sh
# Submit 12B coalition jobs (cid 2..256) for the post-12C-fix re-run, with
# per-coalition --mem/--time from a 3-tier scheme keyed on sector membership.
#
# Resource cost tracks sector FOOTPRINT, not coalition size (empirical, from
# collect_seff.sh -> analyze_seff.R): `roads` (bit 7) dominates; built/crop
# large; dam/mines/oil_gas tiny. Coalition size barely predicts cost.
#
# Tiers (n = cid - 1; bit 7 = roads, bits 0/1 = built/crop):
#   HEAVY  contains roads            (n & 128)        -> 512G / 24:00:00
#   MID    no roads; built/crop or >=4 sectors        -> 320G / 16:00:00
#   LIGHT  otherwise (<=3 small/mid sectors)          -> 224G / 08:00:00
#
# Skip-guard: an array task is skipped if that species' density table already
# exists. Stale buggy-12C tables MUST be deleted before the first run so the
# guard only sees this re-run's output -- the script is then safe to re-run to
# mop up failed / timed-out coalitions (only missing tasks are resubmitted).
#
#   cd /home/mannfred/scratch/impact_assessment/Rscripts
#   DRYRUN=1 bash 12B_submit_tiered.sh    # preview, submit nothing
#   bash 12B_submit_tiered.sh             # submit all coalitions (2..256)
#
# CIDS env var restricts the run to a subset of coalition ids (space-separated).
# For 12F preliminary single-sector results, only the 8 singletons + the
# 6-sector floor are needed -- 9 coalitions instead of 255:
#   CIDS="2 3 5 9 17 33 65 129 64" bash 12B_submit_tiered.sh
set -u
cd "$(dirname "$0")" || exit 1

ACCOUNT="${ACCOUNT:-def-bayne}"
DRYRUN="${DRYRUN:-0}"
CIDS="${CIDS:-$(seq 2 256)}"     # coalition ids to submit (default: full run)
YEAR=2020
SH=12B_repredict_birds.sh
DT_DIR=../data/derived_data/density_tables
SPECIES=(CAWA OVEN)              # array task 1 = CAWA, task 2 = OVEN

[ -f "$SH" ] || { echo "ERROR: $SH not found in $(pwd)"; exit 1; }

nbits() { local n=$1 c=0; while [ "$n" -gt 0 ]; do c=$(( c + n%2 )); n=$(( n/2 )); done; echo "$c"; }

# tier_for <cid> -> echoes "MEM TIME TIER"
tier_for() {
  local n=$(( $1 - 1 )) b
  b=$(nbits "$n")
  if   [ $(( n & 128 )) -ne 0 ];                 then echo "512G 24:00:00 HEAVY"
  elif [ $(( n & 3 )) -ne 0 ] || [ "$b" -ge 4 ]; then echo "320G 16:00:00 MID"
  else                                                echo "224G 08:00:00 LIGHT"
  fi
}

# has_table <species> <cid>  -- final OR bf_only table counts as done
has_table() {
  [ -f "$DT_DIR/$1_${YEAR}_coalition_$2.rds" ]         && return 0
  [ -f "$DT_DIR/$1_${YEAR}_coalition_$2_bf_only.rds" ] && return 0
  return 1
}

echo "=============================================================="
echo " 12B tiered re-run   account=$ACCOUNT  dryrun=$DRYRUN  year=$YEAR"
echo "=============================================================="
n_h=0; n_m=0; n_l=0; n_skip=0
for cid in $CIDS; do
  miss=()
  for t in 1 2; do
    has_table "${SPECIES[$((t-1))]}" "$cid" || miss+=("$t")
  done
  if [ ${#miss[@]} -eq 0 ]; then n_skip=$(( n_skip+1 )); continue; fi
  if [ ${#miss[@]} -eq 2 ]; then arr="1-2"; else arr="${miss[0]}"; fi

  read -r mem tlim tier <<< "$(tier_for "$cid")"
  echo " coal_$cid  $tier  mem=$mem time=$tlim array=$arr"
  cmd="cat $SH | sbatch --account=$ACCOUNT --mem=$mem --time=$tlim --array=$arr --export=ALL,COALITION_ID=$cid --job-name=coal_${cid}"
  if [ "$DRYRUN" = "1" ]; then echo "   DRYRUN: $cmd"; else eval "$cmd"; fi
  case $tier in
    HEAVY) n_h=$(( n_h+1 ));;
    MID)   n_m=$(( n_m+1 ));;
    LIGHT) n_l=$(( n_l+1 ));;
  esac
done
echo "--------------------------------------------------------------"
echo " submitted: HEAVY=$n_h  MID=$n_m  LIGHT=$n_l   skipped(complete)=$n_skip"
echo " track:  squeue -u \$USER -o '%.14i %.12j %.2t %.10M %.10L %R'"
echo "=============================================================="
