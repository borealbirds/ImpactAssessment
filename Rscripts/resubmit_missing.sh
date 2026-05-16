#!/bin/bash
# Submit all missing/incomplete coalition jobs.
# Skips CIDs where both species are done or already in a named coal_{cid} queue slot.
# Memory: 500G if (CID-1) has >=4 bits set, else 192G.

DENSITY_DIR="$HOME/scratch/impact_assessment/data/derived_data/density_tables"
SCRIPT="$HOME/scratch/impact_assessment/Rscripts/12B_repredict_birds.sh"

completed_cids=$(ls "$DENSITY_DIR" | grep '_coalition_' | \
  sed 's/.*_coalition_\([0-9]*\)_bf_only\.rds$/\1/' | sort -n)

done_cids=$(echo "$completed_cids" | uniq -c | awk '$1 == 2 {print $2}')

queued_cids=$(squeue --me -h -o "%j" | grep -E '^coal_[0-9]+$' | \
  sed 's/coal_//' | sort -n | uniq)

submit_cid() {
  local cid=$1
  local n=$((cid - 1))
  local bits=0
  while [ $n -gt 0 ]; do bits=$((bits + n % 2)); n=$((n / 2)); done
  [ $bits -ge 4 ] && mem=500G || mem=192G

  cat "$SCRIPT" | sbatch \
    --account=def-bayne --mem=$mem \
    --export=ALL,COALITION_ID=$cid \
    --job-name=coal_${cid}
  echo "Submitted CID $cid (${bits} sectors, ${mem})"
}

echo "=== Submitting missing/incomplete coalitions ==="
for cid in $(seq 2 256); do
  if echo "$done_cids"   | grep -qx "$cid"; then continue; fi
  if echo "$queued_cids" | grep -qx "$cid"; then continue; fi
  submit_cid "$cid"
done
echo "=== Done ==="
