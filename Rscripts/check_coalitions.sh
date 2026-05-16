#!/bin/bash
DENSITY_DIR="$HOME/scratch/impact_assessment/data/derived_data/density_tables"

completed=$(ls "$DENSITY_DIR" | grep '_coalition_' | \
  sed 's/.*_coalition_\([0-9]*\)_bf_only\.rds$/\1/' | sort -n)

# coal_{cid} named jobs — CID is in the job name
coal_running=$(squeue --me -h -o "%j" | grep -E '^coal_[0-9]+$' | \
  sed 's/coal_//' | sort -n | uniq)

# coalition_repr jobs — CID not in name; look it up via scontrol
repr_running=""
for jid in $(squeue --me -h -o "%i %j" | awk '/coalition_repr/{print $1}' | \
    sed 's/_.*//' | sort -u); do
  cid=$(scontrol show job "$jid" 2>/dev/null | \
    grep -o 'COALITION_ID=[0-9]*' | sed 's/COALITION_ID=//')
  [ -n "$cid" ] && repr_running="${repr_running}"$'\n'"${cid}"
done
repr_running=$(echo "$repr_running" | grep -v '^$' | sort -n | uniq)

# Suspicious job: coal_ with no CID (likely submission error)
bad_job=$(squeue --me -h -o "%i %j" | awk '$2 == "coal_" {print $1}')
if [ -n "$bad_job" ]; then
  echo "=== WARNING: job named 'coal_' (no CID) found — likely submission error ==="
  scontrol show job "$(echo "$bad_job" | sed 's/_.*//' | head -1)" 2>/dev/null | \
    grep -o 'COALITION_ID=[0-9]*'
  echo ""
fi

running=$(printf '%s\n%s' "$coal_running" "$repr_running" | grep -v '^$' | sort -n | uniq)

echo "=== Species found ==="
ls "$DENSITY_DIR" | grep '_coalition_' | sed 's/_.*//' | sort | uniq

echo ""
echo "=== Files per coalition (N_species  CID) ==="
echo "$completed" | uniq -c

echo ""
echo "=== Coalitions with BOTH species done ==="
echo "$completed" | uniq -c | awk '$1 == 2 {print "  CID " $2}'

echo ""
echo "=== Coalitions with only ONE species done ==="
echo "$completed" | uniq -c | awk '$1 == 1 {print "  CID " $2}'

echo ""
echo "=== coal_{cid} jobs in SLURM queue ==="
echo "$coal_running" | tr '\n' ' '
echo ""

echo ""
echo "=== coalition_repr jobs in SLURM queue (CIDs via scontrol) ==="
if [ -z "$repr_running" ]; then
  n_repr=$(squeue --me -h -o "%j" | grep -c '^coalition_repr' || true)
  echo "  (scontrol did not return COALITION_ID — $n_repr coalition_repr array jobs pending)"
else
  echo "$repr_running" | tr '\n' ' '
fi
echo ""

echo ""
echo "=== Truly MISSING (no files, not in queue) ==="
for cid in $(seq 2 256); do
  in_c=$(echo "$completed" | grep -x "$cid")
  in_r=$(echo "$running"   | grep -x "$cid")
  [ -z "$in_c" ] && [ -z "$in_r" ] && printf "  %d" "$cid"
done
echo ""

echo ""
n_done=$(echo "$completed" | sort -n | uniq | wc -l)
n_coal=$(echo "$coal_running" | grep -c . || true)
n_repr_jobs=$(squeue --me -h -o "%j" | grep -c '^coalition_repr' || true)
echo "Summary: $n_done CIDs have files | $n_coal named coal_* in queue | $n_repr_jobs coalition_repr array jobs in queue"
