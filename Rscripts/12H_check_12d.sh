#!/bin/bash
# ---
# 12H_check_12d.sh  —  READ-ONLY progress probe for a running 12D_combine job
# ---
# Reads each array task's slurm log + sacct elapsed and extrapolates whether
# 12D_combine will finish inside its --time limit. 12D's cost is entirely the
# per-BCR zonal-stats block; the 255-iteration loop skips no-bf_only coalitions
# instantly. So runtime per species ~ (bf_only tables) x (BCRs/coalition) x
# (sec/BCR). This estimates that from progress already logged.
#
# Nothing is cancelled or submitted. Run on a Fir login node:
#   cd /home/mannfred/scratch/impact_assessment/Rscripts
#   bash 12H_check_12d.sh            # defaults to job 40540285
#   bash 12H_check_12d.sh 40567123   # a different 12D job id
# ---
set -u

JOB="${1:-40540285}"
D="$HOME/scratch/impact_assessment"
DT="$D/data/derived_data/density_tables"
SPECIES=(CAWA OVEN)            # task 1 = CAWA, task 2 = OVEN

hms() {  # seconds -> H:MM:SS (or D-HH:MM:SS)
  local s=$1 d h m
  d=$(( s/86400 )); s=$(( s%86400 )); h=$(( s/3600 )); s=$(( s%3600 )); m=$(( s/60 )); s=$(( s%60 ))
  if [ "$d" -gt 0 ]; then printf '%d-%02d:%02d:%02d' "$d" "$h" "$m" "$s"
  else printf '%d:%02d:%02d' "$h" "$m" "$s"; fi
}

echo "=============================================================================="
echo " 12H probe   job=$JOB"
echo "=============================================================================="

for t in 1 2; do
  sp="${SPECIES[$((t-1))]}"

  L="$D/Rscripts/slurm-${JOB}_${t}.out"
  [ -f "$L" ] || L="$D/logs/slurm-${JOB}_${t}.out"

  echo "------------------------------------------------------------------------------"
  echo " task $t  ($sp)   log: $L"
  if [ ! -f "$L" ]; then echo "   (log not found — job not started or path differs)"; continue; fi

  # sacct: elapsed seconds + state + time limit (more reliable than squeue %M)
  read -r eraw state tlim < <(
    sacct -j "${JOB}_${t}" -X -n -P -o ElapsedRaw,State,TimelimitRaw 2>/dev/null \
      | head -n1 | awk -F'|' '{print $1, $2, $3}'
  )
  eraw=${eraw:-0}; state=${state:-UNKNOWN}; tlim=${tlim:-0}
  # TimelimitRaw is in MINUTES; fall back to 4h if unparsed
  limsec=$(( tlim * 60 )); [ "$limsec" -le 0 ] && limsec=14400

  entered=$(grep -c 'bf_only table:'        "$L")
  filled=$( grep -c 'obs columns filled'    "$L")
  done=$(   grep -c 'bf_only file removed'  "$L")
  warn=$(   grep -c 'WARNING: obs still NA' "$L")
  ondisk=$( ls "$DT/${sp}_2020_coalition_"*_bf_only.rds 2>/dev/null | wc -l )

  echo "   state=$state  elapsed=$(hms "$eraw")  limit=$(hms "$limsec")"
  echo "   first: $(head -n1 "$L")"
  echo "   last : $(tail -n1 "$L")"
  echo "   coalitions: entered=$entered  finalized(done=$done + warn=$warn)  bf_only_on_disk=$sp:$ondisk"
  echo "   BCRs filled so far: $filled"

  # ---- extrapolate -----------------------------------------------------------
  # total bf_only this run must process = removed(done) + still-on-disk(ondisk)
  # processed so far                    = done + warn
  # remaining                           = ondisk - warn   (warn ones stay on disk but are done)
  processed=$(( done + warn ))
  total=$(( done + ondisk ))
  remaining=$(( ondisk - warn )); [ "$remaining" -lt 0 ] && remaining=0

  if [ "$processed" -le 0 ] || [ "$eraw" -le 0 ]; then
    echo "   ETA: no coalition finalized yet — still in startup/first coalition;"
    echo "        re-run this probe in ~10 min for a rate."
    continue
  fi

  per=$(awk -v e="$eraw" -v p="$processed" 'BEGIN{printf "%.1f", e/p}')
  eta=$(awk -v r="$remaining" -v e="$eraw" -v p="$processed" 'BEGIN{printf "%d", r*e/p}')
  finish=$(( eraw + eta ))
  margin=$(( limsec - finish ))

  echo "   rate: $(hms "$eraw") / $processed coalitions = ${per}s per coalition"
  echo "   remaining: $remaining of $total coalitions  ->  ETA $(hms "$eta")"
  if [ "$margin" -ge 0 ]; then
    echo "   VERDICT: finishes ~$(hms "$finish") into the job  (+$(hms "$margin") under the limit)  OK"
  else
    over=$(( -margin ))
    echo "   VERDICT: needs ~$(hms "$finish")  -> OVER limit by $(hms "$over")  WILL TIME OUT"
    echo "            (check below whether singleton ids 2,3,5,9,17,33,65,129 + 64"
    echo "             are already finalized — those are all 12F needs.)"
  fi
done

# ---- which 12F inputs are already final (final .rds, no _bf_only) ------------
echo "------------------------------------------------------------------------------"
echo " 12F input readiness  (final table present AND bf_only gone = ready):"
for sp in "${SPECIES[@]}"; do
  line="   $sp:"
  for cid in 2 3 5 9 17 33 65 129 64; do
    f="$DT/${sp}_2020_coalition_${cid}.rds"
    b="$DT/${sp}_2020_coalition_${cid}_bf_only.rds"
    if   [ -f "$f" ] && [ ! -f "$b" ]; then tag="$cid:ready"
    elif [ -f "$f" ] && [ -f "$b" ];   then tag="$cid:final+bf?"
    elif [ -f "$b" ];                  then tag="$cid:bf_only"
    else                                    tag="$cid:MISSING"
    fi
    line="$line $tag"
  done
  echo "$line"
done
echo "=============================================================================="
