#!/bin/bash
# usage: bash run_species_batches.sh prep | smoke | plan | wave N   (species: species_todo.txt; MAX=400 tasks per wave)
# years: YEARS=1995,2000 bash run_species_batches.sh ...   (12A/12C/12D/12H read YEARS; default 2020; MAX counts one year)
set -euo pipefail
cd "$(dirname "$0")"
BOOT=${BOOT:-/home/mannfred/projects/def-ecknight/NationalModels/output/06_bootstraps}
LIST=${LIST:-species_todo.txt}
MAX=${MAX:-400}
SPP=$(tr -d '\r' < "$LIST" | grep -v '^#' | awk 'NF')

ncan() { ls "$BOOT/$1" | grep -cE "^$1_can[0-9]+\.Rdata$" || true; }

waves() {
  local w=1 n=0 out="" sp k
  for sp in $SPP; do
    k=$(ncan "$sp")
    if [ "$n" -gt 0 ] && [ $((n + k)) -gt "$MAX" ]; then
      echo "$w $n$out"
      w=$((w + 1)); n=0; out=""
    fi
    n=$((n + k)); out="$out $sp"
  done
  if [ "$n" -gt 0 ]; then echo "$w $n$out"; fi
}

case "${1:-}" in
  prep)
    bash 12A_observed.sh $SPP
    bash 12C_build_prediction_weights.sh $SPP
    ;;
  smoke)
    S=""; c=0
    for sp in $SPP; do
      if [ -f "$BOOT/$sp/${sp}_can42.Rdata" ] && [ $c -lt 3 ]; then S="$S $sp"; c=$((c + 1)); fi
    done
    [ -n "$S" ] || { echo "no species with a can42 model"; exit 1; }
    TEST_BCR=can42 TEST_N_BOOT=2 bash 12D_repredict_all_coalitions.sh $S
    ;;
  plan)
    waves | awk '{printf "wave %s: %s tasks, %d species:", $1, $2, NF - 2; for (i = 3; i <= NF; i++) printf " %s", $i; print ""}'
    echo "total: $(waves | awk '{t += $2} END {print t}') tasks"
    ;;
  wave)
    W=$(waves | awk -v w="${2:?wave number}" '$1 == w {for (i = 3; i <= NF; i++) printf "%s ", $i}')
    [ -n "$W" ] || { echo "no wave $2 (see: bash run_species_batches.sh plan)"; exit 1; }
    bash 12D_repredict_all_coalitions.sh $W
    ;;
  *)
    echo "usage: bash run_species_batches.sh prep | smoke | plan | wave N"; exit 1
    ;;
esac
