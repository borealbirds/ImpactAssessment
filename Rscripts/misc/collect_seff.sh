#!/bin/bash
# collect_seff.sh
# Collect SLURM resource accounting (sacct) for every 12B coalition job that
# left a slurm-*.out file in this directory. Writes seff_summary.psv
# (pipe-delimited) for local parsing into per-coalition-size mem/time stats.
#
# Run on a Fir login node:
#   bash /home/mannfred/scratch/impact_assessment/Rscripts/collect_seff.sh

cd /home/mannfred/scratch/impact_assessment/Rscripts || exit 1

ids=$(ls slurm-*.out 2>/dev/null | sed -E 's/^slurm-([0-9]+)_.*/\1/' | sort -un | paste -sd,)
if [ -z "$ids" ]; then
  echo "ERROR: no slurm-*.out files found in $(pwd)"
  exit 1
fi

# explicit --format overrides any personal SACCT_FORMAT default;
# JobName carries coal_<id>, which the local parser maps to coalition size.
sacct -j "$ids" --parsable2 --units=M \
  --format=JobID,JobName,State,Elapsed,ElapsedRaw,Timelimit,ReqMem,MaxRSS,MaxVMSize,TotalCPU,AllocCPUS,NNodes,Start,End \
  > seff_summary.psv

echo "wrote $(pwd)/seff_summary.psv  ($(wc -l < seff_summary.psv) rows)"
