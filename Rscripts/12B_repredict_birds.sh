#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=12:00:00
#SBATCH --job-name=coalition_repredict
#SBATCH --array=1-2
#SBATCH --mail-user=mannfred@ualberta.ca

# Two-phase submission:
#
#   Phase 1: Canonical observed predictions (run once per species)
#     sbatch 12A_observed.sh
#
#   Phase 2: Coalition counterfactual predictions (one job per coalition)
#     Run as:  bash 12_repredict_birds.sh
#     This submits one SLURM array job per coalition ID (2-256),
#     with a dependency on the Phase 1 job.
#
# To submit Phase 2 with dependency on a Phase 1 job:
#   OBS_JOB_ID=$(sbatch --parsable 12A_observed.sh)
#   bash 12_repredict_birds.sh $OBS_JOB_ID
#
# To submit just a subset of coalitions (e.g. single-sector coalitions only):
#   for CID in 2 3 5 9 17 33 65 129; do
#     sbatch --export=ALL,COALITION_ID=$CID 12B_repredict_birds.sh
#   done

# if run as `bash 12_repredict_birds.sh [OBS_JOB_ID]` on login node,
# submit one job per coalition and exit
if [ -z "$SLURM_JOB_ID" ]; then
  OBS_JOB_ID=${1:-""}
  DEP_FLAG=""
  if [ -n "$OBS_JOB_ID" ]; then
    DEP_FLAG="--dependency=afterok:$OBS_JOB_ID"
  fi

  # 8 sectors → 2^8 = 256 coalitions; skip ID 1 (empty coalition, v=0)
  for CID in $(seq 2 256); do
    sbatch $DEP_FLAG --export=ALL,COALITION_ID=$CID --job-name=coal_${CID} "$0"
  done
  exit 0
fi

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla 12B_repredict_birds.R
