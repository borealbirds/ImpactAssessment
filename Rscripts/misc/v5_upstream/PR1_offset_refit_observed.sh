#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=03:00:00
#SBATCH --job-name=pr1_offsets
#SBATCH --mail-user=mannfred@ualberta.ca
# submit (from Rscripts/misc/v5_upstream): bash PR1_offset_refit_observed.sh   (defaults SPECIES=OVEN,SWTH,REVI BCR=can81 N_BOOT=8)

if [ -z "${SLURM_JOB_ID:-}" ]; then
  set -euo pipefail
  export SPECIES=${SPECIES:-OVEN,SWTH,REVI}
  export BCR=${BCR:-can81}
  export N_BOOT=${N_BOOT:-8}
  LOG=/home/mannfred/scratch/impact_assessment/offset_test/logs
  mkdir -p "$LOG"
  N=$(( $(echo "$SPECIES" | tr ',' '\n' | grep -c .) * 2 ))
  R=$(STEP=run sbatch --parsable --export=ALL --array=1-$N --output=$LOG/pr1_run_%A_%a.out "$0")
  S=$(STEP=summary sbatch --parsable --export=ALL --cpus-per-task=1 --mem=8G --time=00:30:00 --dependency=afterany:$R --output=$LOG/pr1_summary_%A.out "$0")
  echo "PR1 species=$SPECIES bcr=$BCR n_boot=$N_BOOT run=$R (1-$N) summary=$S logs=$LOG"
  exit 0
fi

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla PR1_offset_refit_observed.R
