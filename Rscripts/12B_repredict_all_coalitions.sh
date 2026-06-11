#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=384G
#SBATCH --time=24:00:00
#SBATCH --job-name=coal_all
#SBATCH --array=1-2
#SBATCH --mail-user=mannfred@ualberta.ca

# Restructured Phase 2: ONE job per species (array task 1=CAWA, 2=OVEN) computes
# ALL 255 coalitions internally (DESIGN_12C_restructure.md). Replaces the 255-jobs
# -per-species 12B_repredict_birds.sh + 12B_submit_tiered.sh fan-out.
#
# Resources: --mem/--time are an initial generous envelope (the per-BCR superset
# compute ~ one all-8 coalition + 255 cheap reductions). Re-profile with
# collect_seff.sh -> analyze_seff.R on the new jobs and tighten.
#
#   cd /home/mannfred/scratch/impact_assessment/Rscripts
#   sbatch 12B_repredict_all_coalitions.sh                 # both species
#   sbatch --array=1 12B_repredict_all_coalitions.sh       # CAWA only
#
# Targeted smoke test (one small BCR, 2 bootstraps) before a full run:
#   sbatch --array=1 --time=01:00:00 --mem=192G \
#     --export=ALL,TEST_BCR=can60,TEST_N_BOOT=2 12B_repredict_all_coalitions.sh

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla 12B_repredict_all_coalitions.R
