#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=96G
#SBATCH --time=03:00:00
#SBATCH --job-name=bf_rootcause
#SBATCH --mail-user=mannfred@ualberta.ca

# READ-ONLY backfill ROOT-CAUSE audit (Open Limitation #5, stage 2).
# Input-side audit on the native 07/08 grid: replicates 08A predictor construction
# (no BART) over all 674 subbasins and decomposes where the high-HF NA originates.
# Writes rootcause_{abiotic_na_by_cov,cascade_by_subbasin,rare_response,bcr_rollup}.csv
# to logs/. Single task (internal loop over subbasins), species-independent.

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

# absolute path so the job works regardless of the sbatch submission directory
# (the R script is already pinned to the cluster's absolute ia_dir internally).
Rscript --vanilla /home/mannfred/scratch/impact_assessment/Rscripts/misc/diag_backfill_rootcause.R
