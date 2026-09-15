#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --job-name=bf_cov_audit
#SBATCH --mail-user=mannfred@ualberta.ca

# READ-ONLY backfill-mosaic coverage audit (Open Limitation #5 diagnosis).
# Loops internally over both species (CAWA, OVEN) and all BCRs — no array.
# Writes coverage_stage{0..4}*.csv + coverage_localizer_bcr.csv to logs/.

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

# absolute path so the job works regardless of the sbatch submission directory
# (the R script is already pinned to the cluster's absolute ia_dir internally).
Rscript --vanilla /home/mannfred/scratch/impact_assessment/Rscripts/misc/diag_backfill_coverage_audit.R
