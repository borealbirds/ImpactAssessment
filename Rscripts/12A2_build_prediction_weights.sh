#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --job-name=build_pred_weights
#SBATCH --array=1-2
#SBATCH --mail-user=mannfred@ualberta.ca

# Build per-species x BCR prediction weight rasters (range x not-water x in-data-limit).
# Run ONCE before 12B/12D. Reads source masks from data/raw_data/v5_gis (no G: access).
# Adjust --array to match the number of species in species_vec.

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

Rscript --vanilla 12A2_build_prediction_weights.R
