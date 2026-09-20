#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=128G
#SBATCH --time=08:00:00
#SBATCH --job-name=2020_run_backfill
#SBATCH --array=1-674%30
#SBATCH --mail-user=mannfred@ualberta.ca
# ---------------------------------------------------------------------------
# ONE TIER, ALL 674 SUBBASINS. This replaces the former 07_train_and_backfill.sh
# / 07_train_and_backfill_larger.sh pair (deleted 2026-09-19) and the two
# hand-maintained complementary --array lists that had to partition 1-674
# exactly or race each other writing the same subbasin_{i}_backfill.tif.
#
# WHY THE TIERS ARE GONE. The split existed because ~97 subbasins OOM'd at 64G.
# That looked like a size effect and was not: 08A assembled the output stack
# with `result_raster[[j]] <- v`, and terra copies the ENTIRE stack on every
# such assignment, so filling ~2000 layers churned ~nlyr^2 * ncell * 8 bytes
# through the allocator. More RAM only bought more garbage before a GC, which
# is why 57/62/98 died on a full 750G node exactly as 454/474/583 died at 64G.
# 08A now fills one ncell x nlyr matrix and calls terra::values() once
# (commit 73e8b68). The six reruns (job 60134338, all COMPLETED 0:0) came back:
#
#   S57   50.8 GB  1:31:24   (was 750 GB, OOM)
#   S62   44.4 GB  1:28:10   (was 750 GB, OOM)
#   S98   50.2 GB  0:47:24   (was 750 GB, OOM)
#   S454   7.4 GB  0:57:21   (was  64 GB, OOM)
#   S474   9.5 GB  2:05:33   (was  64 GB, OOM)
#   S583   7.3 GB  0:28:52   (was  64 GB, OOM)
#
# MEMORY. 57/62/98 were the worst cases in the entire set, so 128G is ~2.5x the
# observed post-fix peak. Peak is roughly (BART working set) + 2 * ncell * nlyr
# * 8 bytes, the factor of 2 being the transient copy inside terra::values()<-;
# for the largest subbasin (S98, ncell=746790, ~2000 layers) that middle term is
# ~24 GB. 96G would very likely suffice -- 128G buys headroom for a subbasin
# bigger than any measured, and the whole run still costs LESS than the old two
# tiers (674 x 128G = 86 TB vs 577 x 64G + 97 x 750G = 110 TB).
#
# If one index ever does OOM, give that index more room rather than re-tiering:
#   sbatch --array=<i> --mem=256G 07_train_and_backfill.sh
#
# TIME. Longest observed post-fix run is 2:05:33 (S474), so 08:00:00 is ~4x
# margin and backfills onto idle nodes sooner than the old 12h large tier did.
#
# CPU. Genuinely single-threaded: 08B calls BART::gbart()/BART::mbart(), not the
# mc.* variants, and nothing sets mc.cores, terraOptions(threads=) or
# OMP_NUM_THREADS. More cores would multiply the allocation charge for nothing.
#
# %30 is the concurrency throttle, not a reservation -- it caps how many array
# tasks run at once, and can be raised or dropped freely.
#
# COMPLETION IS PROVED BY ARTIFACTS, NOT BY SLURM STATE OR BY THE LOGS. 07 wraps
# the call in tryCatch, so an R-level failure still exits COMPLETED, and
# Y2020_S*.log appends across runs -- four of the six subbasins above carried
# "writing N layers" lines from a March run while having no output at all. The
# oracle is subbasin_{i}_confusion.rds, which 08A writes last:
#   find ../data/derived_data/bart_models/2020 -name '*_confusion.rds' | wc -l
# ---------------------------------------------------------------------------

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla 07_train_and_backfill.R ${SLURM_ARRAY_TASK_ID}
