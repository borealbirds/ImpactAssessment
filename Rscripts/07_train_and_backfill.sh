#!/bin/bash
#SBATCH --account=def-bayne
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=64G
#SBATCH --time=06:00:00
#SBATCH --job-name=2020_run_backfill
# ---------------------------------------------------------------------------
# ARRAY TIERING. Subbasins vary enormously in size. The ~97 largest OOM at the
# 64G/6h settings below and are run separately by 07_train_and_backfill_larger.sh
# (750G/12h) with its own baked-in index list. THE TWO SCRIPTS MUST PARTITION
# 1-674: every index appears in exactly one of them, or the overlapping tasks
# race each other writing the same subbasin_{i}_backfill.tif.
#
# The active line below is that complement, MINUS subbasins 1-3, which the
# 2026-09-14 A5 smoke (job 59865956) already completed at full production
# settings. For a fresh run from an empty bart_models/2020, use the plain
# complement instead (the second commented line).
#
# Regenerate the complement after editing larger.sh's array, from this dir:
#   python3 -c "import re;b=set(int(x) for x in re.search(r'--array=([0-9,]+)',open('07_train_and_backfill_larger.sh').read()).group(1).split(','));r=[];[r.append([i,i]) if not r or i>r[-1][1]+1 else r[-1].__setitem__(1,i) for i in range(1,675) if i not in b];print(','.join(f'{a}-{c}' if c>a else str(a) for a,c in r))"
#
# The %N suffix is the concurrency throttle. These are single-threaded 64G
# tasks (BART runs at tc=1), so %40 is a light load on 750G nodes and finishes
# ~4x sooner than the %10 this script used to carry.
# ---------------------------------------------------------------------------
#SBATCH --array=4-18,20-23,25-26,28-29,31-47,49-52,59-60,65-76,78-81,84-90,93-97,99,101,109-121,123-178,180-183,185-189,191-202,204-219,221-227,229-234,236-238,240-259,261-271,274-275,277,279-295,297-313,315-316,320-323,325-338,340-346,348-351,353-359,361-375,378,380-385,387-396,398-404,406-426,428-437,439-460,462-466,468-471,473-485,487,489-511,513-515,517-529,531-537,539,541-542,544-545,548-552,554-557,561-564,566-571,573-585,587,589-595,597-598,600-602,604-607,609-620,622-628,631-636,639-641,643-650,652-661,663-665,667-671,674%40
# plain complement of larger.sh, all 577 (use for a run from an empty bart_models/2020):
# #SBATCH --array=1-18,20-23,25-26,28-29,31-47,49-52,59-60,65-76,78-81,84-90,93-97,99,101,109-121,123-178,180-183,185-189,191-202,204-219,221-227,229-234,236-238,240-259,261-271,274-275,277,279-295,297-313,315-316,320-323,325-338,340-346,348-351,353-359,361-375,378,380-385,387-396,398-404,406-426,428-437,439-460,462-466,468-471,473-485,487,489-511,513-515,517-529,531-537,539,541-542,544-545,548-552,554-557,561-564,566-571,573-585,587,589-595,597-598,600-602,604-607,609-620,622-628,631-636,639-641,643-650,652-661,663-665,667-671,674%40
# every subbasin, no tiering (the large ones will OOM at 64G):
# #SBATCH --array=1-674%10
#SBATCH --mail-user=mannfred@ualberta.ca

module load StdEnv/2023
module load gcc/12.3
module load gdal/3.9.1
module load udunits/2.2.28
module load r/4.4.0

export NODELIST=$(echo $(srun hostname))
Rscript --vanilla 07_train_and_backfill.R ${SLURM_ARRAY_TASK_ID}

# these settings worked for the vast majority of subbasins
# I ran:
# grep -L "writing .* layers to .*/subbasin_[0-9]\+" Y2020_*.log | \
# xargs -r grep -L "^done$"
# to identify subbasins that were OOM killed or timed out and re-ran with 07_train_and_backfill_larger.sh