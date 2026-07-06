#!/bin/bash
#SBATCH --job-name=cov_power
#SBATCH --partition=day
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --array=1-9000
#SBATCH --output=/vast/palmer/scratch/narita/hc654/cov_power_out/logs/array_%A_%a.out
#SBATCH --error=/vast/palmer/scratch/narita/hc654/cov_power_out/logs/array_%A_%a.out

# ---------------------------------------------------------------------------
# Coverage/power sweep: 9 combos x 1000 chunks = 9000 array tasks.
#   task_id -> combo = (id-1)/1000 + 1 ,  chunk = (id-1)%1000 + 1  (in run_worker.R)
# 10000 sims/cell; each chunk = SIM_CFG$nsim_per_chunk (=10) sims.
#
# THIS FILE: algo1 + algo2 at n_rand=1000, EACH run both unadjusted and
# covariate-adjusted. Measured on a Xeon Gold 6240 core: ALGO1 DOMINATES at
# ~20 min/sim (algo2 is only ~26 s/sim), so ~3.5 h for a 10-sim chunk and
# ~3.3 GB peak. 12 h / 16 G gives margin for slower cells. The p3/rt3 columns
# (the n_rand=10000 case) are still written but NA-filled (stable schema).
#
# To also run the expensive algo2 @ n_rand=10000 case, set RUN_HI=1 below. That
# adds ~2 h/sim (~7 GB peak) on top -- then bump --time toward 24:00:00 and
# --mem to 64G, as the intersection grid (O(n_rand^2)) dominates that case.
#
# NOTE: each sim's row is written incrementally (run_worker.R), so a timeout no
# longer loses the whole chunk -- only the in-flight sim.
# ---------------------------------------------------------------------------

module load R/4.4.1-foss-2022b

# Toggle the n_rand=10000 case (default off -> config.R run_hi=FALSE).
# export RUN_HI=1

DIR=/vast/palmer/home.grace/hc654/rilate/sim_slurm/sim_coverage_power
# Output (results + logs) goes to scratch: ~150 MB / 18k files would eat most of
# the home inode buffer (already at 86% of the 500k file limit). Scratch has
# effectively unlimited headroom. NOTE: scratch is purged after 60 days and not
# backed up -- aggregate the per-task CSVs and move the summary back to home.
OUTDIR=/vast/palmer/scratch/narita/hc654/cov_power_out
export RESULTS_DIR="$OUTDIR/results"
cd "$DIR"
mkdir -p "$OUTDIR/logs" "$RESULTS_DIR"

echo "task=${SLURM_ARRAY_TASK_ID} host=$(hostname) cores=${SLURM_CPUS_PER_TASK} out=$RESULTS_DIR start=$(date)"
# Record the hardware: node CPU model + current clock (per-core speed varies a
# lot across node generations and directly drives the rt1/rt2 timings).
echo "hw: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//') | $(grep -m1 'cpu MHz' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//') MHz | node=$(hostname)"

# Do NOT let the trailing echo mask an Rscript crash: capture the exit status and
# propagate it so a bus-error / OOM task is reported as FAILED, not COMPLETED.
Rscript "$DIR/run_worker.R"
rc=$?
echo "end=$(date) rc=$rc"
exit $rc
