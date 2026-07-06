#!/bin/bash
#SBATCH --job-name=cov_power_rerun
#SBATCH --partition=day
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --mem=32G
#SBATCH --array=91,100,309,310
#SBATCH --output=/vast/palmer/scratch/narita/hc654/cov_power_out/logs/rerun_%A_%a.out
#SBATCH --error=/vast/palmer/scratch/narita/hc654/cov_power_out/logs/rerun_%A_%a.out

# ---------------------------------------------------------------------------
# Resubmit ONLY the 4 tasks that FAILED (exit 7) in job 59053138.
# These are all in combo01 (task_id -> combo = (id-1)/1000 + 1):
#   91 -> combo01 chunk091 ; 100 -> combo01 chunk100
#   309 -> combo01 chunk309 ; 310 -> combo01 chunk310
# Everything else (grid mapping, module, dirs, worker) is identical to
# submit_array.sbatch -- run_worker.R derives combo/chunk from the task id, so
# the same task ids reproduce the same cells. Output writes into the same
# RESULTS_DIR, filling combo01 back to 1000/1000.
# ---------------------------------------------------------------------------

module load R/4.4.1-foss-2022b

# Toggle the n_rand=10000 case (default off -> config.R run_hi=FALSE).
# export RUN_HI=1

DIR=/vast/palmer/home.grace/hc654/rilate/sim_slurm/sim_coverage_power
OUTDIR=/vast/palmer/scratch/narita/hc654/cov_power_out
export RESULTS_DIR="$OUTDIR/results"
cd "$DIR"
mkdir -p "$OUTDIR/logs" "$RESULTS_DIR"

echo "task=${SLURM_ARRAY_TASK_ID} host=$(hostname) cores=${SLURM_CPUS_PER_TASK} out=$RESULTS_DIR start=$(date)"
echo "hw: $(grep -m1 'model name' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//') | $(grep -m1 'cpu MHz' /proc/cpuinfo | cut -d: -f2- | sed 's/^ *//') MHz | node=$(hostname)"

Rscript "$DIR/run_worker.R"
rc=$?
echo "end=$(date) rc=$rc"
exit $rc
