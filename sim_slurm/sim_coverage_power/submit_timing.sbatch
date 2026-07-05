#!/bin/bash
#SBATCH --job-name=cov_timing
#SBATCH --partition=day
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00
#SBATCH --mem=64G
#SBATCH --array=1,1001,2001,3001,4001
#SBATCH --output=/vast/palmer/scratch/narita/hc654/cov_timing_out/logs/timing_%A_%a.out
#SBATCH --error=/vast/palmer/scratch/narita/hc654/cov_timing_out/logs/timing_%A_%a.out

# ---------------------------------------------------------------------------
# EXPLORATORY: measure how long ONE simulation takes, per component, split into
# unadjusted vs covariate-adjusted. 5 array tasks = combos 1..5, chunk 1 (the
# array ids 1,1001,2001,3001,4001 map to combo=(id-1)/1000+1, chunk=1), so we
# sample the timing across 5 different parameter cells. Each task runs a single
# sim of algo1@1000, algo2@1000, algo2@10000, timing each version separately.
# 64G / 6h: the algo2@10000 grid is the memory/time bottleneck we're measuring.
# ---------------------------------------------------------------------------

module load R/4.4.1-foss-2022b

DIR=/vast/palmer/home.grace/hc654/rilate/sim_slurm/sim_coverage_power
OUTDIR=/vast/palmer/scratch/narita/hc654/cov_timing_out
export RESULTS_DIR="$OUTDIR/results"
cd "$DIR"
mkdir -p "$OUTDIR/logs" "$RESULTS_DIR"

echo "task=${SLURM_ARRAY_TASK_ID} host=$(hostname) cores=${SLURM_CPUS_PER_TASK} out=$RESULTS_DIR start=$(date)"
Rscript "$DIR/run_worker_timing.R"
rc=$?
echo "end=$(date) rc=$rc"
exit $rc
