# TODO

_Last updated: 2026-07-05_

- [x] c. time the two algorithms
      - Measured on a Xeon Gold 6240 core, per sim, each run both unadjusted + covariate-adjusted:
        - **n_rand=1000:** algo1 ≈ **20.6 min/sim**, algo2 ≈ **26 s/sim** → algo1 is ~**47× slower** and dominates the cost.
        - **algo2 @ n_rand=10000:** unadjusted **91.8 min**, covariate-adjusted **110.1 min**; ~7–8 GB peak, single-threaded / CPU-bound (the O(n_rand^2) intersection grid).
      - Cluster is ~6–11× slower per run than an i7-14700F laptop: single-thread (old 2019 Xeon @ ~3 GHz, shared socket) + memory bandwidth. Cluster's value is parallelism, not per-run speed.

- [~] d. assess coverage under the simulated dataset
      - **Running:** SLURM job **59053138** (submitted 2026-07-05), 9000 tasks.
        - 10,000 sims/cell × 9 cells (gamma × tau grid); n_rand=1000 for BOTH algo1 & algo2, each unadjusted + covariate-adjusted; 10 sims/task; 12 h / 32 G.
        - The n_rand=10000 case is deferred (gated behind `RUN_HI=1`; `p3`/`rt3` columns written as NA).
        - Results → `/vast/palmer/scratch/narita/hc654/cov_power_out/results/task_comboNN_chunkNNN.csv` (each row stamps `host` + `cpu_model`, since tasks land on heterogeneous nodes, e.g. 6240 vs 6342).
      - **Next:** when complete, aggregate the per-task CSVs (coverage = mean(coveredK), power = mean(rejectK), for K=1 algo1 / K=2 algo2, unadjusted & `c` = adjusted) and move the summary back to home before scratch is purged (60-day TTL).

- [ ] f. always double-check for stale numbers: whenever the DGP, algorithm, or example changes, re-knit README.Rmd and confirm any hard-coded numbers in the prose (confidence sets, widths, p-values) still match the code output. They drift silently -- ideally pin them by passing `seed=` to `run_rilate()` in the example chunks.
