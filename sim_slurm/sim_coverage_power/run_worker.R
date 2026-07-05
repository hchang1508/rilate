#!/usr/bin/env Rscript
################################################################################
## run_worker.R
## One SLURM array task = one (combo, chunk) pair.
##
##   * 9 parameter cells (combo_id 1..9) x 10 chunks (chunk_id 1..10) = 90 tasks
##   * SLURM_ARRAY_TASK_ID in 1..90 maps to:
##       combo_id = (id - 1) %/% 10 + 1
##       chunk_id = (id - 1) %%  10 + 1
##   * each task runs SIM_CFG$nsim_per_chunk (=100) simulations.
##
## For every simulation we draw ONE observed sample from the cell's fixed true
## potential-outcome table, append SIM_CFG$n_cov pure-noise N(0,1) covariates
## (no predictive power), then run these procedures on that SAME sample with the
## SAME permutation seed:
##   1) Algorithm 1, n_rand = SIM_CFG$n_rand    (=1000)
##   2) Algorithm 2, n_rand = SIM_CFG$n_rand    (=1000)
##   3) Algorithm 2, n_rand = SIM_CFG$n_rand_hi (=10000; high-resolution)
##      -- OPTIONAL, gated by SIM_CFG$run_hi (env RUN_HI=1). OFF by default
##         because it costs ~1 h/sim; when off, the p3/rt3 columns are NA-filled
##         so the CSV schema is unchanged.
## Each procedure is run BOTH unadjusted (columns p1/.., suffix none) and
## covariate-adjusted (columns p1c/.., suffix "c"); run_rilate() produces both
## in one call off the same zsim. We record, per case & version, the
## randomization p-value, whether H0: LATE=0 is rejected (power), whether the
## true LATE is covered, the CI width / #intervals / infiniteness; plus the
## per-call wall-clock runtime (covers both versions). We compare the algo1-vs-
## algo2 (n_rand=1000) confidence sets by Hausdorff distance, for both the
## unadjusted (hausdorff) and adjusted (hausdorff_cov) analyses, flagging
## disagreements above SIM_CFG$flag_tol.
##
## Usage:
##   Rscript run_worker.R                 # reads SLURM_ARRAY_TASK_ID
##   Rscript run_worker.R <combo> <chunk> # explicit (local testing)
##   NSIM_PER_CHUNK=3 Rscript run_worker.R 1 1   # small local smoke test
################################################################################

suppressPackageStartupMessages({
    ## Load the rilate package from source (not installed).
    library(devtools)
})

## --- Locate paths ------------------------------------------------------------
this_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA)
if (is.na(this_dir) || is.null(this_dir)) {
    ## Rscript path: derive from --file= argument.
    args_all <- commandArgs(trailingOnly = FALSE)
    fa <- grep("^--file=", args_all, value = TRUE)
    this_dir <- if (length(fa)) dirname(sub("^--file=", "", fa)) else getwd()
}
this_dir <- normalizePath(this_dir)
pkg_dir  <- normalizePath(file.path(this_dir, "..", ".."))   # rilate/
## Output dir: RESULTS_DIR env (set by the sbatch to point at scratch) if
## present, else the local ./results (handy for interactive smoke tests).
res_dir  <- Sys.getenv("RESULTS_DIR", file.path(this_dir, "results"))
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(this_dir, "config.R"))
suppressMessages(devtools::load_all(pkg_dir, quiet = TRUE))

## --- Resolve (combo_id, chunk_id) --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
    combo_id <- as.integer(args[1])
    chunk_id <- as.integer(args[2])
} else {
    task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", ""))
    if (is.na(task_id)) {
        stop("No combo/chunk args and SLURM_ARRAY_TASK_ID is unset.")
    }
    combo_id <- (task_id - 1L) %/% SIM_CFG$chunks + 1L
    chunk_id <- (task_id - 1L) %%  SIM_CFG$chunks + 1L
}

n_combos <- nrow(SIM_GRID)
if (combo_id < 1 || combo_id > n_combos) {
    stop(sprintf("combo_id %d out of range 1..%d", combo_id, n_combos))
}
if (chunk_id < 1 || chunk_id > SIM_CFG$chunks) {
    stop(sprintf("chunk_id %d out of range 1..%d", chunk_id, SIM_CFG$chunks))
}

## Allow small overrides for local smoke tests (production leaves these unset,
## so they default to the config values).
nsim_per_chunk <- as.integer(Sys.getenv("NSIM_PER_CHUNK",
                                        as.character(SIM_CFG$nsim_per_chunk)))
SIM_CFG$n_rand    <- as.integer(Sys.getenv("N_RAND",
                                           as.character(SIM_CFG$n_rand)))
SIM_CFG$n_rand_hi <- as.integer(Sys.getenv("N_RAND_HI",
                                           as.character(SIM_CFG$n_rand_hi)))
## Whether to run the expensive algo2 @ n_rand_hi (=10000) case. Defaults to the
## config value (FALSE); enable for a run with RUN_HI=1 (or true/yes). When off,
## the p3/rt3 columns are still written (as NA) so the CSV schema is stable.
run_hi_env <- Sys.getenv("RUN_HI", "")
if (nzchar(run_hi_env)) {
    SIM_CFG$run_hi <- toupper(run_hi_env) %in% c("1", "TRUE", "T", "YES", "Y")
}

## --- Unpack this cell's parameters -------------------------------------------
cell      <- SIM_GRID[combo_id, ]
fractions <- c(cell$fa, cell$fc, cell$fn)     # (fa, fc, fn)
taus      <- c(0, cell$tau, 0)                # (tau_a, tau_c, tau_n); LATE = tau
true_late <- cell$tau
sig_level <- 1 - SIM_CFG$alpha
N  <- SIM_CFG$N
N1 <- SIM_CFG$N1
N0 <- N - N1

cat(sprintf(
    "[combo %d/%d  chunk %d/%d]  gamma(fc)=%.2f  tau=%.2f  fractions=(%.2f,%.2f,%.2f)  nsim=%d\n",
    combo_id, n_combos, chunk_id, SIM_CFG$chunks,
    cell$fc, cell$tau, fractions[1], fractions[2], fractions[3], nsim_per_chunk))

## --- Build the cell's ONE fixed true table (identical across all chunks) ------
## gen_data() signature: (N_at, N_nt, N_c, tau_at, tau_nt, tau_c, mode).
N_at <- round(fractions[1] * N)
N_c  <- round(fractions[2] * N)
N_nt <- N - N_at - N_c
set.seed(SIM_CFG$seed_true + combo_id)
true_table <- gen_data(N_at, N_nt, N_c, taus[1], taus[3], taus[2], SIM_CFG$mode)

## Run the AR procedure once for a given algorithm, silently, returning BOTH the
## unadjusted and covariate-adjusted analyses. run_rilate() computes both in a
## single call when covariates are supplied and with_covariate = TRUE, sharing
## the same permutation matrix (zsim) across the two.
run_both <- function(observed, algorithm, zsim_seed, cov_names,
                     n_rand = SIM_CFG$n_rand) {
    utils::capture.output(
        fit <- run_rilate(observed,
                          x              = cov_names,
                          algorithm      = algorithm,
                          n_rand         = n_rand,
                          with_covariate = TRUE,
                          alpha          = SIM_CFG$alpha,
                          seed           = zsim_seed,
                          verbose        = FALSE),
        file = nullfile()
    )
    list(u = fit$results$without_covariates,   # unadjusted
         c = fit$results$with_covariates)      # covariate-adjusted
}

## Extract the per-analysis summary metrics from one algorithm result.
## Reads true_late / sig_level from the enclosing script scope.
cs_metrics <- function(res) {
    cs <- res$confidence_set
    list(p    = res$p_value,
         rej  = as.integer(res$p_value <= sig_level),
         cov  = as.integer(point_in_cs(true_late, cs)),
         nint = length(.as_interval_list(cs)),
         w    = cs_width(cs),
         inf  = as.integer(cs_is_infinite(cs)),
         cs   = cs)
}

## --- One simulation (algo1 + algo2 on the same sample) -----------------------
## Factored into a pure function of `s` so the chunk's sims can be run across
## cores. Each sim seeds its own RNG streams deterministically from sim_global,
## so results are independent of the number of cores and the execution order.
run_sim <- function(s) {
    sim_global <- (chunk_id - 1L) * SIM_CFG$nsim_per_chunk + s   # 1..1000

    ## Draw one observed sample (same for all procedures).
    set.seed(SIM_CFG$seed_obs + combo_id * 10000L + sim_global)
    observed <- as.data.frame(gen_data_onesim(true_table, N1, N0))

    ## Append pure-noise N(0,1) covariates (no predictive power): drawn
    ## independently of Y/D/Z from a dedicated seed stream, so covariate
    ## adjustment should be near-neutral (a validity/robustness check).
    set.seed(SIM_CFG$seed_cov + combo_id * 10000L + sim_global)
    cov_names <- paste0("X", seq_len(SIM_CFG$n_cov))
    for (nm in cov_names) observed[[nm]] <- stats::rnorm(nrow(observed))

    ## Same permutation matrix for all procedures (identical zsim seed).
    zsim_seed <- SIM_CFG$seed_zsim + combo_id * 100000L + sim_global

    ## Each call returns BOTH the unadjusted (u) and covariate-adjusted (c)
    ## analyses; rt* times the whole call (both versions together).
    t1 <- system.time(r1 <- run_both(observed, "algo1", zsim_seed, cov_names))["elapsed"]
    t2 <- system.time(r2 <- run_both(observed, "algo2", zsim_seed, cov_names))["elapsed"]

    m1u <- cs_metrics(r1$u); m1c <- cs_metrics(r1$c)
    m2u <- cs_metrics(r2$u); m2c <- cs_metrics(r2$c)

    ## Third case: algo2 with a high-resolution permutation count (n_rand=10000).
    ## OPTIONAL (SIM_CFG$run_hi): it dominates runtime (~1 h/sim), so it is off
    ## by default. When off, emit NA metrics/timing so the CSV schema is stable.
    na_metrics <- list(p = NA_real_, rej = NA_integer_, cov = NA_integer_,
                       nint = NA_integer_, w = NA_real_, inf = NA_integer_,
                       cs = NULL)
    if (isTRUE(SIM_CFG$run_hi)) {
        t3 <- system.time(
            r3 <- run_both(observed, "algo2", zsim_seed, cov_names,
                           n_rand = SIM_CFG$n_rand_hi)
        )["elapsed"]
        m3u <- cs_metrics(r3$u); m3c <- cs_metrics(r3$c)
    } else {
        t3  <- NA_real_
        m3u <- na_metrics; m3c <- na_metrics
    }

    hd   <- hausdorff_cs(m1u$cs, m2u$cs)   # unadjusted: algo1 vs algo2
    hd_c <- hausdorff_cs(m1c$cs, m2c$cs)   # adjusted:   algo1 vs algo2

    data.frame(
        combo_id   = combo_id,
        gamma      = cell$fc,
        fc         = cell$fc,
        tau        = cell$tau,
        chunk_id   = chunk_id,
        sim_global = sim_global,
        true_late  = true_late,
        n_cov      = SIM_CFG$n_cov,
        # -- Algorithm 1 @ n_rand=1000 -----------------------------------------
        p1   = m1u$p, reject1  = m1u$rej, covered1  = m1u$cov,
        nint1  = m1u$nint, width1  = m1u$w, inf1  = m1u$inf,
        p1c  = m1c$p, reject1c = m1c$rej, covered1c = m1c$cov,
        nint1c = m1c$nint, width1c = m1c$w, inf1c = m1c$inf,
        rt1  = as.numeric(t1),
        # -- Algorithm 2 @ n_rand=1000 -----------------------------------------
        p2   = m2u$p, reject2  = m2u$rej, covered2  = m2u$cov,
        nint2  = m2u$nint, width2  = m2u$w, inf2  = m2u$inf,
        p2c  = m2c$p, reject2c = m2c$rej, covered2c = m2c$cov,
        nint2c = m2c$nint, width2c = m2c$w, inf2c = m2c$inf,
        rt2  = as.numeric(t2),
        # -- Algorithm 2 @ n_rand=10000 (high-resolution) ----------------------
        p3   = m3u$p, reject3  = m3u$rej, covered3  = m3u$cov,
        nint3  = m3u$nint, width3  = m3u$w, inf3  = m3u$inf,
        p3c  = m3c$p, reject3c = m3c$rej, covered3c = m3c$cov,
        nint3c = m3c$nint, width3c = m3c$w, inf3c = m3c$inf,
        rt3  = as.numeric(t3),
        # -- Comparisons (algo1 vs algo2) --------------------------------------
        hausdorff     = hd,
        flag_diff     = as.integer(is.infinite(hd)   || hd   > SIM_CFG$flag_tol),
        hausdorff_cov = hd_c,
        flag_diff_cov = as.integer(is.infinite(hd_c) || hd_c > SIM_CFG$flag_tol),
        stringsAsFactors = FALSE
    )
}

## --- Simulation loop (serial: one sim after another on a single core) --------
## Each sim runs 3 procedures x {unadjusted, covariate-adjusted}; the algo2@10k
## pair dominates (~1 h, ~30 min/run). With nsim_per_chunk = 1 that is ~1.3 h/task.
## Write each sim's row to the CSV as soon as it finishes, so a crash (OOM /
## bus error) partway through the chunk keeps the sims already completed instead
## of losing all of them. Accumulate into a .part file and rename to the final
## name only on success, so a consumer never picks up a truncated file.
out_file  <- file.path(res_dir,
                       sprintf("task_combo%02d_chunk%03d.csv", combo_id, chunk_id))
part_file <- paste0(out_file, ".part")
cat(sprintf("Running %d sims serially ...\n", nsim_per_chunk))
rows <- vector("list", nsim_per_chunk)
t_loop <- system.time({
    for (s in seq_len(nsim_per_chunk)) {
        row <- run_sim(s)
        rows[[s]] <- row
        write.table(row, part_file, sep = ",", row.names = FALSE,
                    col.names = (s == 1L), append = (s != 1L))
    }
})["elapsed"]

out <- do.call(rbind, rows)
invisible(file.rename(part_file, out_file))   # atomic promotion to final name
cat(sprintf("All %d sims done in %.1fs wall.\n", nsim_per_chunk, t_loop))

cat(sprintf(
    "DONE combo %d chunk %d -> %s\n  [unadj] cov1=%.3f cov2=%.3f cov3=%.3f  pow1=%.3f pow2=%.3f pow3=%.3f\n  [adj]   cov1=%.3f cov2=%.3f cov3=%.3f  pow1=%.3f pow2=%.3f pow3=%.3f\n  flags(u/c)=%d/%d  rt: algo1=%.1fs algo2=%.1fs algo2@10k=%.1fs  (each rt incl. both versions)\n",
    combo_id, chunk_id, out_file,
    mean(out$covered1), mean(out$covered2), mean(out$covered3),
    mean(out$reject1),  mean(out$reject2),  mean(out$reject3),
    mean(out$covered1c), mean(out$covered2c), mean(out$covered3c),
    mean(out$reject1c),  mean(out$reject2c),  mean(out$reject3c),
    sum(out$flag_diff), sum(out$flag_diff_cov),
    sum(out$rt1), sum(out$rt2), sum(out$rt3)))
