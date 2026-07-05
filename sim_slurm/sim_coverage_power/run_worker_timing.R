#!/usr/bin/env Rscript
################################################################################
## run_worker_timing.R  --  EXPLORATORY timing worker.
##
## Runs exactly ONE simulation for a given (combo, chunk) and reports the
## wall-clock cost of each component, split into unadjusted vs covariate-
## adjusted:
##   1) algo1 @ n_rand
##   2) algo2 @ n_rand
##   3) algo2 @ n_rand_hi   (the suspected bottleneck)
##
## Timing comes from run_rilate()'s per-run `runtimes` (which times the
## without_covariates and with_covariates analyses separately in one call), so
## nothing is recomputed. Also records the whole-call elapsed (incl. setup /
## zsim generation) per component.
##
## Uses the SAME sample/seed construction as run_worker.R so these timings are
## representative of a real production sim (n_rand / n_rand_hi default to the
## config's production values -- do NOT shrink them here).
##
## Usage:
##   Rscript run_worker_timing.R                 # reads SLURM_ARRAY_TASK_ID
##   Rscript run_worker_timing.R <combo> <chunk> # explicit (local)
################################################################################

suppressPackageStartupMessages(library(devtools))

## --- Locate paths ------------------------------------------------------------
this_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) NA)
if (is.na(this_dir) || is.null(this_dir)) {
    args_all <- commandArgs(trailingOnly = FALSE)
    fa <- grep("^--file=", args_all, value = TRUE)
    this_dir <- if (length(fa)) dirname(sub("^--file=", "", fa)) else getwd()
}
this_dir <- normalizePath(this_dir)
pkg_dir  <- normalizePath(file.path(this_dir, "..", ".."))
res_dir  <- Sys.getenv("RESULTS_DIR", file.path(this_dir, "timing_results"))
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

source(file.path(this_dir, "config.R"))
suppressMessages(devtools::load_all(pkg_dir, quiet = TRUE))

## --- Resolve (combo_id, chunk_id) --------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) >= 2) {
    combo_id <- as.integer(args[1]); chunk_id <- as.integer(args[2])
} else {
    task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", ""))
    if (is.na(task_id)) stop("No combo/chunk args and SLURM_ARRAY_TASK_ID is unset.")
    combo_id <- (task_id - 1L) %/% SIM_CFG$chunks + 1L
    chunk_id <- (task_id - 1L) %%  SIM_CFG$chunks + 1L
}
n_combos <- nrow(SIM_GRID)
stopifnot(combo_id >= 1, combo_id <= n_combos,
          chunk_id >= 1, chunk_id <= SIM_CFG$chunks)

## Production n_rand values by default (overridable for a quick local check).
SIM_CFG$n_rand    <- as.integer(Sys.getenv("N_RAND",    as.character(SIM_CFG$n_rand)))
SIM_CFG$n_rand_hi <- as.integer(Sys.getenv("N_RAND_HI", as.character(SIM_CFG$n_rand_hi)))

## --- Unpack this cell's parameters -------------------------------------------
cell      <- SIM_GRID[combo_id, ]
fractions <- c(cell$fa, cell$fc, cell$fn)
taus      <- c(0, cell$tau, 0)
N  <- SIM_CFG$N; N1 <- SIM_CFG$N1; N0 <- N - N1

cat(sprintf(
    "[TIMING combo %d/%d chunk %d]  fc=%.2f tau=%.2f  n_rand=%d  n_rand_hi=%d\n",
    combo_id, n_combos, chunk_id, cell$fc, cell$tau,
    SIM_CFG$n_rand, SIM_CFG$n_rand_hi))

## --- Build the cell's fixed true table (same as production) ------------------
N_at <- round(fractions[1] * N); N_c <- round(fractions[2] * N)
N_nt <- N - N_at - N_c
set.seed(SIM_CFG$seed_true + combo_id)
true_table <- gen_data(N_at, N_nt, N_c, taus[1], taus[3], taus[2], SIM_CFG$mode)

## --- Draw ONE observed sample (sim s = 1), same seeds as production ----------
sim_global <- (chunk_id - 1L) * SIM_CFG$nsim_per_chunk + 1L
set.seed(SIM_CFG$seed_obs + combo_id * 10000L + sim_global)
observed <- as.data.frame(gen_data_onesim(true_table, N1, N0))
set.seed(SIM_CFG$seed_cov + combo_id * 10000L + sim_global)
cov_names <- paste0("X", seq_len(SIM_CFG$n_cov))
for (nm in cov_names) observed[[nm]] <- stats::rnorm(nrow(observed))
zsim_seed <- SIM_CFG$seed_zsim + combo_id * 100000L + sim_global

## --- Time one component: run_rilate() times unadj vs adj internally ----------
time_component <- function(label, algorithm, n_rand) {
    t_call <- system.time(
        utils::capture.output(
            fit <- run_rilate(observed, x = cov_names, algorithm = algorithm,
                              n_rand = n_rand, with_covariate = TRUE,
                              alpha = SIM_CFG$alpha, seed = zsim_seed,
                              verbose = FALSE),
            file = nullfile())
    )[["elapsed"]]
    rt <- fit$runtimes                       # named: without_covariates / with_covariates
    t_unadj <- as.numeric(rt[["without_covariates"]])
    t_adj   <- as.numeric(rt[["with_covariates"]])
    cat(sprintf("  %-14s n_rand=%-5d  unadj=%8.1fs  adj=%8.1fs  (adj/unadj=%.2fx)  call=%8.1fs\n",
                label, n_rand, t_unadj, t_adj,
                if (t_unadj > 0) t_adj / t_unadj else NA_real_, t_call))
    data.frame(component = label, n_rand = n_rand,
               t_unadj = t_unadj, t_adj = t_adj, t_call = as.numeric(t_call),
               stringsAsFactors = FALSE)
}

cat("Timing 1 sim across components (each times unadjusted vs covariate-adjusted)...\n")
t_total <- system.time({
    rows <- rbind(
        time_component("algo1",      "algo1", SIM_CFG$n_rand),
        time_component("algo2",      "algo2", SIM_CFG$n_rand),
        time_component("algo2_hi",   "algo2", SIM_CFG$n_rand_hi)
    )
})[["elapsed"]]

rows$combo_id   <- combo_id
rows$chunk_id   <- chunk_id
rows$fc         <- cell$fc
rows$tau        <- cell$tau
rows$sim_global <- sim_global
rows$host       <- Sys.info()[["nodename"]]
rows <- rows[, c("combo_id", "chunk_id", "fc", "tau", "sim_global", "host",
                 "component", "n_rand", "t_unadj", "t_adj", "t_call")]

out_file <- file.path(res_dir, sprintf("timing_combo%02d_chunk%03d.csv",
                                       combo_id, chunk_id))
write.csv(rows, out_file, row.names = FALSE)

cat(sprintf(
    "DONE 1 sim in %.1fs wall -> %s\n  total unadj=%.1fs  total adj=%.1fs  grand=%.1fs\n",
    t_total, out_file,
    sum(rows$t_unadj), sum(rows$t_adj), sum(rows$t_unadj) + sum(rows$t_adj)))
