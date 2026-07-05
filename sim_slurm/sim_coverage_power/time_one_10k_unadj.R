#!/usr/bin/env Rscript
## Time ONLY the UNADJUSTED algo2 @ n_rand=10000, on the same sample production
## would draw for combo 1, chunk 1, sim 1. No covariates are used -- we call
## AR_algo2 on the core (Y, D, Z) design directly, mirroring run_rilate()'s
## without_covariates path. Companion to time_one_10k.R (covariate-adjusted).
suppressPackageStartupMessages(library(devtools))
DIR <- "/vast/palmer/home.grace/hc654/rilate/sim_slurm/sim_coverage_power"
pkg_dir <- "/vast/palmer/home.grace/hc654/rilate"
source(file.path(DIR, "config.R"))
suppressMessages(devtools::load_all(pkg_dir, quiet = TRUE))

combo_id <- 1L; chunk_id <- 1L
cell <- SIM_GRID[combo_id, ]
fractions <- c(cell$fa, cell$fc, cell$fn)
taus <- c(0, cell$tau, 0)
N <- SIM_CFG$N; N1 <- SIM_CFG$N1; N0 <- N - N1
N_at <- round(fractions[1]*N); N_c <- round(fractions[2]*N); N_nt <- N - N_at - N_c
set.seed(SIM_CFG$seed_true + combo_id)
true_table <- gen_data(N_at, N_nt, N_c, taus[1], taus[3], taus[2], SIM_CFG$mode)

sim_global <- (chunk_id-1L)*SIM_CFG$nsim_per_chunk + 1L
set.seed(SIM_CFG$seed_obs + combo_id*10000L + sim_global)
observed <- as.data.frame(gen_data_onesim(true_table, N1, N0))
## (Draw covariates too, only to keep the RNG stream identical to production so
## the zsim seed / sample match the covariate-adjusted run exactly.)
set.seed(SIM_CFG$seed_cov + combo_id*10000L + sim_global)
cov_names <- paste0("X", seq_len(SIM_CFG$n_cov))
for (nm in cov_names) observed[[nm]] <- stats::rnorm(nrow(observed))
zsim_seed <- SIM_CFG$seed_zsim + combo_id*100000L + sim_global

## --- Prepare the WITHOUT-covariates (core) design ---------------------------
std_core <- data.frame(Y_observed = observed$Y_observed,
                       D_observed = observed$D_observed,
                       assignment = observed$assignment)
n_rand <- 10000L
set.seed(zsim_seed)
zsim <- sapply(seq_len(n_rand), gen_assignment_CR_index, N1 = N1, N0 = N0)

cat(sprintf("START %s  host=%s  UNADJUSTED algo2  n_rand=%d\n",
            format(Sys.time()), Sys.info()[["nodename"]], n_rand))
t_unadj <- system.time(
  utils::capture.output(
    res <- AR_algo2(std_core, N1 = N1, N0 = N0, zsim = zsim,
                    tol = 1e-8, alpha = SIM_CFG$alpha),
    file = nullfile())
)[["elapsed"]]

cat(sprintf("DONE  %s\n", format(Sys.time())))
cat(sprintf("  UNADJUSTED algo2 @ n_rand=10000 : %.1f s  (%.2f min)\n",
            t_unadj, t_unadj/60))
cat(sprintf("  peak R mem (gc max used): %.1f GB\n",
            sum(gc()[, "max used"] * c(8, 8)) / 1e9))
