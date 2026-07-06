#!/usr/bin/env Rscript
################################################################################
## summary.R
## Summarise the coverage/power simulation results, one row per parameter combo.
##
## For each of the 9 combos it reports:
##   1. coverage and power (rejection rate), for the unadjusted and the
##      covariate-adjusted analyses, rounded to 3 digits;
##   2. median and max per-sim runtime (seconds) for Algorithm 2 (the default,
##      columns *2) and Algorithm 1 (columns *1).
##
## Usage (from this aggregated/ folder):
##   Rscript summary.R
################################################################################

## Read every combo CSV, in combo order.
files <- sort(list.files(".", pattern = "^combo.*\\.csv$"))

## Summarise one combo file into a single-row data frame.
summarise_combo <- function(file) {
    d <- read.csv(file)
    data.frame(
        combo = d$combo_id[1],
        gamma = d$gamma[1],
        tau   = d$tau[1],
        nsim  = nrow(d),

        ## 1. Coverage and power (unadjusted, then covariate-adjusted).
        cov_unadj  = round(mean(d$covered1),  3),  # coverage of true LATE
        pow_unadj  = round(mean(d$reject1),   3),  # rejection rate of H0: LATE = 0
        cov_adj    = round(mean(d$covered1c), 3),
        pow_adj    = round(mean(d$reject1c),  3),

        ## 2. Runtime in seconds: Algorithm 2 (default) and Algorithm 1.
        algo2_median_s = round(median(d$rt2), 3),
        algo2_max_s    = round(max(d$rt2),    3),
        algo1_median_s = round(median(d$rt1), 3),
        algo1_max_s    = round(max(d$rt1),    3)
    )
}

summary_table <- do.call(rbind, lapply(files, summarise_combo))

print(summary_table, row.names = FALSE)
