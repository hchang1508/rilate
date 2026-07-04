################################################################################
## assess_coverage.R
## Monte-Carlo assessment of test size and confidence-set coverage.
##
## PURPOSE:
## Fix one "true" potential-outcome table, then repeatedly draw a completely
## randomized experiment from it (gen_data_onesim), run the Anderson-Rubin
## procedure (run_rilate) on each draw, and record:
##   * the randomization p-value for H0: LATE = 0, and
##   * whether the true LATE (the complier effect) lies in the confidence set.
## Aggregating over the draws gives the empirical size of the test and the
## empirical coverage of the confidence set.
################################################################################

#' Is a point contained in a confidence set?
#'
#' @param point Scalar value to test for membership.
#' @param confidence_set A list of `c(lower, upper)` intervals (as returned in
#'   `run_rilate()$results$*$confidence_set`), or a single length-2 interval.
#' @param tol Numerical tolerance applied to each interval's endpoints.
#' @return `TRUE` if `point` lies in any interval, otherwise `FALSE`.
#' @noRd
point_in_confidence_set = function(point, confidence_set, tol = 1e-8) {

    if (is.null(confidence_set) || length(confidence_set) == 0) {
        return(FALSE)
    }
    # A degenerate result (e.g. the full real line) may arrive as a bare
    # length-2 numeric vector rather than a list of intervals.
    if (is.numeric(confidence_set) && length(confidence_set) == 2) {
        confidence_set = list(confidence_set)
    }

    for (interval in confidence_set) {
        lo = interval[1]
        hi = interval[2]
        if (point >= lo - tol && point <= hi + tol) {
            return(TRUE)
        }
    }
    FALSE
}

#' Assess test size and confidence-set coverage by simulation
#'
#' Builds a single fixed "true" potential-outcome table for a population split
#' into always-reporters, compliers, and never-reporters, then repeats the
#' following `nsim` times: draw one completely randomized experiment from that
#' table with [gen_data_onesim()], run [run_rilate()] on it, and record the
#' randomization p-value (for the null LATE = 0) and whether the true LATE lies
#' in the confidence set. The true LATE is the complier treatment effect,
#' `taus[2]`.
#'
#' Aggregating across the `nsim` draws yields the empirical **size** of the test
#' (the fraction of draws that reject at level `1 - alpha`) and the empirical
#' **coverage** of the confidence set (the fraction of draws whose set contains
#' the true LATE). When the true LATE is non-zero the reported "size" is really
#' the rejection rate / power of the LATE = 0 test.
#'
#' @param N Total number of units.
#' @param N1 Number of units assigned to treatment (`Z = 1`); the remaining
#'   `N - N1` units are controls.
#' @param fractions Numeric vector or list `c(fa, fc, fn)` of population
#'   fractions for always-reporters, compliers, and never-reporters. Must sum to
#'   1. Counts are `round(fa * N)` and `round(fc * N)`, with the never-reporter
#'   count taken as the remainder so the three counts sum to `N` exactly.
#' @param taus Numeric vector or list `c(tau_a, tau_c, tau_n)` of
#'   principal-stratum treatment effects for always-reporters, compliers, and
#'   never-reporters. The complier effect `tau_c` is the true LATE whose coverage
#'   is assessed.
#' @param mode Either `"constant"` (default) to give every unit in a stratum
#'   exactly its `tau`, or `"heterogeneous"` for unit-level effects calibrated to
#'   each stratum's mean.
#' @param nsim Number of simulated experiments to draw from the fixed true table.
#' @param algorithm Which algorithm [run_rilate()] should use: `"algo2"`
#'   (default, fast) or `"algo1"`.
#' @param n_rand Number of randomizations (permuted assignments) per call to
#'   [run_rilate()]. Default `1000`.
#' @param alpha Confidence level for the confidence set; the test significance
#'   level is `1 - alpha`. Default `0.95`.
#' @param seed Optional integer seed for reproducibility. When supplied, both the
#'   true table and all simulation draws are reproducible.
#' @param verbose If `TRUE`, print a one-line progress update every `10%` of the
#'   run. The (very verbose) internal output of [run_rilate()] is always
#'   suppressed. Default `TRUE`.
#' @return A list with elements: `size` (empirical rejection rate of H0: LATE = 0
#'   at level `1 - alpha`), `coverage` (empirical coverage of the confidence set),
#'   `true_late`, `nsim`, `alpha`, `sig_level`, `counts` (the `N_at`/`N_c`/`N_nt`
#'   stratum sizes), and the per-simulation vectors `p_values`, `covered`, and
#'   `rejected`.
#' @export
assess_coverage = function(N, N1, fractions, taus,
                           mode = "constant",
                           nsim = 100,
                           algorithm = c("algo2", "algo1"),
                           n_rand = 1000,
                           alpha = 0.95,
                           seed = NULL,
                           verbose = TRUE) {

    # --- Unpack and validate inputs ------------------------------------------
    algorithm = match.arg(algorithm)
    mode      = match.arg(mode, c("constant", "heterogeneous"))
    fractions = unlist(fractions)
    taus      = unlist(taus)

    if (length(fractions) != 3) {
        stop("`fractions` must have three elements: (fa, fc, fn).")
    }
    if (length(taus) != 3) {
        stop("`taus` must have three elements: (tau_a, tau_c, tau_n).")
    }
    if (N1 < 0 || N1 > N) {
        stop("`N1` must satisfy 0 <= N1 <= N.")
    }
    if (abs(sum(fractions) - 1) > 1e-8) {
        stop("`fractions` (fa, fc, fn) must sum to 1; got sum = ",
             sum(fractions), ".")
    }
    if (nsim < 1) {
        stop("`nsim` must be a positive integer.")
    }

    # --- Convert fractions to integer counts that sum to N exactly -----------
    N_at = round(fractions[1] * N)
    N_c  = round(fractions[2] * N)
    N_nt = N - N_at - N_c
    if (N_nt < 0) {
        stop("Rounding produced a negative never-reporter count; ",
             "check `fractions` and `N`.")
    }

    if (!is.null(seed)) {
        set.seed(seed)
    }

    # --- Build ONE fixed true table ------------------------------------------
    # gen_data() expects the order (N_at, N_nt, N_c, tau_at, tau_nt, tau_c).
    true_table = gen_data(N_at, N_nt, N_c, taus[1], taus[3], taus[2], mode)

    # The estimand is the complier average treatment effect (LATE).
    true_late = taus[2]
    sig_level = 1 - alpha

    # Silence pbapply progress bars from the algorithm during the loop.
    if (requireNamespace("pbapply", quietly = TRUE)) {
        old_pbo = pbapply::pboptions(type = "none")
        on.exit(pbapply::pboptions(old_pbo), add = TRUE)
    }

    # --- Simulation loop ------------------------------------------------------
    p_values = numeric(nsim)
    covered  = logical(nsim)

    report_every = max(1L, floor(nsim / 10))

    for (s in seq_len(nsim)) {

        # Draw one observed sample from the fixed true table.
        observed = as.data.frame(gen_data_onesim(true_table, N1, N - N1))

        # Run the AR procedure with run_rilate()'s own console output muted.
        utils::capture.output(
            fit <- run_rilate(observed,
                              algorithm      = algorithm,
                              n_rand         = n_rand,
                              with_covariate = FALSE,
                              alpha          = alpha),
            file = nullfile()
        )

        res = fit$results$without_covariates
        p_values[s] = res$p_value
        covered[s]  = point_in_confidence_set(true_late, res$confidence_set)

        if (verbose && (s %% report_every == 0 || s == nsim)) {
            cat(sprintf("  assess_coverage: %d / %d draws done\n", s, nsim))
        }
    }

    rejected = p_values <= sig_level

    # --- Aggregate ------------------------------------------------------------
    size     = mean(rejected)
    coverage = mean(covered)

    if (verbose) {
        cat("========================================\n")
        cat("assess_coverage(): summary\n")
        cat("========================================\n")
        cat("  Simulations (nsim)   :", nsim, "\n")
        cat("  True LATE (tau_c)     :", true_late, "\n")
        cat("  Stratum counts       : at =", N_at, " c =", N_c,
            " nt =", N_nt, "\n")
        cat("  Confidence level     :", alpha,
            "(sig level =", sig_level, ")\n")
        cat("  Size / rejection rate:", round(size, 4),
            "(P[reject H0: LATE = 0])\n")
        cat("  CS coverage          :", round(coverage, 4),
            "(P[true LATE in CS])\n")
        cat("\n")
    }

    list(
        size      = size,
        coverage  = coverage,
        true_late = true_late,
        nsim      = nsim,
        alpha     = alpha,
        sig_level = sig_level,
        counts    = c(N_at = N_at, N_c = N_c, N_nt = N_nt),
        p_values  = p_values,
        covered   = covered,
        rejected  = rejected
    )
}
