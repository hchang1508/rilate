# ------------------------------------------------------------------------------
# Integration tests for R/AR_algo1.R
#
# NOTE: These unit tests were created by Claude.
#
# AR_algo1() runs the full Anderson-Rubin permutation procedure and returns a
# named list: `confidence_set` (a list of [lower, upper] intervals) and
# `p_value` (the randomization-based p-value for the null beta = 0, i.e.
# LATE = 0, using the finite-sample-valid (1 + count)/(1 + n) convention).
#
# The deterministic p-value/interval MATH is unit-tested against hand-built
# coefficients in test-AR_helpers.R (ar_value_at, ar_pvalue_at, merge_*, ...).
# These tests exercise the assembled algorithm end-to-end.
# ------------------------------------------------------------------------------

# Silence pbapply progress bars so test output stays clean.
if (requireNamespace("pbapply", quietly = TRUE)) {
  pbapply::pboptions(type = "none")
}

# --- Fixtures ----------------------------------------------------------------

# Run an expression, swallowing its (verbose) console output, and return its
# value.
quiet <- function(expr) {
  utils::capture.output(v <- expr)
  v
}

# A small experiment with noncompliance (20% always-takers, 60% compliers,
# 20% never-takers) and a complier effect of `tau_c`.
make_ar_data <- function(N = 60, N1 = 30, seed = 1, tau_c = 1) {
  set.seed(seed)
  sim <- gen_sim_data(N = N, N1 = N1,
                      fractions = c(0.2, 0.6, 0.2),
                      taus = c(0, tau_c, 0))
  as.data.frame(sim$observed)
}

# Prepare a reproducible (data_without, N1, N0, zsim) bundle via run_rilate's gate.
prep_inputs <- function(df, seed = 123, n_rand = 150) {
  quiet(run_rilate(df, seed = seed, n_rand = n_rand))
}

# --- Return structure --------------------------------------------------------

test_that("AR_algo1 returns a named list with confidence_set and p_value", {
  p <- prep_inputs(make_ar_data())
  res <- quiet(AR_algo1(p$inputs$without_covariates, N1 = p$N1, N0 = p$N0,
                        zsim = p$zsim, tol = 1e-8, alpha = 0.95))
  expect_type(res, "list")
  expect_named(res, c("confidence_set", "p_value"))
  expect_true(is.list(res$confidence_set))
  expect_true(is.numeric(res$p_value) && length(res$p_value) == 1)
})

# --- p-value correctness -----------------------------------------------------

test_that("p_value matches an independent (1 + count)/(1 + n) recomputation", {
  p <- prep_inputs(make_ar_data())
  res <- quiet(AR_algo1(p$inputs$without_covariates, N1 = p$N1, N0 = p$N0,
                        zsim = p$zsim, tol = 1e-8, alpha = 0.95))

  # Recompute AR(0) for observed data and every permutation, independently.
  dt <- p$inputs$without_covariates
  obs <- solve_coef_01(dt, p$N1, p$N0)
  sim <- sapply(seq_len(ncol(p$zsim)), function(i) {
    d <- dt
    d$assignment <- p$zsim[, i]
    solve_coef_01(d, p$N1, p$N0)
  })
  obs0 <- obs[3] / obs[6]
  sim0 <- sim[3, ] / sim[6, ]
  expected <- (1 + sum(sim0 >= obs0)) / (1 + ncol(p$zsim))

  expect_equal(res$p_value, expected)
})

test_that("p_value lies in (0, 1]", {
  p <- prep_inputs(make_ar_data())
  res <- quiet(AR_algo1(p$inputs$without_covariates, N1 = p$N1, N0 = p$N0,
                        zsim = p$zsim, tol = 1e-8, alpha = 0.95))
  expect_gt(res$p_value, 0)
  expect_lte(res$p_value, 1)
})

# --- CS / p-value duality ----------------------------------------------------

test_that("duality: 0 is in the CS  <=>  p_value > alpha_sig", {
  # n_rand chosen so alpha_sig*(n+1) is non-integer => no boundary ties.
  p <- prep_inputs(make_ar_data(tau_c = 2), n_rand = 150)
  res <- quiet(AR_algo1(p$inputs$without_covariates, N1 = p$N1, N0 = p$N0,
                        zsim = p$zsim, tol = 1e-8, alpha = 0.95))
  alpha_sig <- 1 - 0.95
  zero_in_cs <- any(vapply(res$confidence_set,
                           function(iv) iv[1] <= 0 && 0 <= iv[2],
                           logical(1)))
  expect_equal(zero_in_cs, res$p_value > alpha_sig)
})

# --- Confidence set ----------------------------------------------------------

test_that("every confidence interval is ordered lower <= upper", {
  p <- prep_inputs(make_ar_data())
  res <- quiet(AR_algo1(p$inputs$without_covariates, N1 = p$N1, N0 = p$N0,
                        zsim = p$zsim, tol = 1e-8, alpha = 0.95))
  for (iv in res$confidence_set) {
    expect_length(iv, 2)
    expect_lte(iv[1], iv[2])
  }
})

# --- Determinism -------------------------------------------------------------

test_that("same inputs -> identical p_value and confidence_set", {
  p <- prep_inputs(make_ar_data())
  a <- quiet(AR_algo1(p$inputs$without_covariates, N1 = p$N1, N0 = p$N0,
                      zsim = p$zsim, tol = 1e-8, alpha = 0.95))
  b <- quiet(AR_algo1(p$inputs$without_covariates, N1 = p$N1, N0 = p$N0,
                      zsim = p$zsim, tol = 1e-8, alpha = 0.95))
  expect_identical(a$p_value, b$p_value)
  expect_identical(a$confidence_set, b$confidence_set)
})

# --- Signal vs. null (statistical sanity) ------------------------------------

test_that("a strong complier effect yields a smaller p-value than the null", {
  # Same zsim/seed for both so the comparison is apples-to-apples.
  p_null   <- prep_inputs(make_ar_data(tau_c = 0,  seed = 7), seed = 7)
  p_strong <- prep_inputs(make_ar_data(tau_c = 8,  seed = 7), seed = 7)

  r_null <- quiet(AR_algo1(p_null$inputs$without_covariates, N1 = p_null$N1,
                           N0 = p_null$N0, zsim = p_null$zsim,
                           tol = 1e-8, alpha = 0.95))
  r_strong <- quiet(AR_algo1(p_strong$inputs$without_covariates, N1 = p_strong$N1,
                             N0 = p_strong$N0, zsim = p_strong$zsim,
                             tol = 1e-8, alpha = 0.95))
  expect_lt(r_strong$p_value, r_null$p_value)
})
