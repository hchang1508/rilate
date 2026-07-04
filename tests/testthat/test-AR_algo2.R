# ------------------------------------------------------------------------------
# Integration tests for R/AR_algo2.R
#
# NOTE: These unit tests were created by Claude.
#
# AR_algo2() is the fast ("jumping") implementation of the Anderson-Rubin
# procedure. Like AR_algo1() it returns a named list with `confidence_set` and
# `p_value` (randomization p-value for beta = 0, LATE = 0, finite-sample-valid
# (1 + count)/(1 + n) convention). These tests mirror the AR_algo1 integration
# tests and additionally check that the two algorithms AGREE on the same inputs.
# ------------------------------------------------------------------------------

# Silence pbapply progress bars so test output stays clean.
if (requireNamespace("pbapply", quietly = TRUE)) {
  pbapply::pboptions(type = "none")
}

# --- Fixtures ----------------------------------------------------------------

quiet2 <- function(expr) {
  utils::capture.output(v <- expr)
  v
}

make_ar_data2 <- function(N = 60, N1 = 30, seed = 1, tau_c = 1) {
  set.seed(seed)
  sim <- gen_sim_data(N = N, N1 = N1,
                      fractions = c(0.2, 0.6, 0.2),
                      taus = c(0, tau_c, 0))
  as.data.frame(sim$observed)
}

prep_inputs2 <- function(df, seed = 123, n_rand = 150) {
  quiet2(run_rilate(df, seed = seed, n_rand = n_rand))
}

run_algo2 <- function(p) {
  quiet2(AR_algo2(p$inputs$without_covariates, N1 = p$N1, N0 = p$N0,
                  zsim = p$zsim, tol = 1e-8, alpha = 0.95))
}

# Represent a confidence set as a matrix for tolerant comparison.
cs_to_matrix <- function(cs) {
  if (length(cs) == 0) matrix(nrow = 0, ncol = 2) else do.call(rbind, cs)
}

# --- Return structure --------------------------------------------------------

test_that("AR_algo2 returns a named list with confidence_set and p_value", {
  res <- run_algo2(prep_inputs2(make_ar_data2()))
  expect_type(res, "list")
  expect_named(res, c("confidence_set", "p_value"))
  expect_true(is.list(res$confidence_set))
  expect_true(is.numeric(res$p_value) && length(res$p_value) == 1)
})

# --- p-value correctness -----------------------------------------------------

test_that("p_value matches an independent (1 + count)/(1 + n) recomputation", {
  p <- prep_inputs2(make_ar_data2())
  res <- run_algo2(p)

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
  res <- run_algo2(prep_inputs2(make_ar_data2()))
  expect_gt(res$p_value, 0)
  expect_lte(res$p_value, 1)
})

# --- Confidence set ----------------------------------------------------------

test_that("every confidence interval is ordered lower <= upper", {
  res <- run_algo2(prep_inputs2(make_ar_data2()))
  for (iv in res$confidence_set) {
    expect_length(iv, 2)
    expect_lte(iv[1], iv[2])
  }
})

# --- Determinism -------------------------------------------------------------

test_that("same inputs -> identical p_value and confidence_set", {
  p <- prep_inputs2(make_ar_data2())
  a <- run_algo2(p)
  b <- run_algo2(p)
  expect_identical(a$p_value, b$p_value)
  expect_identical(a$confidence_set, b$confidence_set)
})

# --- Agreement with AR_algo1 -------------------------------------------------

test_that("AR_algo2 agrees with AR_algo1 on the same inputs", {
  p <- prep_inputs2(make_ar_data2(tau_c = 2))
  r1 <- quiet2(AR_algo1(p$inputs$without_covariates, N1 = p$N1, N0 = p$N0,
                        zsim = p$zsim, tol = 1e-8, alpha = 0.95))
  r2 <- run_algo2(p)
  # identical p-value (same helper, same inputs) ...
  expect_equal(r2$p_value, r1$p_value)
  # ... and the same confidence set (algo2 is the fast equivalent of algo1)
  expect_equal(cs_to_matrix(r2$confidence_set),
               cs_to_matrix(r1$confidence_set),
               tolerance = 1e-6)
})
