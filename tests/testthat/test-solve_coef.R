# ------------------------------------------------------------------------------
# Unit tests for R/solve_coef.R
#
# NOTE: These unit tests were created by Claude.
#
# Exercises the internal (unexported) function solve_coef_01(), which builds the
# 6 coefficients [a, b, c, d, e, f] of the covariate-adjusted AR rational
# function AR(beta) = (a*beta^2 + b*beta + c) / (d*beta^2 + e*beta + f).
#
# The tests lean on two structural invariants that must hold for ANY input:
#   * the numerator is the perfect square (td*beta - ty)^2, so its discriminant
#     is exactly zero:                     coef[2]^2 == 4 * coef[1] * coef[3]
#   * the denominator is positive-semidefinite (a sum of squares over N^2), so
#     by Cauchy-Schwarz:                    coef[5]^2 <= 4 * coef[4] * coef[6]
#
# Run automatically via devtools::test() / R CMD check, or on their own with
#   testthat::test_file("tests/testthat/test-solve_coef.R")
# ------------------------------------------------------------------------------

# --- Fixtures ----------------------------------------------------------------

# No covariates: with only Y/D/assignment present, "Y ~ . - assignment" reduces
# to an intercept-only model, so intercepts are plain group means and every
# coefficient is computable by hand.
anchor_df <- data.frame(
  Y_observed = c(4, 6, 1, 3),
  D_observed = c(1, 1, 0, 1),
  assignment = c(1, 1, 0, 0)
)
# Treated (rows 1-2): meanY = 5, meanD = 1   -> resy = (-1, 1),  resd = (0, 0)
# Control (rows 3-4): meanY = 2, meanD = 0.5 -> resy = (-1, 1),  resd = (-.5, .5)
#   ty = 5 - 2 = 3 ;  td = 1 - 0.5 = 0.5
#   numerator:   c(td^2, -2*ty*td, ty^2)        = c(0.25, -3, 9)
#   denom(N1=N0=2): coef4 = 0.5/4 = 0.125
#                   coef5 = -2*(0 + 1)/4 = -0.5
#                   coef6 = (2 + 2)/4    = 1
anchor_expected <- c(0.25, -3, 9, 0.125, -0.5, 1)

# A richer, non-degenerate fixture with an informative covariate, used for the
# structural-invariant checks (values are not hand-computed).
general_df <- data.frame(
  Y_observed = c(2.0, 3.5, 3.0, 5.5, 6.0,  1.0, 1.5, 2.0, 2.0, 3.0),
  D_observed = c(0,   1,   0,   1,   1,    0,   0,   1,   0,   1),
  assignment = c(1,   1,   1,   1,   1,    0,   0,   0,   0,   0),
  x          = c(1,   2,   3,   4,   5,    1,   2,   3,   4,   5)
)

# --- Basic shape / type ------------------------------------------------------

test_that("returns a plain (unnamed) numeric vector of length 6", {
  res <- solve_coef_01(anchor_df, N1 = 2, N0 = 2)
  expect_type(res, "double")
  expect_length(res, 6)
  expect_null(names(res))
  expect_true(all(is.finite(res)))
})

# --- Exact values, no covariates ---------------------------------------------

test_that("no-covariate case matches the hand-computed coefficients", {
  res <- solve_coef_01(anchor_df, N1 = 2, N0 = 2)
  expect_equal(res, anchor_expected)
})

# --- Structural invariant: numerator is a perfect square ---------------------

test_that("numerator is the perfect square (td*beta - ty)^2", {
  for (df in list(anchor_df, general_df)) {
    n1 <- sum(df$assignment == 1)
    n0 <- sum(df$assignment == 0)
    res <- solve_coef_01(df, N1 = n1, N0 = n0)
    # a >= 0, c >= 0
    expect_gte(res[1], 0)
    expect_gte(res[3], 0)
    # discriminant of the numerator is exactly zero: b^2 = 4ac
    expect_equal(res[2]^2, 4 * res[1] * res[3])
  }
})

test_that("numerator's double root is the Wald ratio ty/td", {
  res <- solve_coef_01(anchor_df, N1 = 2, N0 = 2)
  root <- -res[2] / (2 * res[1])                 # vertex of the parabola
  expect_equal(root, 3 / 0.5)                     # ty/td = 6
  # numerator evaluated at its double root is ~0
  expect_equal(res[1] * root^2 + res[2] * root + res[3], 0)
})

# --- Structural invariant: denominator is positive-semidefinite --------------

test_that("denominator is PSD: nonneg diagonal and Cauchy-Schwarz on e^2", {
  for (df in list(anchor_df, general_df)) {
    n1 <- sum(df$assignment == 1)
    n0 <- sum(df$assignment == 0)
    res <- solve_coef_01(df, N1 = n1, N0 = n0)
    expect_gte(res[4], 0)
    expect_gte(res[6], 0)
    # e^2 <= 4 d f  (tiny slack for floating point)
    expect_lte(res[5]^2, 4 * res[4] * res[6] + 1e-8)
  }
})

# --- N1 / N0 enter only as an exact 1/N^2 scaling ----------------------------

test_that("denominator scales as 1/N1^2 when the control arm is degenerate", {
  # Control outcomes/compliance are constant -> zero control residuals, so the
  # entire denominator comes from the treated arm and scales cleanly with N1.
  scale_df <- data.frame(
    Y_observed = c(4, 6, 2, 2),
    D_observed = c(1, 0, 0, 0),
    assignment = c(1, 1, 0, 0)
  )
  base   <- solve_coef_01(scale_df, N1 = 2, N0 = 2)
  scaled <- solve_coef_01(scale_df, N1 = 4, N0 = 2)   # double N1 -> quarter denom
  expect_equal(scaled[1:3], base[1:3])                # numerator is N-independent
  expect_equal(scaled[4:6], base[4:6] / 4)
})

test_that("coefficients use the passed N1/N0, not the data's row counts", {
  # Same data, different declared arm sizes -> different denominator.
  a <- solve_coef_01(anchor_df, N1 = 2, N0 = 2)
  b <- solve_coef_01(anchor_df, N1 = 3, N0 = 2)
  expect_equal(a[1:3], b[1:3])                         # numerator unaffected
  expect_false(isTRUE(all.equal(a[4:6], b[4:6])))      # denominator changes
})

# --- Covariate adjustment behaviour ------------------------------------------

test_that("a group-centered covariate leaves the numerator unchanged but shrinks the denominator", {
  base_df <- data.frame(
    Y_observed = c(2, 5, 6, 11,  1, 2, 4, 5),
    D_observed = c(0, 1, 1, 1,   0, 0, 1, 1),
    assignment = c(1, 1, 1, 1,   0, 0, 0, 0)
  )
  # Covariate has mean 0 within each arm, so intercepts stay equal to the group
  # means => identical numerator; but it explains variance => smaller residuals.
  cov_df <- cbind(base_df, x = c(-3, -1, 1, 3,  -3, -1, 1, 3))

  no_cov <- solve_coef_01(base_df, N1 = 4, N0 = 4)
  w_cov  <- solve_coef_01(cov_df,  N1 = 4, N0 = 4)

  expect_equal(w_cov[1:3], no_cov[1:3])               # numerator identical
  expect_lte(w_cov[4], no_cov[4] + 1e-12)             # var(resd) does not grow
  expect_lte(w_cov[6], no_cov[6] + 1e-12)             # var(resy) does not grow
})

# --- Order invariance --------------------------------------------------------

test_that("result is invariant to column order", {
  reordered <- general_df[, c("x", "assignment", "D_observed", "Y_observed")]
  expect_equal(
    solve_coef_01(reordered, N1 = 5, N0 = 5),
    solve_coef_01(general_df, N1 = 5, N0 = 5)
  )
})

test_that("result is invariant to row order", {
  shuffled <- general_df[c(10, 1, 6, 3, 8, 2, 9, 4, 7, 5), ]
  expect_equal(
    solve_coef_01(shuffled, N1 = 5, N0 = 5),
    solve_coef_01(general_df, N1 = 5, N0 = 5)
  )
})

# --- Degenerate inputs -------------------------------------------------------

test_that("constant outcomes/compliance within arms give a zero denominator", {
  const_df <- data.frame(
    Y_observed = c(5, 5, 2, 2),
    D_observed = c(1, 1, 0, 0),
    assignment = c(1, 1, 0, 0)
  )
  res <- solve_coef_01(const_df, N1 = 2, N0 = 2)
  expect_equal(res[4:6], c(0, 0, 0))                  # no residual variance
  expect_equal(res[1:3], c(1, -6, 9))                 # ty = 3, td = 1
})

test_that("N1 = 0 yields non-finite denominator coefficients (no guard)", {
  # Documents current behaviour: division by N1^2 = 0 with nonzero treated
  # residuals produces Inf/NaN rather than an error.
  scale_df <- data.frame(
    Y_observed = c(4, 6, 2, 2),
    D_observed = c(1, 0, 0, 0),
    assignment = c(1, 1, 0, 0)
  )
  res <- solve_coef_01(scale_df, N1 = 0, N0 = 2)
  expect_false(all(is.finite(res)))
})

test_that("a missing required column raises an error", {
  no_d <- data.frame(
    Y_observed = c(4, 6, 1, 3),
    assignment = c(1, 1, 0, 0)
  )
  expect_error(solve_coef_01(no_d, N1 = 2, N0 = 2))
})
