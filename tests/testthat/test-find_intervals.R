# ------------------------------------------------------------------------------
# Unit tests for R/find_intervals.R
#
# NOTE: These unit tests were created by Claude.
#
# They exercise the internal (unexported) functions that build acceptance
# intervals from a pair of AR rational functions:
#   - compare_AR()        TRUE iff AR1(beta) < AR2(beta) at a point
#   - find_intervals_03() intervals where AR1(beta) < AR2(beta) (a 2xk matrix,
#                         row 1 = left endpoints, row 2 = right endpoints)
#
# testthat makes these available directly because the package is loaded for the
# test run. Run automatically via devtools::test() or R CMD check, or on their
# own with
#   testthat::test_file("tests/testthat/test-find_intervals.R")
# ------------------------------------------------------------------------------

# --- Small helpers -----------------------------------------------------------
# AR value of coef = c(a, b, c, d, e, f) at a given beta:
#   AR(beta) = (a*beta^2 + b*beta + c) / (d*beta^2 + e*beta + f)
ar_value <- function(coef, beta) {
  (coef[1] * beta^2 + coef[2] * beta + coef[3]) /
    (coef[4] * beta^2 + coef[5] * beta + coef[6])
}

# Coefficient builders (numerator over a denominator of 1):
CONST  <- function(k)    c(0, 0, k, 0, 0, 1)  # AR(beta) = k
LINEAR <- function(m, k) c(0, m, k, 0, 0, 1)  # AR(beta) = m*beta + k
SQUARE <- c(1, 0, 0, 0, 0, 1)                 # AR(beta) = beta^2

# --- compare_AR --------------------------------------------------------------

test_that("compare_AR is TRUE when the first AR is below the second", {
  expect_true(compare_AR(0, CONST(0), CONST(1)))   # 0 < 1
  expect_false(compare_AR(0, CONST(2), CONST(1)))  # 2 < 1 is false
})

test_that("compare_AR uses a strict inequality at a crossing point", {
  # AR1(beta) = beta and AR2(beta) = 1 both equal 1 at beta = 1.
  expect_false(compare_AR(1, LINEAR(1, 0), CONST(1)))
})

test_that("compare_AR agrees with a direct evaluation on the standard branch", {
  beta <- 3
  direct <- ar_value(SQUARE, beta) < ar_value(CONST(1), beta)
  # compare_AR returns a 1x1 matrix (from %*%); coerce before comparing.
  expect_equal(as.vector(compare_AR(beta, SQUARE, CONST(1))), direct)
})

test_that("compare_AR handles very large |beta| via the 1/beta branch", {
  # |beta| > 1e8 triggers the overflow-safe transformation.
  # AR1(beta) = beta, AR2(beta) = 2*beta.
  expect_true(compare_AR(1e9, LINEAR(1, 0), LINEAR(2, 0)))    # beta < 2 beta
  expect_false(compare_AR(-1e9, LINEAR(1, 0), LINEAR(2, 0)))  # -b < -2b false
})

# --- find_intervals_03 -------------------------------------------------------

# Helper: check that a returned matrix is a well-formed 2-row endpoint matrix.
expect_endpoint_matrix <- function(m) {
  expect_true(is.matrix(m))
  expect_equal(nrow(m), 2L)
  expect_true(all(m[1, ] <= m[2, ]))  # left endpoint <= right endpoint
}

test_that("no crossing, first AR below everywhere -> the whole real line", {
  # AR1 = 0 < AR2 = 1 for every beta, and the curves never cross.
  m <- find_intervals_03(CONST(0), CONST(1))
  expect_endpoint_matrix(m)
  expect_equal(unname(m), matrix(c(-Inf, Inf), nrow = 2))
})

test_that("no crossing, first AR above everywhere -> empty acceptance set", {
  # AR1 = 2 is never below AR2 = 1; nothing is accepted.
  expect_null(find_intervals_03(CONST(2), CONST(1)))
})

test_that("a single crossing produces one half-line", {
  # AR1(beta) = beta, AR2(beta) = 1, crossing at beta = 1.
  # beta < 1  ->  beta < 1 (accept);  beta > 1  ->  reject.
  m <- find_intervals_03(LINEAR(1, 0), CONST(1))
  expect_endpoint_matrix(m)
  expect_equal(unname(m), matrix(c(-Inf, 1), nrow = 2))
})

test_that("two crossings with a low middle -> one bounded interval", {
  # beta^2 < 1  iff  -1 < beta < 1.
  m <- find_intervals_03(SQUARE, CONST(1))
  expect_endpoint_matrix(m)
  expect_equal(unname(m), matrix(c(-1, 1), nrow = 2))
})

test_that("two crossings with a high middle -> two disjoint half-lines", {
  # 1 < beta^2  iff  beta < -1 or beta > 1.
  m <- find_intervals_03(CONST(1), SQUARE)
  expect_endpoint_matrix(m)
  expect_equal(ncol(m), 2L)
  expect_equal(unname(m), matrix(c(-Inf, -1, 1, Inf), nrow = 2))
})

test_that("interior points of returned intervals satisfy AR1 < AR2", {
  m <- find_intervals_03(SQUARE, CONST(1))
  for (j in seq_len(ncol(m))) {
    lo <- m[1, j]
    hi <- m[2, j]
    mid <- if (is.infinite(lo) && is.infinite(hi)) 0
           else if (is.infinite(lo)) hi - 1
           else if (is.infinite(hi)) lo + 1
           else (lo + hi) / 2
    expect_true(ar_value(SQUARE, mid) < ar_value(CONST(1), mid))
  }
})
