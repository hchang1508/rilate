# ------------------------------------------------------------------------------
# Unit tests for R/calculate_intersections.R
#
# NOTE: These unit tests were created by Claude.
#
# They exercise the internal (unexported) functions, which testthat makes
# available directly because the package is loaded for the test run:
#   - AR_intersection()            intersections of two AR rational functions
#   - calculate_intersections_02() all pairwise intersections across a matrix
#   - check_intersection()         verify a candidate intersection point
#
# Run automatically via devtools::test() or R CMD check, or on their own with
#   testthat::test_file("tests/testthat/test-calculate_intersections.R")
# ------------------------------------------------------------------------------

# Silence the pbapply progress bar during tests.
pbapply::pboptions(type = "none")

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

# --- AR_intersection ---------------------------------------------------------

test_that("finds the single crossing of a line and a constant (beta = 1)", {
  roots <- AR_intersection(LINEAR(1, 0), CONST(1))  # beta = 1
  expect_equal(sort(roots), 1)
})

test_that("finds both roots of beta^2 = 1", {
  roots <- sort(AR_intersection(SQUARE, CONST(1)))
  expect_equal(roots, c(-1, 1))
})

test_that("finds both roots of beta^2 = beta", {
  roots <- sort(AR_intersection(SQUARE, LINEAR(1, 0)))
  expect_equal(roots, c(0, 1))
})

test_that("identical AR functions yield no isolated intersections", {
  expect_length(AR_intersection(SQUARE, SQUARE), 0)
})

test_that("two distinct constants never intersect", {
  expect_length(AR_intersection(CONST(0), CONST(1)), 0)
})

test_that("returned roots actually satisfy AR1(beta) == AR2(beta)", {
  c1 <- SQUARE
  c2 <- CONST(1)
  for (b in AR_intersection(c1, c2)) {
    expect_equal(ar_value(c1, b), ar_value(c2, b))
  }
})

test_that("errors when the constant term v0 is NA", {
  bad <- c(0, 0, NA, 0, 0, 1)  # NA in the numerator constant -> v0 is NA
  expect_error(AR_intersection(bad, CONST(1)), "NA")
})

# --- check_intersection ------------------------------------------------------

test_that("accepts a true crossing and rejects a false one", {
  true_pt  <- check_intersection(1,   SQUARE, CONST(1))
  false_pt <- check_intersection(0.5, SQUARE, CONST(1))
  expect_true(true_pt$is_valid)
  expect_false(false_pt$is_valid)
  expect_equal(true_pt$AR1, 1)
  expect_equal(true_pt$AR2, 1)
})

# --- calculate_intersections_02 ---------------------------------------------

test_that("returns the sorted, unique union of all pairwise intersections", {
  # Columns: beta^2, the constant 1, and the line beta.
  AR_sim_coef <- cbind(SQUARE, CONST(1), LINEAR(1, 0))
  # Pairwise crossings:
  #   beta^2 = 1     -> {-1, 1}
  #   beta^2 = beta  -> { 0, 1}
  #   1      = beta  -> {    1}
  # Union, sorted & de-duplicated:
  expect_equal(calculate_intersections_02(AR_sim_coef), c(-1, 0, 1))
})

test_that("output is sorted and contains no duplicates", {
  AR_sim_coef <- cbind(SQUARE, CONST(1), LINEAR(1, 0))
  result <- calculate_intersections_02(AR_sim_coef)
  expect_false(is.unsorted(result))
  expect_equal(anyDuplicated(result), 0L)
})
