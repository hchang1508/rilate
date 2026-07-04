# ------------------------------------------------------------------------------
# Unit tests for R/AR_helpers.R
#
# NOTE: These unit tests were created by Claude.
#
# AR_helpers.R holds the quantile/interval utilities shared by AR_algo1() and
# AR_algo2():
#   find_quantile_index, find_quantile_at_point, select_intervals_for_region,
#   unlist_interval_list, merge_CS_regions
#
# The AR functions are rational in beta:
#   AR_j(beta) = (a1 beta^2 + a2 beta + a3) / (a4 beta^2 + a5 beta + a6)
# stored as a 6 x n coefficient matrix (rows 1:3 numerator, rows 4:6
# denominator). The fixture below makes each AR_j a CONSTANT in beta (numerator
# = a3, denominator = 1) so the quantile arithmetic is fully predictable and
# independent of the evaluation point.
# ------------------------------------------------------------------------------

# --- Fixtures ----------------------------------------------------------------

# Build a 6 x n coefficient matrix whose AR functions are the constants `ar`
# (for every beta, in both the standard and large-|beta| evaluation branches).
const_coef <- function(ar) {
  m <- matrix(0, nrow = 6, ncol = length(ar))
  m[3, ] <- ar   # numerator constant term  -> numerator = ar
  m[6, ] <- 1    # denominator constant term -> denominator = 1
  m
}

# ==============================================================================
# ar_value_at
# ==============================================================================

test_that("ar_value_at evaluates a single AR function (vector coef)", {
  coef <- c(1, 2, 3, 4, 5, 6)   # (b^2 + 2b + 3) / (4b^2 + 5b + 6)
  expect_equal(ar_value_at(coef, 0), 3 / 6)                 # constants only
  expect_equal(ar_value_at(coef, 2), (4 + 4 + 3) / (16 + 10 + 6))  # 11/32
})

test_that("ar_value_at evaluates all AR functions (matrix coef)", {
  coef <- cbind(c(1, 2, 3, 4, 5, 6),   # -> 3/6 = 0.5 at beta = 0
                c(0, 0, 2, 0, 0, 1))   # -> 2/1 = 2   at beta = 0
  expect_equal(ar_value_at(coef, 0), c(0.5, 2))
})

# ==============================================================================
# ar_pvalue_at
# ==============================================================================

test_that("ar_pvalue_at uses the (1 + count)/(1 + n) upper-tail convention", {
  # observed AR(0) = 3; permutation AR(0) values = {1,2,3,4,5}
  obs <- c(0, 0, 3, 0, 0, 1)
  sim <- const_coef(c(1, 2, 3, 4, 5))
  # #{sim >= 3} = 3 (values 3,4,5, tie counted), n = 5 -> (1 + 3)/(1 + 5)
  expect_equal(ar_pvalue_at(0, obs, sim), 4 / 6)
})

test_that("ar_pvalue_at never returns 0 (finite-sample-valid floor)", {
  # observed is strictly more extreme than every permutation -> smallest p
  obs <- c(0, 0, 6, 0, 0, 1)          # AR(0) = 6 > max(sim)
  sim <- const_coef(c(1, 2, 3, 4, 5))
  expect_equal(ar_pvalue_at(0, obs, sim), 1 / 6)   # (1 + 0)/(1 + 5)
  expect_gt(ar_pvalue_at(0, obs, sim), 0)
})

test_that("ar_pvalue_at returns 1 when observed is the least extreme", {
  obs <- c(0, 0, 0, 0, 0, 1)          # AR(0) = 0 <= all sim
  sim <- const_coef(c(1, 2, 3, 4, 5))
  expect_equal(ar_pvalue_at(0, obs, sim), 1)       # (1 + 5)/(1 + 5)
})

# ==============================================================================
# find_quantile_at_point
# ==============================================================================

test_that("find_quantile_at_point returns the index at the (1-alpha) quantile", {
  coef <- const_coef(c(5, 1, 3, 2, 4))
  # alpha = 0.05 -> 0.95 quantile of {1,2,3,4,5} (type 1) = 5 -> AR index 1
  expect_equal(find_quantile_at_point(0, coef, alpha = 0.05), 1L)
})

test_that("find_quantile_at_point is invariant to the evaluation point beta", {
  coef <- const_coef(c(5, 1, 3, 2, 4))
  at_zero  <- find_quantile_at_point(0,    coef, alpha = 0.05)
  at_small <- find_quantile_at_point(3.14, coef, alpha = 0.05)
  # large |beta| exercises the alternative 1/beta evaluation branch
  at_large <- find_quantile_at_point(1e9,  coef, alpha = 0.05)
  expect_equal(at_zero, at_small)
  expect_equal(at_zero, at_large)
})

test_that("find_quantile_at_point returns every tied index at the quantile", {
  coef <- const_coef(c(1, 2, 2, 3))
  # alpha = 0.5 -> median (type 1) of {1,2,2,3} = 2 -> AR indices 2 and 3
  expect_equal(find_quantile_at_point(0, coef, alpha = 0.5), c(2L, 3L))
})

test_that("find_quantile_at_point tolerance controls quantile matching", {
  coef <- const_coef(c(1, 2, 3, 4, 5))
  # An absurdly tiny tol still matches the exact quantile value (no rounding)
  expect_equal(find_quantile_at_point(0, coef, alpha = 0.05, tol = 1e-15), 5L)
})

# ==============================================================================
# find_quantile_index  (interval-aware wrapper around find_quantile_at_point)
# ==============================================================================

test_that("find_quantile_index handles the first, middle and last intervals", {
  coef <- const_coef(c(5, 1, 3, 2, 4))
  intersects <- c(0, 10, 20)          # 3 crossings -> 4 intervals
  # AR is constant in beta, so every interval yields the same quantile index (1)
  expect_equal(find_quantile_index(1, coef, intersects, alpha = 0.05), 1L)  # first
  expect_equal(find_quantile_index(2, coef, intersects, alpha = 0.05), 1L)  # middle
  expect_equal(find_quantile_index(4, coef, intersects, alpha = 0.05), 1L)  # last
})

test_that("find_quantile_index returns the minimum tied index", {
  coef <- const_coef(c(1, 2, 2, 3))    # median tie at indices 2 and 3
  intersects <- c(-1, 1)
  expect_equal(find_quantile_index(1, coef, intersects, alpha = 0.5), 2L)
})

# ==============================================================================
# unlist_interval_list
# ==============================================================================

test_that("unlist_interval_list returns an empty 0 x 2 matrix for an empty list", {
  out <- unlist_interval_list(list())
  expect_true(is.matrix(out))
  expect_equal(dim(out), c(0L, 2L))
})

test_that("unlist_interval_list row-binds interval vectors into a matrix", {
  out <- unlist_interval_list(list(c(0, 1), c(2, 3), c(-5, -4)))
  expect_equal(dim(out), c(3L, 2L))
  expect_equal(out[, 1], c(0, 2, -5))
  expect_equal(out[, 2], c(1, 3, -4))
})

# ==============================================================================
# select_intervals_for_region
# ==============================================================================

# region 1 -> quantile 10, region 2 -> quantile 20, region 3 -> quantile 10
sel_mapping     <- cbind(c(1, 2, 3), c(10, 20, 10))
sel_all_indices <- c(10, 20)
# intervals parallel to sel_all_indices; each is a 2 x k matrix (row1 left,
# row2 right). Quantile 10 owns (0,5) and (8,12); quantile 20 owns (100,110).
sel_intervals <- list(
  matrix(c(0, 5, 8, 12), nrow = 2),   # quantile 10
  matrix(c(100, 110),    nrow = 2)    # quantile 20
)

test_that("select_intervals_for_region returns whole intervals inside the region", {
  out <- select_intervals_for_region(-Inf, Inf, 1, sel_intervals,
                                      sel_mapping, sel_all_indices)
  expect_equal(out, list(c(0, 5), c(8, 12)))
})

test_that("select_intervals_for_region trims intervals to the region boundaries", {
  out <- select_intervals_for_region(3, 10, 1, sel_intervals,
                                      sel_mapping, sel_all_indices)
  expect_equal(out, list(c(3, 5), c(8, 10)))
})

test_that("select_intervals_for_region drops intervals outside the region", {
  out <- select_intervals_for_region(20, 30, 1, sel_intervals,
                                      sel_mapping, sel_all_indices)
  expect_length(out, 0)
})

test_that("select_intervals_for_region follows the region -> quantile mapping", {
  # region 2 maps to quantile 20 -> the (100,110) interval
  out <- select_intervals_for_region(-Inf, Inf, 2, sel_intervals,
                                      sel_mapping, sel_all_indices)
  expect_equal(out, list(c(100, 110)))
})

test_that("select_intervals_for_region returns NULL for an unmapped region", {
  out <- select_intervals_for_region(-Inf, Inf, 99, sel_intervals,
                                      sel_mapping, sel_all_indices)
  expect_null(out)
})

# ==============================================================================
# merge_CS_regions
# ==============================================================================

# helper: one interval as a 1 x 2 matrix, matching unlist_interval_list() output
iv <- function(l, r) matrix(c(l, r), nrow = 1)

test_that("merge_CS_regions returns an empty list for empty input", {
  expect_equal(merge_CS_regions(list()), list())
})

test_that("merge_CS_regions merges overlapping intervals", {
  out <- merge_CS_regions(list(iv(0, 1), iv(0.5, 2)))
  expect_equal(out, list(c(0, 2)))
})

test_that("merge_CS_regions keeps disjoint intervals separate", {
  out <- merge_CS_regions(list(iv(0, 1), iv(3, 4)))
  expect_equal(out, list(c(0, 1), c(3, 4)))
})

test_that("merge_CS_regions sorts unsorted input by left endpoint", {
  out <- merge_CS_regions(list(iv(3, 4), iv(0, 1)))
  expect_equal(out, list(c(0, 1), c(3, 4)))
})

test_that("merge_CS_regions merges intervals that touch within tol", {
  merged     <- merge_CS_regions(list(iv(0, 1), iv(1 + 5e-9, 2)), tol = 1e-8)
  not_merged <- merge_CS_regions(list(iv(0, 1), iv(1 + 5e-9, 2)), tol = 1e-12)
  expect_equal(merged, list(c(0, 2)))
  expect_length(not_merged, 2)
})

test_that("merge_CS_regions returns a single interval unchanged", {
  expect_equal(merge_CS_regions(list(iv(0, 1))), list(c(0, 1)))
})
