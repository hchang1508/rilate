# ------------------------------------------------------------------------------
# Unit tests for R/run_rilate.R
#
# NOTE: These unit tests were created by Claude.
#
# run_rilate() is the user-facing wrapper: it standardizes columns, (optionally)
# demeans covariates, guards the per-treatment-group covariate design against
# ill-conditioning, builds zsim, dispatches to the selected AR algorithm, and
# returns each version's confidence set and randomization p-value.
#
# The preparation tests below assert on the wrapper's prepared outputs
# (covariates, inputs, guard, zsim, setup messages); they use a small n_rand so
# the (always-run) algorithm call stays cheap. The "wiring" section checks the
# returned algorithm results directly.
#
# Covariates are OPT-IN via `x`: when `x = NULL` (default) no covariates are
# used and the covariate-adjusted version is skipped; only the columns named in
# `x` are treated as covariates, and any other columns of `data` are ignored.
# ------------------------------------------------------------------------------

# Silence pbapply progress bars (used by the algorithms) so output stays clean.
if (requireNamespace("pbapply", quietly = TRUE)) {
  pbapply::pboptions(type = "none")
}

# Small permutation count to keep the always-run algorithm cheap in prep tests.
NR <- 40

# --- Fixtures ----------------------------------------------------------------

# A small well-conditioned data set with two informative covariates, balanced
# across a completely randomized assignment.
make_df <- function(n = 40, seed = 1) {
  set.seed(seed)
  z <- rep(c(0, 1), each = n / 2)
  data.frame(
    Y_observed = rnorm(n),
    D_observed = rbinom(n, 1, 0.5),
    assignment = z,
    x1 = rnorm(n),
    x2 = rnorm(n)
  )
}

# A realistic experiment with noncompliance and a complier effect, for the
# wiring tests.
make_late_df <- function(N = 60, N1 = 30, seed = 1, tau_c = 2) {
  set.seed(seed)
  sim <- gen_sim_data(N = N, N1 = N1,
                      fractions = c(0.2, 0.6, 0.2),
                      taus = c(0, tau_c, 0))
  as.data.frame(sim$observed)
}

quiet <- function(expr) suppressWarnings(utils::capture.output(res <- expr))

# --- Column mapping ----------------------------------------------------------

test_that("run_rilate maps custom y/d/z names", {
  df <- make_df()
  names(df) <- c("out", "took", "z", "x1", "x2")
  quiet(res <- run_rilate(df, y = "out", d = "took", z = "z", seed = 1, n_rand = NR))
  expect_equal(colnames(res$inputs$without),
               c("Y_observed", "D_observed", "assignment"))
})

test_that("run_rilate errors on a missing mapped column", {
  df <- make_df()
  expect_error(run_rilate(df, y = "nope"), "not found")
})

test_that("run_rilate errors when y/d/z are not distinct", {
  df <- make_df()
  expect_error(run_rilate(df, y = "Y_observed", d = "Y_observed"), "distinct")
})

# --- Input validation --------------------------------------------------------

test_that("run_rilate requires a binary assignment column", {
  df <- make_df()
  df$assignment[1] <- 2
  expect_error(run_rilate(df), "binary")
})

test_that("run_rilate requires both treatment arms to be non-empty", {
  df <- make_df()
  df$assignment <- 1
  expect_error(run_rilate(df), "non-empty")
})

# --- Covariates are opt-in via `x` -------------------------------------------

test_that("x = NULL (default) uses no covariates, even with extra columns", {
  # make_df() has x1/x2 columns, but without `x` they must be ignored.
  quiet(res <- run_rilate(make_df(), seed = 1, n_rand = NR))
  expect_length(res$covariates, 0)
  expect_equal(res$runs, "without")
  expect_null(res$inputs$with)
  expect_null(res$guard)
})

test_that("x names covariates -> both without- and with-covariates versions", {
  quiet(res <- run_rilate(make_df(), x = c("x1", "x2"), seed = 1, n_rand = NR))
  expect_equal(res$covariates, c("x1", "x2"))
  expect_equal(res$runs, c("without", "with"))
  expect_false(is.null(res$inputs$with))
  # with-cov table keeps the covariates; without-cov table drops them
  expect_true(all(c("x1", "x2") %in% colnames(res$inputs$with)))
  expect_false(any(c("x1", "x2") %in% colnames(res$inputs$without)))
})

test_that("x can select a subset of columns as covariates", {
  quiet(res <- run_rilate(make_df(), x = "x1", seed = 1, n_rand = NR))
  expect_equal(res$covariates, "x1")
  expect_true("x1" %in% colnames(res$inputs$with))
  expect_false("x2" %in% colnames(res$inputs$with))
})

test_that("with_covariate = FALSE prepares only the without-cov version", {
  quiet(res <- run_rilate(make_df(), x = c("x1", "x2"),
                      with_covariate = FALSE, seed = 1, n_rand = NR))
  expect_equal(res$runs, "without")
  expect_null(res$inputs$with)
})

# --- `x` validation ----------------------------------------------------------

test_that("x = character(0) is treated as no covariates", {
  quiet(res <- run_rilate(make_df(), x = character(0), seed = 1, n_rand = NR))
  expect_length(res$covariates, 0)
  expect_equal(res$runs, "without")
})

test_that("non-character x is rejected", {
  expect_error(run_rilate(make_df(), x = 1:2), "character vector")
})

test_that("duplicate names in x are rejected", {
  expect_error(run_rilate(make_df(), x = c("x1", "x1")), "duplicate")
})

test_that("x naming a missing column errors", {
  expect_error(run_rilate(make_df(), x = c("x1", "nope")), "not found")
})

test_that("x overlapping the y/d/z columns errors", {
  expect_error(run_rilate(make_df(), x = c("x1", "assignment")), "must not name")
})

# --- Demeaning ---------------------------------------------------------------

test_that("covariates are demeaned to (numerically) zero mean", {
  quiet(res <- run_rilate(make_df(), x = c("x1", "x2"), seed = 1, n_rand = NR))
  expect_equal(mean(res$inputs$with$x1), 0, tolerance = 1e-12)
  expect_equal(mean(res$inputs$with$x2), 0, tolerance = 1e-12)
})

test_that("a demeaning message is printed and names the covariates", {
  out <- utils::capture.output(
    res <- run_rilate(make_df(), x = c("x1", "x2"), seed = 1, n_rand = NR))
  # message uses sep = "" so spacing is tight: "Demeaned2covariate(s): x1, x2"
  expect_true(any(grepl("Demeaned\\s*2\\s*covariate", out)))
  expect_true(any(grepl("x1, x2", out)))
})

test_that("non-numeric covariates are rejected", {
  df <- make_df()
  df$grp <- rep(letters[1:2], length.out = nrow(df))
  expect_error(run_rilate(df, x = "grp"), "numeric")
})

# --- Setup messages ----------------------------------------------------------

test_that("setup message reports N, N1, N0, draws and compliance", {
  out <- utils::capture.output(res <- run_rilate(make_df(n = 40), n_rand = 50,
                                             seed = 1))
  expect_true(any(grepl("Sample size N +: 40", out)))
  expect_true(any(grepl("N1 \\(Z = 1\\) : 20", out)))
  expect_true(any(grepl("N0 \\(Z = 0\\) : 20", out)))
  expect_true(any(grepl("simulated assignment draws +: 50", out)))
  expect_true(any(grepl("Compliance rate", out)))
  expect_true(any(grepl("first stage", out)))
  expect_true(any(grepl("P\\(D=1\\|Z=1\\) - P\\(D=1\\|Z=0\\)", out)))
  expect_true(any(grepl("Will return", out)))
})

# --- Invertibility guard -----------------------------------------------------

test_that("guard stops on a perfectly collinear covariate (rank deficient)", {
  df <- make_df()
  df$x2 <- 2 * df$x1          # exact collinearity within every group
  expect_error(run_rilate(df, x = c("x1", "x2")),
               "rank-deficient|collinear|ill-conditioned")
})

test_that("guard stops on a near-collinear covariate via condition number", {
  df <- make_df()
  df$x2 <- df$x1 + 1e-9 * rnorm(nrow(df))  # nearly collinear
  expect_error(run_rilate(df, x = c("x1", "x2"), cond_threshold = 1e6),
               "ill-conditioned|rank-deficient")
})

test_that("guard records rank and condition number for both groups", {
  quiet(res <- run_rilate(make_df(), x = c("x1", "x2"), seed = 1, n_rand = NR))
  expect_named(res$guard, c("treated", "control"))
  expect_equal(res$guard$treated$rank, 3)   # intercept + x1 + x2
  expect_true(is.finite(res$guard$treated$cond))
})

test_that("guard is skipped (NULL) when no with-covariates run happens", {
  quiet(res <- run_rilate(make_df(), seed = 1, n_rand = NR))  # x = NULL
  expect_null(res$guard)
})

# --- zsim generation ---------------------------------------------------------

test_that("zsim has shape N x n_rand with exactly N1 treated per column", {
  quiet(res <- run_rilate(make_df(n = 40), n_rand = 30, seed = 1))
  expect_equal(dim(res$zsim), c(40, 30))
  expect_true(all(colSums(res$zsim) == res$N1))
  expect_true(all(res$zsim %in% c(0, 1)))
})

test_that("seed makes zsim reproducible", {
  quiet(a <- run_rilate(make_df(), n_rand = 20, seed = 42))
  quiet(b <- run_rilate(make_df(), n_rand = 20, seed = 42))
  expect_identical(a$zsim, b$zsim)
})

# --- Wiring: run_rilate runs the algorithm and returns CS + p_value --------------

test_that("run_rilate returns CS + p_value per version", {
  quiet(res <- run_rilate(make_late_df(), algorithm = "algo2", n_rand = 60,
                      seed = 1))
  expect_named(res$results$without, c("confidence_set", "p_value"))
  expect_true(is.list(res$results$without$confidence_set))
  expect_gt(res$results$without$p_value, 0)
  expect_lte(res$results$without$p_value, 1)
})

test_that("run_rilate result matches calling the algorithm directly", {
  # run_rilate returns the prepared inputs and zsim it used, so a direct call on
  # those must reproduce the embedded result exactly.
  quiet(res <- run_rilate(make_late_df(), algorithm = "algo2", n_rand = 60,
                      seed = 1))
  quiet(direct <- AR_algo2(res$inputs$without, N1 = res$N1, N0 = res$N0,
                           zsim = res$zsim, tol = res$tol, alpha = res$alpha))
  expect_equal(res$results$without$p_value, direct$p_value)
  expect_equal(res$results$without$confidence_set, direct$confidence_set)
})

test_that("run_rilate dispatches to algo1 when requested", {
  quiet(res <- run_rilate(make_late_df(), algorithm = "algo1", n_rand = 60,
                      seed = 1))
  expect_equal(res$algorithm, "algo1")
  expect_named(res$results$without, c("confidence_set", "p_value"))
})

test_that("run_rilate with covariates returns results for both versions", {
  df <- make_late_df()
  df$w <- rnorm(nrow(df))   # a covariate
  quiet(res <- run_rilate(df, x = "w", algorithm = "algo2", n_rand = 60, seed = 1))
  expect_equal(res$runs, c("without", "with"))
  expect_named(res$results$without, c("confidence_set", "p_value"))
  expect_named(res$results$with, c("confidence_set", "p_value"))
})
