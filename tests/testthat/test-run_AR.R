# ------------------------------------------------------------------------------
# Unit tests for R/run_AR.R
#
# NOTE: These unit tests were created by Claude.
#
# run_AR() is the user-facing wrapper: it standardizes columns, demeans
# covariates, guards the per-treatment-group covariate design against
# ill-conditioning, builds zsim, and dispatches to the AR algorithms behind a
# gate (execute = FALSE). These tests exercise the preparation and gating logic
# WITHOUT invoking the (not-yet-integrated) algorithms.
# ------------------------------------------------------------------------------

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

quiet <- function(expr) suppressWarnings(utils::capture.output(res <- expr))

# --- Column mapping ----------------------------------------------------------

test_that("run_AR maps custom y/d/z names and treats the rest as covariates", {
  df <- make_df()
  names(df) <- c("out", "took", "z", "x1", "x2")
  quiet(res <- run_AR(df, y = "out", d = "took", z = "z", seed = 1))
  expect_equal(res$covariates, c("x1", "x2"))
  expect_equal(colnames(res$inputs$without),
               c("Y_observed", "D_observed", "assignment"))
})

test_that("run_AR errors on a missing mapped column", {
  df <- make_df()
  expect_error(run_AR(df, y = "nope"), "not found")
})

test_that("run_AR errors when y/d/z are not distinct", {
  df <- make_df()
  expect_error(run_AR(df, y = "Y_observed", d = "Y_observed"), "distinct")
})

# --- Input validation --------------------------------------------------------

test_that("run_AR requires a binary assignment column", {
  df <- make_df()
  df$assignment[1] <- 2
  expect_error(run_AR(df), "binary")
})

test_that("run_AR requires both treatment arms to be non-empty", {
  df <- make_df()
  df$assignment <- 1
  expect_error(run_AR(df), "non-empty")
})

# --- Which versions are prepared ---------------------------------------------

test_that("no covariates -> only the without-covariates version is prepared", {
  df <- make_df()[, c("Y_observed", "D_observed", "assignment")]
  quiet(res <- run_AR(df, seed = 1))
  expect_equal(res$runs, "without")
  expect_null(res$inputs$with)
  expect_length(res$covariates, 0)
})

test_that("covariates + default with_covariate=TRUE -> both versions", {
  quiet(res <- run_AR(make_df(), seed = 1))
  expect_equal(res$runs, c("without", "with"))
  expect_false(is.null(res$inputs$with))
  # with-cov table keeps the covariates; without-cov table drops them
  expect_true(all(c("x1", "x2") %in% colnames(res$inputs$with)))
  expect_false(any(c("x1", "x2") %in% colnames(res$inputs$without)))
})

test_that("with_covariate=FALSE prepares only the without-covariates version", {
  quiet(res <- run_AR(make_df(), with_covariate = FALSE, seed = 1))
  expect_equal(res$runs, "without")
  expect_null(res$inputs$with)
})

# --- Demeaning ---------------------------------------------------------------

test_that("covariates are demeaned to (numerically) zero mean", {
  quiet(res <- run_AR(make_df(), seed = 1))
  expect_equal(mean(res$inputs$with$x1), 0, tolerance = 1e-12)
  expect_equal(mean(res$inputs$with$x2), 0, tolerance = 1e-12)
})

test_that("a demeaning message is printed and names the covariates", {
  out <- utils::capture.output(res <- run_AR(make_df(), seed = 1))
  expect_true(any(grepl("Demeaned 2 covariate", out)))
  expect_true(any(grepl("x1, x2", out)))
})

test_that("non-numeric covariates are rejected", {
  df <- make_df()
  df$grp <- rep(letters[1:2], length.out = nrow(df))
  expect_error(run_AR(df), "numeric")
})

# --- Setup messages ----------------------------------------------------------

test_that("setup message reports N, N1, N0, randomizations and compliance", {
  out <- utils::capture.output(res <- run_AR(make_df(n = 40), n_rand = 250,
                                             seed = 1))
  expect_true(any(grepl("Sample size N +: 40", out)))
  expect_true(any(grepl("N1 \\(Z = 1\\) : 20", out)))
  expect_true(any(grepl("N0 \\(Z = 0\\) : 20", out)))
  expect_true(any(grepl("Randomizations +: 250", out)))
  expect_true(any(grepl("Compliance rate", out)))
  expect_true(any(grepl("Will return", out)))
})

# --- Invertibility guard -----------------------------------------------------

test_that("guard stops on a perfectly collinear covariate (rank deficient)", {
  df <- make_df()
  df$x2 <- 2 * df$x1          # exact collinearity within every group
  expect_error(run_AR(df), "rank-deficient|collinear|ill-conditioned")
})

test_that("guard stops on a near-collinear covariate via condition number", {
  df <- make_df()
  df$x2 <- df$x1 + 1e-9 * rnorm(nrow(df))  # nearly collinear
  expect_error(run_AR(df, cond_threshold = 1e6),
               "ill-conditioned|rank-deficient")
})

test_that("guard records rank and condition number for both groups", {
  quiet(res <- run_AR(make_df(), seed = 1))
  expect_named(res$guard, c("treated", "control"))
  expect_equal(res$guard$treated$rank, 3)   # intercept + x1 + x2
  expect_true(is.finite(res$guard$treated$cond))
})

test_that("guard is skipped (NULL) when no with-covariates run happens", {
  quiet(res <- run_AR(make_df(), with_covariate = FALSE, seed = 1))
  expect_null(res$guard)
})

# --- zsim generation ---------------------------------------------------------

test_that("zsim has shape N x n_rand with exactly N1 treated per column", {
  quiet(res <- run_AR(make_df(n = 40), n_rand = 30, seed = 1))
  expect_equal(dim(res$zsim), c(40, 30))
  expect_true(all(colSums(res$zsim) == res$N1))
  expect_true(all(res$zsim %in% c(0, 1)))
})

test_that("seed makes zsim reproducible", {
  quiet(a <- run_AR(make_df(), n_rand = 20, seed = 42))
  quiet(b <- run_AR(make_df(), n_rand = 20, seed = 42))
  expect_identical(a$zsim, b$zsim)
})

# --- Gate --------------------------------------------------------------------

test_that("execute = FALSE returns prepared inputs without running algorithms", {
  quiet(res <- run_AR(make_df(), algorithm = "algo2", seed = 1))
  expect_equal(res$algorithm, "algo2")
  expect_true(all(vapply(res$results, is.null, logical(1))))
})

test_that("gate prints a not-executed notice naming the algorithm", {
  out <- utils::capture.output(res <- run_AR(make_df(), algorithm = "algo1",
                                             seed = 1))
  expect_true(any(grepl("\\[GATE\\].*AR_algo1_custom", out)))
  expect_true(any(grepl("not executed", out)))
})

test_that("execute = TRUE errors clearly when the algorithm is not sourced", {
  # AR_algo2_custom is not part of the package namespace, so the gate should
  # refuse to run rather than fail obscurely.
  if (!exists("AR_algo2_custom", mode = "function")) {
    expect_error(
      quiet(run_AR(make_df(), execute = TRUE, seed = 1)),
      "not found"
    )
  } else {
    succeed("AR_algo2_custom is available; skipping not-sourced check.")
  }
})
