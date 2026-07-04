# ------------------------------------------------------------------------------
# Unit tests for R/gen_sim_data.R
#
# NOTE: These unit tests were created by Claude.
#
# Exercises the data-generating helpers. Rows are always ordered
# always-takers (AT), then never-takers (NT), then compliers (C), with
# compliance encoded by (d1, d0):
#   AT (1, 1)   NT (0, 0)   C (1, 0)
#
# The tests lean on invariants that hold for ANY input regardless of the
# random draws:
#   * gen_outcome_raw_normal: y0_raw == y1_raw (no effect yet) and the (d1, d0)
#     columns match the fixed compliance-type ordering.
#   * add_effects_to_raw_outcomes: each group's mean effect mean(y1 - y0) is
#     calibrated to equal its target tau EXACTLY (both modes); "constant" mode
#     additionally makes every unit-level effect exactly tau.
#   * gen_assignment_CR: exactly N1 of N units are treated, values are 0/1.
#   * gen_data_onesim: D = d1*Z + d0*(1-Z) and Y = y1*D + y0*(1-D), so AT always
#     take treatment, NT never do, and compliers take it iff assigned.
#
# Run via devtools::test() / R CMD check, or on their own with
#   testthat::test_file("tests/testthat/test-gen_sim_data.R")
# ------------------------------------------------------------------------------

# --- gen_outcome_raw_normal --------------------------------------------------

test_that("gen_outcome_raw_normal has the right shape and column names", {
  out <- gen_outcome_raw_normal(2, 3, 4)
  expect_equal(nrow(out), 9)
  expect_equal(colnames(out), c("indices", "y1_raw", "y0_raw", "d1", "d0"))
  expect_equal(out[, "indices"], 1:9)
})

test_that("gen_outcome_raw_normal starts with no treatment effect (y0 == y1)", {
  out <- gen_outcome_raw_normal(2, 3, 4)
  expect_identical(out[, "y1_raw"], out[, "y0_raw"])
})

test_that("gen_outcome_raw_normal encodes AT/NT/C in the fixed row order", {
  Nat <- 2; Nnt <- 3; Nc <- 4
  out <- gen_outcome_raw_normal(Nat, Nnt, Nc)
  expect_equal(out[, "d1"], c(rep(1, Nat), rep(0, Nnt), rep(1, Nc)))
  expect_equal(out[, "d0"], c(rep(1, Nat), rep(0, Nnt), rep(0, Nc)))
})

# --- add_effects_to_raw_outcomes_additive ------------------------------------

test_that("additive effect sets y1 - y0 == tau for every unit", {
  raw <- gen_outcome_raw_normal(2, 3, 4)
  tau <- 1.5
  out <- add_effects_to_raw_outcomes_additive(raw, tau)
  expect_equal(colnames(out), c("indices", "y1", "y0", "d1", "d0"))
  expect_equal(out[, "y1"] - out[, "y0"], rep(tau, 9))
})

test_that("additive effect leaves y0 and compliance columns untouched", {
  raw <- gen_outcome_raw_normal(2, 3, 4)
  out <- add_effects_to_raw_outcomes_additive(raw, 2)
  expect_identical(out[, 3], raw[, 3])   # y0 unchanged
  expect_identical(out[, 4], raw[, 4])   # d1 unchanged
  expect_identical(out[, 5], raw[, 5])   # d0 unchanged
})

# --- add_effects_to_raw_outcomes ---------------------------------------------

test_that("heterogeneous mode calibrates each group's mean effect to its tau", {
  set.seed(1)
  Nat <- 30; Nnt <- 40; Nc <- 50
  raw <- gen_outcome_raw_normal(Nat, Nnt, Nc)
  out <- add_effects_to_raw_outcomes(raw, 1, -2, 3, Nat, Nnt, Nc,
                                     mode = "heterogeneous")
  eff <- out[, "y1"] - out[, "y0"]
  at <- seq_len(Nat)
  nt <- (Nat + 1):(Nat + Nnt)
  c_ <- (Nat + Nnt + 1):(Nat + Nnt + Nc)
  expect_equal(mean(eff[at]), 1, tolerance = 1e-8)
  expect_equal(mean(eff[nt]), -2, tolerance = 1e-8)
  expect_equal(mean(eff[c_]), 3, tolerance = 1e-8)
})

test_that("heterogeneous mode actually varies effects within a group", {
  set.seed(2)
  Nat <- 30; Nnt <- 40; Nc <- 50
  raw <- gen_outcome_raw_normal(Nat, Nnt, Nc)
  out <- add_effects_to_raw_outcomes(raw, 1, 2, 3, Nat, Nnt, Nc,
                                     mode = "heterogeneous")
  eff_at <- (out[, "y1"] - out[, "y0"])[seq_len(Nat)]
  expect_gt(stats::sd(eff_at), 0)
})

test_that("constant mode gives every unit exactly its group's tau", {
  set.seed(3)
  Nat <- 5; Nnt <- 6; Nc <- 7
  raw <- gen_outcome_raw_normal(Nat, Nnt, Nc)
  out <- add_effects_to_raw_outcomes(raw, 1, -2, 3, Nat, Nnt, Nc,
                                     mode = "constant")
  eff <- out[, "y1"] - out[, "y0"]
  expect_equal(eff, c(rep(1, Nat), rep(-2, Nnt), rep(3, Nc)),
               tolerance = 1e-8)
})

test_that("a zero group tau calibrates that group's mean effect to zero", {
  set.seed(4)
  Nat <- 20; Nnt <- 20; Nc <- 20
  raw <- gen_outcome_raw_normal(Nat, Nnt, Nc)
  out <- add_effects_to_raw_outcomes(raw, 0, 0, 5, Nat, Nnt, Nc,
                                     mode = "heterogeneous")
  eff <- out[, "y1"] - out[, "y0"]
  expect_equal(mean(eff[seq_len(Nat)]), 0, tolerance = 1e-8)
  expect_equal(mean(eff[(Nat + 1):(Nat + Nnt)]), 0, tolerance = 1e-8)
})

test_that("add_effects_to_raw_outcomes rejects an unknown mode", {
  raw <- gen_outcome_raw_normal(2, 2, 2)
  expect_error(
    add_effects_to_raw_outcomes(raw, 1, 1, 1, 2, 2, 2, mode = "bogus")
  )
})

# --- gen_assignment_CR / gen_assignment_CR_index -----------------------------

test_that("gen_assignment_CR treats exactly N1 of N units", {
  set.seed(5)
  N1 <- 7; N0 <- 13
  a <- gen_assignment_CR(N1, N0)
  expect_equal(colnames(a), c("indices", "assignment"))
  expect_equal(nrow(a), N1 + N0)
  expect_equal(sum(a[, "assignment"]), N1)
  expect_true(all(a[, "assignment"] %in% c(0, 1)))
})

test_that("gen_assignment_CR_index returns a plain 0/1 vector of length N", {
  set.seed(6)
  N1 <- 4; N0 <- 6
  v <- gen_assignment_CR_index(N1, N0, i = 99)
  expect_null(dim(v))
  expect_length(v, N1 + N0)
  expect_equal(sum(v), N1)
  expect_true(all(v %in% c(0, 1)))
})

test_that("gen_assignment_CR handles the all-treated / all-control extremes", {
  set.seed(7)
  expect_equal(sum(gen_assignment_CR(5, 0)[, "assignment"]), 5)
  expect_equal(sum(gen_assignment_CR(0, 5)[, "assignment"]), 0)
})

# --- gen_data / gen_data_additive wrappers -----------------------------------

test_that("gen_data returns a calibrated table with the y1/y0 column names", {
  set.seed(8)
  tbl <- gen_data(10, 10, 10, 1, 2, 3, mode = "heterogeneous")
  expect_equal(colnames(tbl), c("indices", "y1", "y0", "d1", "d0"))
  eff <- tbl[, "y1"] - tbl[, "y0"]
  expect_equal(mean(eff[1:10]), 1, tolerance = 1e-8)         # always-takers
  expect_equal(mean(eff[11:20]), 2, tolerance = 1e-8)        # never-takers
  expect_equal(mean(eff[21:30]), 3, tolerance = 1e-8)        # compliers
})

test_that("gen_data_additive applies one homogeneous effect to all units", {
  set.seed(9)
  tbl <- gen_data_additive(4, 5, 6, 2.5)
  expect_equal(tbl[, "y1"] - tbl[, "y0"], rep(2.5, 15))
})

# --- gen_data_onesim ---------------------------------------------------------

test_that("gen_data_onesim realizes D and Y from the potential-outcome rules", {
  set.seed(10)
  # One unit of each type with distinct y1/y0 so mis-selection is detectable.
  #                     indices  y1  y0  d1 d0
  tbl <- rbind(
    c(1, 10, 0, 1, 1),   # always-taker: D == 1 always -> Y == y1 == 10
    c(2, 20, 5, 0, 0),   # never-taker:  D == 0 always -> Y == y0 == 5
    c(3, 30, 7, 1, 0)    # complier:     D == Z        -> Y depends on Z
  )
  colnames(tbl) <- c("indices", "y1", "y0", "d1", "d0")

  sim <- gen_data_onesim(tbl, N1 = 2, N0 = 1)
  expect_equal(colnames(sim), c("Y_observed", "D_observed", "assignment"))

  Z <- sim[, "assignment"]
  D <- sim[, "D_observed"]
  Y <- sim[, "Y_observed"]

  # Compliance rules independent of the random assignment.
  expect_equal(D[1], 1)     # always-taker
  expect_equal(D[2], 0)     # never-taker
  expect_equal(D[3], Z[3])  # complier takes treatment iff assigned

  # Observed outcome selects the potential outcome matching realized D.
  expect_equal(Y[1], 10)
  expect_equal(Y[2], 5)
  expect_equal(Y[3], ifelse(D[3] == 1, 30, 7))
})

# --- gen_sim_data (exported entry point) -------------------------------------

test_that("gen_sim_data splits N into strata that sum to N exactly", {
  set.seed(11)
  N <- 101  # odd + non-round fractions to stress the rounding remainder
  res <- gen_sim_data(N, N1 = 50,
                      fractions = c(0.3, 0.3, 0.4),
                      taus = c(1, 2, 3))
  expect_named(res, c("true_table", "observed"))
  expect_equal(nrow(res$true_table), N)
  expect_equal(nrow(res$observed), N)
})

test_that("gen_sim_data maps fractions/taus onto the AT/NT/C ordering", {
  set.seed(12)
  N <- 300
  fractions <- c(0.2, 0.5, 0.3)   # (fa, fc, fn)
  taus <- c(1, 2, 3)              # (tau_a, tau_c, tau_n)
  res <- gen_sim_data(N, N1 = 150, fractions = fractions, taus = taus,
                      mode = "constant")

  N_at <- round(fractions[1] * N)   # 60 always-takers
  N_c  <- round(fractions[2] * N)   # 150 compliers
  eff  <- res$true_table[, "y1"] - res$true_table[, "y0"]

  # Rows are AT, then NT, then C. tau_a applies to AT, tau_n to NT, tau_c to C.
  expect_equal(eff[1], 1)                        # always-taker gets tau_a
  expect_equal(eff[N_at + 1], 3)                 # never-taker gets tau_n
  expect_equal(eff[N - N_c + 1], 2)              # complier gets tau_c
})

test_that("gen_sim_data validates its inputs", {
  expect_error(gen_sim_data(100, 50, fractions = c(0.5, 0.5),
                            taus = c(1, 2, 3)),
               "three elements")
  expect_error(gen_sim_data(100, 50, fractions = c(0.3, 0.3, 0.4),
                            taus = c(1, 2)),
               "three elements")
  expect_error(gen_sim_data(100, 150, fractions = c(0.3, 0.3, 0.4),
                            taus = c(1, 2, 3)),
               "N1")
  expect_error(gen_sim_data(100, 50, fractions = c(0.3, 0.3, 0.3),
                            taus = c(1, 2, 3)),
               "sum to 1")
})

test_that("gen_sim_data observed treatment respects compliance types", {
  set.seed(13)
  N <- 200
  res <- gen_sim_data(N, N1 = 100,
                      fractions = c(0.25, 0.5, 0.25),
                      taus = c(1, 2, 3), mode = "constant")
  d1 <- res$true_table[, "d1"]
  d0 <- res$true_table[, "d0"]
  Z  <- res$observed[, "assignment"]
  D  <- res$observed[, "D_observed"]
  expect_equal(D, d1 * Z + d0 * (1 - Z))
})
