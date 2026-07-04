# ==============================================================================
# AR_helpers.R
# Shared helper functions for the Anderson-Rubin confidence-set algorithms.
#
# These utilities are used by BOTH AR_algo1() (R/AR_algo1.R) and AR_algo2()
# (R/AR_algo2.R). They were previously copy-pasted (with identical bodies) into
# each algorithm file, which risked the two copies silently diverging; they now
# live here as a single source of truth.
#
# `find_quantile_index()` was historically also duplicated under the name
# `find_quantile_R()` in AR_algo2.R -- the two were the same routine, so they are
# unified here under `find_quantile_index()`.
# ==============================================================================

#' Evaluate the AR rational function(s) at a point beta
#'
#' `AR(beta) = (a1 beta^2 + a2 beta + a3) / (a4 beta^2 + a5 beta + a6)`, with the
#' six coefficients stored as rows 1:3 (numerator) and 4:6 (denominator).
#'
#' @param coef Either a length-6 coefficient vector (one AR function) or a
#'   6 x n matrix (one AR function per column).
#' @param beta Point at which to evaluate.
#' @return A scalar when `coef` is a vector; a length-n vector when it is a
#'   matrix.
#' @noRd
ar_value_at <- function(coef, beta) {
  if (is.matrix(coef)) {
    num <- coef[1, ] * beta^2 + coef[2, ] * beta + coef[3, ]
    den <- coef[4, ] * beta^2 + coef[5, ] * beta + coef[6, ]
  } else {
    num <- coef[1] * beta^2 + coef[2] * beta + coef[3]
    den <- coef[4] * beta^2 + coef[5] * beta + coef[6]
  }
  num / den
}

#' Randomization-based (upper-tail) p-value for the AR test at a point beta
#'
#' The Anderson-Rubin statistic is a squared quantity, so the permutation test
#' is one-sided (upper tail). Uses the finite-sample-valid convention that
#' counts the observed assignment as one draw under the null:
#'   p = (1 + #{AR_sim(beta) >= AR_obs(beta)}) / (1 + n_permutations).
#' This is dual to the confidence set built by `find_quantile_at_point()`, which
#' uses the matching (n + 1) critical-value convention.
#'
#' @param beta Null value of the parameter to test (e.g. 0 for LATE = 0).
#' @param AR_obs_coef Length-6 coefficient vector for the observed data.
#' @param AR_sim_coef 6 x n matrix of coefficients, one column per permutation.
#' @return A single p-value in (0, 1].
#' @noRd
ar_pvalue_at <- function(beta, AR_obs_coef, AR_sim_coef) {
  obs <- ar_value_at(AR_obs_coef, beta)
  sim <- ar_value_at(AR_sim_coef, beta)
  (1 + sum(sim >= obs)) / (1 + length(sim))
}

#' Find which permutation gives the alpha quantile in a given interval
#'
#' For the interval indexed by `index` (partitioned by the intersection grid),
#' evaluate the AR functions just inside the interval and return the index of the
#' permutation that sits at the alpha quantile.
#'
#' @param index Which interval (1, 2, 3, ...).
#' @param AR_sim_coef Matrix of AR coefficients (6 x n_permutations).
#' @param intersects Sorted grid of intersection points.
#' @param alpha Significance level (e.g. 0.05 for a 95% CI).
#' @param tol Numerical tolerance.
#' @return Index of the AR function that is the alpha quantile in the interval.
#' @noRd
find_quantile_index <- function(index, AR_sim_coef, intersects, alpha = 0.05,
                                tol = 1e-8) {
  # Find which permutation gives the alpha quantile in a given interval

  if (index == 1) {
    # First interval: (-Inf, intersects[1])
    beta <- intersects[1] - 1
    candidates <- find_quantile_at_point(beta, AR_sim_coef, alpha, tol)
    return(min(candidates))

  } else if (index > 1 && index <= length(intersects)) {
    # Middle interval: (intersects[index-1], intersects[index])
    beta_left <- intersects[index - 1]
    beta_right <- intersects[index]

    # Move inside the interval
    beta_left <- beta_left + 0.1 * (beta_right - beta_left)
    beta_right <- beta_right - 0.1 * (beta_right - beta_left)

    candidates_left <- find_quantile_at_point(beta_left, AR_sim_coef, alpha, tol)
    candidates_right <- find_quantile_at_point(beta_right, AR_sim_coef, alpha, tol)

    # Increase tolerance if no intersection found
    tol_temp <- tol
    while (min(intersect(candidates_left, candidates_right)) == Inf) {
      tol_temp <- tol_temp * 10
      candidates_left <- find_quantile_at_point(beta_left, AR_sim_coef, alpha, tol_temp)
      candidates_right <- find_quantile_at_point(beta_right, AR_sim_coef, alpha, tol_temp)
    }

    return(min(intersect(candidates_left, candidates_right)))

  } else if (index == (length(intersects) + 1)) {
    # Last interval: (intersects[length], Inf)
    beta <- intersects[length(intersects)] + 1
    candidates <- find_quantile_at_point(beta, AR_sim_coef, alpha, tol)
    return(min(candidates))
  }
}

#' Evaluate all AR functions at a point and return those at the alpha quantile
#'
#' The critical value uses the finite-sample-valid (n + 1) convention: instead
#' of the (1 - alpha) quantile of the n permutation values, it takes the
#' ceiling((1 - alpha)(n + 1))-th order statistic (equivalently, the type-1
#' quantile at probability (1 - alpha)(n + 1)/n, capped at 1). This counts the
#' observed assignment as one draw under the null and is dual to the p-value in
#' `ar_pvalue_at()`: at any beta, `AR_obs(beta) <= critical` iff that p-value
#' exceeds alpha.
#'
#' @param beta Point at which to evaluate the AR functions.
#' @param AR_sim_coef Matrix of AR coefficients (6 x n_permutations).
#' @param alpha Significance level.
#' @param tol Numerical tolerance for matching the quantile value.
#' @return Integer indices of the AR functions whose value at `beta` equals the
#'   alpha quantile (within `tol`).
#' @noRd
find_quantile_at_point <- function(beta, AR_sim_coef, alpha, tol = 1e-8) {
  # Evaluate all AR functions at beta and find which ones are at the alpha quantile

  if (abs(beta) > 1e+8) {
    # Use alternative form for large beta
    numerator <- t(AR_sim_coef[1:3, ]) %*% c(1, 1/beta, 1/beta^2)
    denominator <- t(AR_sim_coef[4:6, ]) %*% c(1, 1/beta, 1/beta^2)
  } else {
    # Standard form
    numerator <- t(AR_sim_coef[1:3, ]) %*% c(beta^2, beta, 1)
    denominator <- t(AR_sim_coef[4:6, ]) %*% c(beta^2, beta, 1)
  }

  AR <- numerator / denominator

  # (n + 1) convention: shift the quantile position to include the observed
  # assignment as one draw under the null (cap at 1 so tiny alpha -> full line).
  n_perm <- length(AR)
  p_star <- min((1 - alpha) * (n_perm + 1) / n_perm, 1)
  alpha_quantile <- stats::quantile(AR, p_star, type = 1)

  return(which(abs(AR - alpha_quantile) < tol))
}

#' Select and trim the intervals that fall within a region `[left, right]`
#'
#' @param left,right Region boundaries.
#' @param index Region index into `mapping`.
#' @param intervals List of interval matrices, one per unique quantile index.
#' @param mapping Two-column matrix mapping region index -> quantile index.
#' @param all_indices Unique quantile indices (parallel to `intervals`).
#' @return A list of trimmed `c(left, right)` intervals, or `NULL`.
#' @noRd
select_intervals_for_region <- function(left, right, index, intervals, mapping, all_indices) {
  # Select intervals that fall within [left, right]

  # Find which quantile index corresponds to this region
  quantile_idx <- mapping[mapping[, 1] == index, 2]

  if (length(quantile_idx) == 0) {
    return(NULL)
  }

  # Find which intervals correspond to this quantile
  interval_list_idx <- which(all_indices == quantile_idx)

  if (length(interval_list_idx) == 0) {
    return(NULL)
  }

  # Get the intervals
  selected_intervals <- intervals[[interval_list_idx]]

  if (is.null(selected_intervals) || ncol(selected_intervals) == 0) {
    return(NULL)
  }

  # Filter intervals that overlap with [left, right]
  valid_intervals <- list()
  for (j in 1:ncol(selected_intervals)) {
    int_left <- selected_intervals[1, j]
    int_right <- selected_intervals[2, j]

    # Check if interval overlaps with region
    if (int_right >= left && int_left <= right) {
      # Trim to region boundaries
      trimmed_left <- max(int_left, left)
      trimmed_right <- min(int_right, right)
      valid_intervals[[length(valid_intervals) + 1]] <- c(trimmed_left, trimmed_right)
    }
  }

  return(valid_intervals)
}

#' Convert a list of intervals to a two-column matrix
#'
#' @param interval_list List of length-2 numeric interval vectors.
#' @return A matrix with two columns (empty 0 x 2 matrix if the list is empty).
#' @noRd
unlist_interval_list <- function(interval_list) {
  # Convert list of intervals to matrix
  if (length(interval_list) == 0) {
    return(matrix(nrow = 0, ncol = 2))
  }

  do.call(rbind, interval_list)
}

#' Merge overlapping (or near-adjacent) confidence-set intervals
#'
#' @param CS_region List of interval matrices to combine and merge.
#' @param tol Tolerance within which touching intervals are merged.
#' @return A list of merged `c(left, right)` intervals.
#' @noRd
merge_CS_regions <- function(CS_region, tol = 1e-8) {
  # Merge overlapping intervals

  if (length(CS_region) == 0) {
    return(list())
  }

  # Combine all intervals
  all_intervals <- do.call(rbind, CS_region)

  if (nrow(all_intervals) == 0) {
    return(list())
  }

  # Sort by left endpoint (drop = FALSE keeps a single-row set as a matrix)
  all_intervals <- all_intervals[order(all_intervals[, 1]), , drop = FALSE]

  # Merge overlapping intervals
  merged <- list()
  current <- all_intervals[1, ]

  # seq_len(nrow)[-1] is empty when there is only one interval, avoiding the
  # `2:1` reversal that would otherwise index a non-existent row.
  for (i in seq_len(nrow(all_intervals))[-1]) {
    next_interval <- all_intervals[i, ]

    # Check if intervals overlap or are adjacent
    if (next_interval[1] <= current[2] + tol) {
      # Merge
      current[2] <- max(current[2], next_interval[2])
    } else {
      # No overlap, save current and start new
      merged[[length(merged) + 1]] <- current
      current <- next_interval
    }
  }

  # Add the last interval
  merged[[length(merged) + 1]] <- current

  return(merged)
}
