# ==============================================================================
# B_AR_algo2.R
# Complete Algorithm 2 Implementation (Fast Version with Jumping)
# ==============================================================================
# 
# This implements the Anderson-Rubin permutation test using the FAST algorithm
# that skips intervals by jumping between crossings.
#
# Uses the isolated functions:
#   1. B_solve_coef.R - Calculate AR coefficients
#   2. B_calculate_intersections.R - Find where AR functions cross
#   3. B_find_intervals.R - Find confidence intervals
#
# Key difference from Algorithm 1: Instead of checking EVERY interval,
# Algorithm 2 jumps to the next crossing point (20-25x faster!)
#
# Shared quantile/interval helpers (find_quantile_index, find_quantile_at_point,
# select_intervals_for_region, unlist_interval_list, merge_CS_regions) live in
# R/AR_helpers.R and are shared with AR_algo1(). The only helper unique to
# Algorithm 2 is find_crossings_with_all() (below), which powers the jumping.
#
# ==============================================================================



#' Anderson-Rubin permutation test -- Algorithm 2 (fast, jumping)
#'
#' Computes the same Anderson-Rubin confidence set and randomization p-value as
#' `AR_algo1()`, but instead of checking every interval it jumps between
#' crossing points (via `find_crossings_with_all()`), which is substantially
#' faster. This is the default algorithm. Internal -- dispatched to by
#' `run_rilate()`, not called directly.
#'
#' @param data_table Data frame with columns `Y_observed`, `D_observed`,
#'   `assignment`, and any (demeaned) covariates.
#' @param N1,N0 Number of treated and control units.
#' @param zsim An `N x n_permutations` matrix of permuted assignments.
#' @param tol Numerical tolerance for comparisons.
#' @param alpha Confidence level (default `0.95` for a 95% confidence set).
#' @return A list with `confidence_set` (a list of `[lower, upper]` intervals)
#'   and `p_value` (the randomization p-value for the null `beta = 0`, i.e.
#'   LATE = 0, using the finite-sample-valid `(1 + count)/(1 + n)` convention).
#' @noRd
AR_algo2 <- function(data_table, N1, N0, zsim, tol=1e-8, alpha=0.95) {

  cat("\n========================================\n")
  cat("AR ALGORITHM 2 - FAST IMPLEMENTATION\n")
  cat("========================================\n\n")
  
  nsim_permu <- ncol(zsim)
  cat("Number of permutations:", nsim_permu, "\n")
  cat("Sample size:", nrow(data_table), "\n")
  cat("Treated units:", N1, "\n")
  cat("Control units:", N0, "\n\n")
  
  # ==========================================================================
  # STEP 1: Calculate AR coefficients for observed data
  # ==========================================================================
  
  cat("Step 1: Calculating AR coefficients for observed data...\n")
  AR_obs_coef <- solve_coef_01(data_table, N1, N0)
  cat("  - Observed coefficients computed\n\n")
  
  # ==========================================================================
  # STEP 2: Calculate AR coefficients for all permutations
  # ==========================================================================
  
  cat("Step 2: Calculating AR coefficients for", nsim_permu, "permutations...\n")
  
  AR_sim_coef <- pbapply::pbsapply(1:nsim_permu, function(i) {
    data_temp <- data_table
    data_temp[, "assignment"] <- zsim[, i]
    solve_coef_01(data_temp, N1, N0)
  })
  
  cat("  - Permutation coefficients computed\n\n")

  # ==========================================================================
  # p-value: randomization test of the null beta = 0 (LATE = 0)
  # ==========================================================================
  # AR is upper-tailed; evaluate every AR function at beta = 0 and compare the
  # observed statistic to the permutation distribution.
  p_value <- ar_pvalue_at(0, AR_obs_coef, AR_sim_coef)

  # ==========================================================================
  # STEP 3+4a: Calculate the intersection grid AND the per-curve crossing
  # lookup table in a SINGLE pass.
  # ==========================================================================
  # Previously the grid (Step 3) and the per-curve lookup (Step 4a) each did the
  # pairwise root-finding separately -- the grid over the upper triangle, then a
  # per-curve loop `find_crossings_with_all()` over the FULL square -- so every
  # pairwise `polyroot` was computed ~3 times. calculate_intersections_lut()
  # finds each pair's roots once and distributes them to both owning curves,
  # producing byte-identical `grid` and `lut` with a third of the work. This is
  # the dominant cost of AR_algo2 at large n_permutations.
  cat("Step 3: Calculating intersection grid and lookup table\n")
  ci <- calculate_intersections_lut(AR_sim_coef)
  intersections_grid <- ci$grid
  # Per-curve crossings, stored as GRID INDICES (see calculate_intersections_lut).
  # Working in index space lets the jump loop below locate each crossing on the
  # grid in O(1) instead of searching the ~G-element grid every iteration.
  crossing_idx <- ci$lut_idx

  cat("  - Found", length(intersections_grid), "unique intersection points\n")

  # Check if no intersections
  if (length(intersections_grid) == 0) {
    cat("No intersections found - returning full real line\n")
    return(list(c(-Inf, Inf)))
  }

  # ==========================================================================
  # STEP 4: THE FAST LOOP - Jump between crossings (OPTIMIZED)
  # ==========================================================================

  cat("Step 4: Fast loop - jumping between crossings...\n")
  cat("  - Lookup table built\n")
  
  # Initialize. The threshold is tracked as a GRID INDEX (v_idx) rather than a
  # value: a crossing at grid index k has value grid[k], and grid is ascending,
  # so "value >= grid[v_idx]" is exactly "k >= v_idx". Starting at v_idx = 1
  # keeps all crossings (matches the old threshold v = intersections_grid[1]).
  s <- find_quantile_index(1, AR_sim_coef, intersections_grid, alpha = 1 - alpha, tol = tol)
  s <- as.integer(s)
  v_idx <- 1L                 # Current threshold, as a grid index
  ind_int_old <- 0L           # Track previous index to prevent going backwards

  MaxIter <- length(intersections_grid)
  iter <- 0

  e <- c()        # Endpoints we care about
  indices <- c()  # Quantile indices at each endpoint

  # The jumping loop
  while (iter < MaxIter) {

    # ========================================================================
    # STEP 4.1: Find where CURRENT quantile crosses ALL other AR functions
    # ========================================================================

    # Grid indices where curve s crosses another curve (pre-computed, sorted).
    idx_current <- crossing_idx[[s]]

    # Only keep crossings to the right of the current threshold. Because grid is
    # ascending, filtering by grid index (>= v_idx) is equivalent to the old
    # value filter (crossing value >= grid[v_idx]).
    idx_current <- idx_current[idx_current >= v_idx]

    if (length(idx_current) == 0) {
      # No more crossings - we're done!
      indices <- c(indices, s)
      break
    }

    # ========================================================================
    # STEP 4.2: JUMP to the next crossing
    # ========================================================================
    # The next crossing is the smallest grid index among the remaining crossings
    # (smallest index == smallest value, grid being ascending). We already hold
    # its grid position, so there is NO grid search -- O(1) per iteration, versus
    # the previous min(which(abs(grid - e_i) <= 2*tol)) full-grid scan.
    e_idx <- idx_current[1]                       # crossings are sorted -> min
    e_i   <- intersections_grid[e_idx]            # endpoint value for the CS
    e <- c(e, e_i)
    indices <- c(indices, s)

    # Guard Against Floating Point Error: Move ind_int back to include
    # any grid points that are within tolerance of the current endpoint. This
    # ensures we don't skip over crossings due to tiny numerical differences.
    ind_int <- e_idx
    while (ind_int > 1L &&
           (e_i - intersections_grid[ind_int - 1L]) <= 2 * tol) {
      ind_int <- ind_int - 1L
    }
    ind_int <- max(ind_int, ind_int_old)          # never move backwards

    # Check if we reached the end
    if (ind_int == length(intersections_grid)) {
      indices <- c(indices, s)
      break
    }

    # Jump to the next interval
    v_idx <- ind_int + 1L                         # New threshold (grid index)
    s <- find_quantile_index(ind_int + 1, AR_sim_coef, intersections_grid,
                         alpha = 1 - alpha, tol = tol)
    s <- as.integer(s)

    ind_int_old <- ind_int + 1L  # Update: next ind_int must be >= this

    # Print progress every 100 iterations (like Haoge)
    if (iter %% 100 == 0) {
      cat("  Iteration:", iter, "| Threshold:", round(intersections_grid[v_idx], 4),
          "| Quantile index:", s, "\n")
    }

    iter <- iter + 1
  }
  
  cat("  - Completed in", iter, "iterations (skipped", 
      length(intersections_grid) - iter, "intervals!)\n\n")
  
  # ==========================================================================
  # STEP 5: Find intervals where observed AR < quantile AR
  # ==========================================================================
  
  cat("Step 5: Finding confidence intervals...\n")
  
  # Create mapping
  interval_indice_mapping <- cbind(1:length(indices), indices)
  all_alpha_indices <- unique(indices)
  
  # For each unique quantile, find where observed crosses it
  intervals <- pbapply::pblapply(1:length(all_alpha_indices), function(ind) {
    find_intervals_03(AR_obs_coef, AR_sim_coef[, all_alpha_indices[ind]])
  })
  
  cat("  - Intervals computed\n\n")
  
  # ==========================================================================
  # STEP 6: Select and merge intervals
  # ==========================================================================
  
  cat("Step 6: Constructing confidence set...\n")
  
  # Define endpoints for each region
  left_end_points <- c(-Inf, e)
  right_end_points <- c(e, Inf)
  
  # Select intervals for each region
  CS_list <- pbapply::pblapply(1:length(left_end_points), function(i) {
    select_intervals_for_region(
      left_end_points[i], 
      right_end_points[i], 
      i, 
      intervals, 
      interval_indice_mapping, 
      all_alpha_indices
    )
  })
  
  # Unlist
  CS_region <- lapply(CS_list, unlist_interval_list)
  CS_region <- CS_region[lengths(CS_region) != 0]
  
  cat("  - Confidence regions selected\n\n")
  
  # ==========================================================================
  # STEP 7: Merge overlapping intervals
  # ==========================================================================
  
  cat("Step 7: Merging overlapping intervals...\n")
  CS_region_output <- merge_CS_regions(CS_region, tol)
  cat("  - Final confidence set constructed\n\n")
  
  # ==========================================================================
  # RESULTS
  # ==========================================================================
  
  cat("========================================\n")
  cat("RESULTS\n")
  cat("========================================\n\n")
  
  if (length(CS_region_output) == 0) {
    cat("Confidence set is EMPTY\n\n")
  } else {
    cat(alpha*100, "% Confidence Set:\n", sep="")
    for (i in 1:length(CS_region_output)) {
      interval <- CS_region_output[[i]]
      cat("  [", interval[1], ", ", interval[2], "]\n", sep="")
    }
    cat("\n")
  }

  cat("Randomization p-value (H0: beta = 0):", round(p_value, 4), "\n\n")

  return(list(
    confidence_set = CS_region_output,
    p_value        = p_value
  ))
}

#' Find every crossing of one AR function with all the others
#'
#' Solves `AR_one(beta) = AR_all[, j](beta)` for each permutation `j` and returns
#' the sorted vector of all crossing points. `AR_algo2()` calls this once per AR
#' function to build the lookup table that lets its fast loop jump between
#' crossings instead of scanning every interval.
#'
#' @param AR_one Length-6 coefficient vector for one AR function.
#' @param AR_all A 6 x n matrix of AR coefficients (one function per column).
#' @return A sorted numeric vector of all crossing points.
#' @noRd
find_crossings_with_all <- function(AR_one, AR_all) {

  nsim <- ncol(AR_all)
  
  # Find intersections with each AR function
  intersects <- sapply(1:nsim, function(j) {
    AR_intersection(AR_one, AR_all[, j])
  })
  
  # Unlist and sort
  intersects <- sort(unlist(intersects))

  return(intersects)
}
