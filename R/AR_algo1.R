# ==============================================================================
# AR_algo1_custom.R
# Complete Algorithm 1 Implementation using Isolated Functions
# ==============================================================================
#
# This implements the Anderson-Rubin permutation test using pure R
# Uses the 3 isolated functions:
#   1. solve_coef_01.R - Calculate AR coefficients
#   2. calculate_intersections_02.R - Find where AR functions cross
#   3. find_intervals_03.R - Find confidence intervals
#
# Shared quantile/interval helpers (find_quantile_index, find_quantile_at_point,
# select_intervals_for_region, unlist_interval_list, merge_CS_regions) live in
# R/AR_helpers.R and are shared with AR_algo2().
#
# ==============================================================================

#' Anderson-Rubin permutation test -- Algorithm 1 (exhaustive, pure R)
#'
#' Computes the Anderson-Rubin confidence set for the LATE and the randomization
#' p-value for the null `beta = 0` by scanning every interval of the permutation
#' grid. This is the reference implementation; `AR_algo2()` returns the same
#' result faster by jumping between crossings. Internal -- dispatched to by
#' `run_AR()`, not called directly.
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
AR_algo1 <- function(data_table, N1, N0, zsim, tol=1e-8, alpha=0.95) {

  cat("\n========================================\n")
  cat("AR ALGORITHM 1 - CUSTOM IMPLEMENTATION\n")
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
    # Create temporary data with permuted assignment
    data_temp <- data_table
    data_temp[, "assignment"] <- zsim[, i]
    
    # Calculate coefficients for this permutation
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
  # STEP 3: Calculate intersections between all pairs of AR functions
  # ==========================================================================
  
  cat("Step 3: Calculating intersections between AR functions...\n")
  intersections_grid <- calculate_intersections_02(AR_sim_coef)
  cat("  - Found", length(intersections_grid), "unique intersection points\n\n")
  
  # ==========================================================================
  # STEP 4: For each interval, find the alpha-quantile AR function
  # ==========================================================================
  
  cat("Step 4: Finding", alpha*100, "% quantile for each interval...\n")
  
  # For each region between intersections, find which permutation gives the alpha quantile
  indices <- pbapply::pblapply(1:(length(intersections_grid) + 1), function(index) {
    find_quantile_index(index, AR_sim_coef, intersections_grid, 1 - alpha, tol)
  })
  
  cat("  - Quantile indices found\n\n")
  
  # ==========================================================================
  # STEP 5: Find intervals where observed AR < quantile AR
  # ==========================================================================
  
  cat("Step 5: Finding confidence intervals...\n")
  
  # Create mapping between intervals and quantile indices
  interval_indice_mapping <- cbind(1:length(unlist(indices)), unlist(indices))
  
  # Find all unique quantile indices
  all_alpha_indices <- unique(unlist(indices))
  
  # For each unique quantile, find where observed AR crosses it
  intervals <- pbapply::pblapply(1:length(all_alpha_indices), function(ind) {
    find_intervals_03(AR_obs_coef, AR_sim_coef[, all_alpha_indices[ind]])
  })
  
  cat("  - Intervals computed\n\n")
  
  # ==========================================================================
  # STEP 6: Select and merge confidence intervals
  # ==========================================================================
  
  cat("Step 6: Constructing confidence set...\n")
  
  # Define left and right endpoints for each region
  left_end_points <- c(-Inf, intersections_grid)
  right_end_points <- c(intersections_grid, Inf)
  
  # For each region, select the appropriate intervals
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
  
  # Unlist intervals
  CS_region <- lapply(CS_list, unlist_interval_list)
  
  # Remove empty intervals
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

# ------------------------------------------------------------------------------
# Helper functions (find_quantile_index, find_quantile_at_point,
# select_intervals_for_region, unlist_interval_list, merge_CS_regions) are
# defined in R/AR_helpers.R and shared with AR_algo2().
# ------------------------------------------------------------------------------