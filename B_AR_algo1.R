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
# ==============================================================================

# Load required libraries
library(pbapply)

# Source the 3 component functions
source('B_solve_coef.R')
source('B_calculate_intersections.R')
source('B_find_intervals.R')
# ==============================================================================
# Main Function: AR Algorithm 1
# ==============================================================================

AR_algo1_custom <- function(data_table, N1, N0, zsim, tol=1e-8, alpha=0.95) {
  # Anderson-Rubin Permutation Test - Algorithm 1 (Pure R Implementation)
  #
  # Arguments:
  #   data_table: data frame with columns 'Y_observed', 'D_observed', 'assignment', and covariates
  #   N1: number of treated units
  #   N0: number of control units
  #   zsim: matrix of permuted assignments (N x n_simulations)
  #   tol: tolerance for numerical comparisons
  #   alpha: confidence level (default 0.95 for 95% CI)
  #
  # Returns:
  #   Confidence set as a list of intervals
  
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
  cat("  ✓ Observed coefficients computed\n\n")
  
  # ==========================================================================
  # STEP 2: Calculate AR coefficients for all permutations
  # ==========================================================================
  
  cat("Step 2: Calculating AR coefficients for", nsim_permu, "permutations...\n")
  
  AR_sim_coef <- pbsapply(1:nsim_permu, function(i) {
    # Create temporary data with permuted assignment
    data_temp <- data_table
    data_temp[, "assignment"] <- zsim[, i]
    
    # Calculate coefficients for this permutation
    solve_coef_01(data_temp, N1, N0)
  })
  
  cat("  ✓ Permutation coefficients computed\n\n")
  
  # ==========================================================================
  # STEP 3: Calculate intersections between all pairs of AR functions
  # ==========================================================================
  
  cat("Step 3: Calculating intersections between AR functions...\n")
  intersections_grid <- calculate_intersections_02(AR_sim_coef)
  cat("  ✓ Found", length(intersections_grid), "unique intersection points\n\n")
  
  # ==========================================================================
  # STEP 4: For each interval, find the alpha-quantile AR function
  # ==========================================================================
  
  cat("Step 4: Finding", alpha*100, "% quantile for each interval...\n")
  
  # For each region between intersections, find which permutation gives the alpha quantile
  indices <- pblapply(1:(length(intersections_grid) + 1), function(index) {
    find_quantile_index(index, AR_sim_coef, intersections_grid, 1 - alpha, tol)
  })
  
  cat("  ✓ Quantile indices found\n\n")
  
  # ==========================================================================
  # STEP 5: Find intervals where observed AR < quantile AR
  # ==========================================================================
  
  cat("Step 5: Finding confidence intervals...\n")
  
  # Create mapping between intervals and quantile indices
  interval_indice_mapping <- cbind(1:length(unlist(indices)), unlist(indices))
  
  # Find all unique quantile indices
  all_alpha_indices <- unique(unlist(indices))
  
  # For each unique quantile, find where observed AR crosses it
  intervals <- pblapply(1:length(all_alpha_indices), function(ind) {
    find_intervals_03(AR_obs_coef, AR_sim_coef[, all_alpha_indices[ind]])
  })
  
  cat("  ✓ Intervals computed\n\n")
  
  # ==========================================================================
  # STEP 6: Select and merge confidence intervals
  # ==========================================================================
  
  cat("Step 6: Constructing confidence set...\n")
  
  # Define left and right endpoints for each region
  left_end_points <- c(-Inf, intersections_grid)
  right_end_points <- c(intersections_grid, Inf)
  
  # For each region, select the appropriate intervals
  CS_list <- pblapply(1:length(left_end_points), function(i) {
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
  
  cat("  ✓ Confidence regions selected\n\n")
  
  # ==========================================================================
  # STEP 7: Merge overlapping intervals
  # ==========================================================================
  
  cat("Step 7: Merging overlapping intervals...\n")
  CS_region_output <- merge_CS_regions(CS_region, tol)
  cat("  ✓ Final confidence set constructed\n\n")
  
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
  
  return(CS_region_output)
}

# ==============================================================================
# Helper Functions
# ==============================================================================

find_quantile_index <- function(index, AR_sim_coef, intersects, alpha, tol=1e-8) {
  # Find which permutation gives the alpha quantile in a given interval
  
  if (index == 1) {
    # First interval: (-Inf, intersects[1])
    beta <- intersects[1] - 1
    candidates <- find_quantile_at_point(beta, AR_sim_coef, alpha, tol)
    return(min(candidates))
    
  } else if (index > 1 && index <= length(intersects)) {
    # Middle interval: (intersects[index-1], intersects[index])
    beta_left <- intersects[index-1]
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

find_quantile_at_point <- function(beta, AR_sim_coef, alpha, tol=1e-8) {
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
  alpha_quantile <- quantile(AR, 1 - alpha, type = 1)
  
  return(which(abs(AR - alpha_quantile) < tol))
}

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

unlist_interval_list <- function(interval_list) {
  # Convert list of intervals to matrix
  if (length(interval_list) == 0) {
    return(matrix(nrow=0, ncol=2))
  }
  
  do.call(rbind, interval_list)
}

merge_CS_regions <- function(CS_region, tol=1e-8) {
  # Merge overlapping intervals
  
  if (length(CS_region) == 0) {
    return(list())
  }
  
  # Combine all intervals
  all_intervals <- do.call(rbind, CS_region)
  
  if (nrow(all_intervals) == 0) {
    return(list())
  }
  
  # Sort by left endpoint
  all_intervals <- all_intervals[order(all_intervals[, 1]), ]
  
  # Merge overlapping intervals
  merged <- list()
  current <- all_intervals[1, ]
  
  for (i in 2:nrow(all_intervals)) {
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