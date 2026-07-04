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
# ==============================================================================

# Load required libraries
library(pbapply)

# Source the component functions (assumes working directory is set correctly)
# Load required libraries
library(pbapply)

# Source the component functions with full paths
source('/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run/B_solve_coef.R')
source('/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run/B_calculate_intersections.R')
source('/Users/ag5276/Documents/Github/randomization_noncompliance/arya_test_run/B_find_intervals.R')

# ==============================================================================
# Main Function: AR Algorithm 2 (Fast Version)
# ==============================================================================

AR_algo2_custom <- function(data_table, N1, N0, zsim, tol=1e-8, alpha=0.95) {
  # Anderson-Rubin Permutation Test - Algorithm 2 (Fast with Jumping)
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
  cat("  ✓ Observed coefficients computed\n\n")
  
  # ==========================================================================
  # STEP 2: Calculate AR coefficients for all permutations
  # ==========================================================================
  
  cat("Step 2: Calculating AR coefficients for", nsim_permu, "permutations...\n")
  
  AR_sim_coef <- pbsapply(1:nsim_permu, function(i) {
    data_temp <- data_table
    data_temp[, "assignment"] <- zsim[, i]
    solve_coef_01(data_temp, N1, N0)
  })
  
  cat("  ✓ Permutation coefficients computed\n\n")
  
  # ==========================================================================
  # STEP 3: Calculate ALL intersections (for the grid)
  # ==========================================================================
  
  cat("Step 3: Calculating intersections between AR functions...\n")
  intersections_grid <- calculate_intersections_02(AR_sim_coef)
  
  cat("  ✓ Found", length(intersections_grid), "unique intersection points\n\n")
  
  # Check if no intersections
  if (length(intersections_grid) == 0) {
    cat("No intersections found - returning full real line\n")
    return(list(c(-Inf, Inf)))
  }
  
  # ==========================================================================
  # STEP 4: THE FAST LOOP - Jump between crossings (OPTIMIZED)
  # ==========================================================================
  
  cat("Step 4: Fast loop - jumping between crossings...\n")
  
  # Pre-compute intersection table for fast lookup
  cat("  Building intersection lookup table...\n")
  intersections_with_indices <- list()
  for (i in 1:ncol(AR_sim_coef)) {
    intersections_with_indices[[i]] <- find_crossings_with_all(AR_sim_coef[, i], AR_sim_coef)
  }
  cat("  ✓ Lookup table built\n")
  
  # Initialize
  s <- find_quantile_R(1, AR_sim_coef, intersections_grid, alpha = 1 - alpha, tol = tol)
  s <- as.integer(s)
  v <- intersections_grid[1]  # Current threshold
  ind_int_old <- 0  # Track previous index to prevent going backwards
  
  MaxIter <- length(intersections_grid)
  iter <- 0
  
  e <- c()        # Endpoints we care about
  indices <- c()  # Quantile indices at each endpoint
  
  # The jumping loop
  while (iter < MaxIter) {
    
    # ========================================================================
    # STEP 4.1: Find where CURRENT quantile crosses ALL other AR functions
    # ========================================================================
    
    # Use pre-computed intersections (FAST!)
    int_current <- intersections_with_indices[[s]]
    
    # Only keep intersections to the right of current threshold
    int_current <- int_current[int_current >= v]
    
    if (length(int_current) == 0) {
      # No more crossings - we're done!
      indices <- c(indices, s)
      break
    }
    
    # Find the minimum (next crossing point)
    e_i <- min(int_current)
    e <- c(e, e_i)
    indices <- c(indices, s)
    
    # ========================================================================
    # STEP 4.2: JUMP to the next crossing
    # ========================================================================
    
    # Find where this crossing is on the grid
    ind_int <- min(which(abs(intersections_grid - e_i) <= 2*tol))
    ind_int <- max(ind_int, ind_int_old)  # Ensure we don't go backwards!
    
    # Handle case where no match found (increase tolerance)
    tol_temp <- tol
    while (ind_int == Inf) {
      tol_temp <- tol_temp * 2
      ind_int <- min(which(abs(intersections_grid - e_i) <= tol_temp))
      ind_int <- max(ind_int, ind_int_old)
    }
    
    # Check if we reached the end
    if (ind_int == length(intersections_grid)) {
      indices <- c(indices, s)
      break
    }
    
    # Jump to the next interval
    v <- intersections_grid[ind_int + 1]  # New threshold
    s <- find_quantile_R(ind_int + 1, AR_sim_coef, intersections_grid, 
                         alpha = 1 - alpha, tol = tol)
    s <- as.integer(s)
    
    ind_int_old <- ind_int + 1  # Update: next ind_int must be >= this
    
    # Print progress every 100 iterations (like Haoge)
    if (iter %% 100 == 0) {
      cat("  Iteration:", iter, "| Threshold:", round(v, 4), "| Quantile index:", s, "\n")
    }
    
    iter <- iter + 1
  }
  
  cat("  ✓ Completed in", iter, "iterations (skipped", 
      length(intersections_grid) - iter, "intervals!)\n\n")
  
  # ==========================================================================
  # STEP 5: Find intervals where observed AR < quantile AR
  # ==========================================================================
  
  cat("Step 5: Finding confidence intervals...\n")
  
  # Create mapping
  interval_indice_mapping <- cbind(1:length(indices), indices)
  all_alpha_indices <- unique(indices)
  
  # For each unique quantile, find where observed crosses it
  intervals <- pblapply(1:length(all_alpha_indices), function(ind) {
    find_intervals_03(AR_obs_coef, AR_sim_coef[, all_alpha_indices[ind]])
  })
  
  cat("  ✓ Intervals computed\n\n")
  
  # ==========================================================================
  # STEP 6: Select and merge intervals
  # ==========================================================================
  
  cat("Step 6: Constructing confidence set...\n")
  
  # Define endpoints for each region
  left_end_points <- c(-Inf, e)
  right_end_points <- c(e, Inf)
  
  # Select intervals for each region
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
  
  # Unlist
  CS_region <- lapply(CS_list, unlist_interval_list)
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
# Helper Function: Find crossings of one AR with all others
# ==============================================================================

find_crossings_with_all <- function(AR_one, AR_all) {
  # Find where one specific AR function crosses all other AR functions
  #
  # Arguments:
  #   AR_one: vector of 6 coefficients for one AR function
  #   AR_all: matrix of 6 coefficients for all AR functions (6 x n)
  #
  # Returns:
  #   Sorted vector of all crossing points
  
  nsim <- ncol(AR_all)
  
  # Find intersections with each AR function
  intersects <- sapply(1:nsim, function(j) {
    AR_intersection(AR_one, AR_all[, j])
  })
  
  # Unlist and sort
  intersects <- sort(unlist(intersects))
  
  return(intersects)
}

# ==============================================================================
# Helper Function: Find quantile index
# ==============================================================================

find_quantile_R <- function(index, AR_sim_coef, intersects, alpha = 0.05, tol = 1e-8) {
  # Find which permutation gives the (1-alpha) quantile in a given interval
  #
  # Arguments:
  #   index: which interval (1, 2, 3, ...)
  #   AR_sim_coef: matrix of AR coefficients
  #   intersects: grid of intersection points
  #   alpha: significance level (0.05 for 95% CI)
  #   tol: numerical tolerance
  #
  # Returns:
  #   Index of the AR function that is the (1-alpha) quantile
  
  if (index == 1) {
    # First interval: (-Inf, intersects[1])
    beta <- intersects[1] - 1
    candidates <- find_quantile_at_point(beta, AR_sim_coef, alpha, tol)
    return(min(candidates))
    
  } else if (index > 1 && index <= length(intersects)) {
    # Middle interval
    beta_left <- intersects[index - 1]
    beta_right <- intersects[index]
    beta_left <- beta_left + 0.1 * (beta_right - beta_left)
    beta_right <- beta_right - 0.1 * (beta_right - beta_left)
    
    candidates_left <- find_quantile_at_point(beta_left, AR_sim_coef, alpha, tol)
    candidates_right <- find_quantile_at_point(beta_right, AR_sim_coef, alpha, tol)
    
    tol_temp <- tol
    while (min(intersect(candidates_left, candidates_right)) == Inf) {
      tol_temp <- tol_temp * 10
      candidates_left <- find_quantile_at_point(beta_left, AR_sim_coef, alpha, tol_temp)
      candidates_right <- find_quantile_at_point(beta_right, AR_sim_coef, alpha, tol_temp)
    }
    
    return(min(intersect(candidates_left, candidates_right)))
    
  } else if (index == (length(intersects) + 1)) {
    # Last interval
    beta <- intersects[length(intersects)] + 1
    candidates <- find_quantile_at_point(beta, AR_sim_coef, alpha, tol)
    return(min(candidates))
  }
}

find_quantile_at_point <- function(beta, AR_sim_coef, alpha, tol = 1e-8) {
  # Evaluate all AR functions at beta and find which are at the quantile
  
  if (abs(beta) > 1e+8) {
    numerator <- t(AR_sim_coef[1:3, ]) %*% c(1, 1/beta, 1/beta^2)
    denominator <- t(AR_sim_coef[4:6, ]) %*% c(1, 1/beta, 1/beta^2)
  } else {
    numerator <- t(AR_sim_coef[1:3, ]) %*% c(beta^2, beta, 1)
    denominator <- t(AR_sim_coef[4:6, ]) %*% c(beta^2, beta, 1)
  }
  
  AR <- numerator / denominator
  alpha_quantile <- quantile(AR, 1 - alpha, type = 1)
  
  return(which(abs(AR - alpha_quantile) < tol))
}

# ==============================================================================
# Helper Functions (from Algorithm 1)
# ==============================================================================

select_intervals_for_region <- function(left, right, index, intervals, mapping, all_indices) {
  quantile_idx <- mapping[mapping[, 1] == index, 2]
  if (length(quantile_idx) == 0) return(NULL)
  
  interval_list_idx <- which(all_indices == quantile_idx)
  if (length(interval_list_idx) == 0) return(NULL)
  
  selected_intervals <- intervals[[interval_list_idx]]
  if (is.null(selected_intervals) || ncol(selected_intervals) == 0) return(NULL)
  
  valid_intervals <- list()
  for (j in 1:ncol(selected_intervals)) {
    int_left <- selected_intervals[1, j]
    int_right <- selected_intervals[2, j]
    
    if (int_right >= left && int_left <= right) {
      trimmed_left <- max(int_left, left)
      trimmed_right <- min(int_right, right)
      valid_intervals[[length(valid_intervals) + 1]] <- c(trimmed_left, trimmed_right)
    }
  }
  
  return(valid_intervals)
}

unlist_interval_list <- function(interval_list) {
  if (length(interval_list) == 0) {
    return(matrix(nrow = 0, ncol = 2))
  }
  do.call(rbind, interval_list)
}

merge_CS_regions <- function(CS_region, tol = 1e-8) {
  if (length(CS_region) == 0) return(list())
  
  all_intervals <- do.call(rbind, CS_region)
  if (nrow(all_intervals) == 0) return(list())
  
  all_intervals <- all_intervals[order(all_intervals[, 1]), ]
  
  merged <- list()
  current <- all_intervals[1, ]
  
  for (i in 2:nrow(all_intervals)) {
    next_interval <- all_intervals[i, ]
    
    if (next_interval[1] <= current[2] + tol) {
      current[2] <- max(current[2], next_interval[2])
    } else {
      merged[[length(merged) + 1]] <- current
      current <- next_interval
    }
  }
  
  merged[[length(merged) + 1]] <- current
  return(merged)
}