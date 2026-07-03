#####################################################################################################
## This script contains functions to calculate intersections between pairs of AR functions, #########
## including the main function `calculate_intersections_02`, the helper function            #########
## `AR_intersection`, and a verification function `check_intersection`.                     #########
#####################################################################################################
#####################################################################################################

#' Calculate all intersection points between pairs of AR functions
#'
#' @param AR_sim_coef Matrix of AR coefficients (6 rows x n_simulations columns).
#' @return Sorted vector of unique intersection points.
#' @noRd
calculate_intersections_02 <- function(AR_sim_coef) {
  nsim_permu <- ncol(AR_sim_coef)

  # Find intersections between all pairs of AR functions
  intersects <- unlist(pbapply::pbsapply(1:(nsim_permu-1), function(i) {
    # For each function i, find intersections with all functions j > i
    sapply((i+1):nsim_permu, function(j) {
      AR_intersection(AR_sim_coef[, i], AR_sim_coef[, j])
    })
  }))
  
  # Return sorted unique intersections
  return(sort(unique(intersects)))
}
################################################################################
## AR_intersection Function
## Extracted from: 2_AR_inversion_algo1.R (lines 249-318)
################################################################################

#' Calculate the intersections of two AR functions
#'
#' @param coef1 Vector of 6 coefficients `[a, b, c, d, e, f]` for the first AR
#'   function, `AR1(beta) = (a*beta^2 + b*beta + c) / (d*beta^2 + e*beta + f)`.
#' @param coef2 Vector of 6 coefficients for the second AR function.
#' @return Vector of real intersection points (beta values where
#'   `AR1(beta) = AR2(beta)`).
#' @noRd
AR_intersection <- function(coef1, coef2){

  # Extract coefficients
  a1 <- coef1[1]
  b1 <- coef1[2]
  c1 <- coef1[3]
  d1 <- coef1[4]
  e1 <- coef1[5]
  f1 <- coef1[6]
  
  a2 <- coef2[1]
  b2 <- coef2[2]
  c2 <- coef2[3]
  d2 <- coef2[4]
  e2 <- coef2[5]
  f2 <- coef2[6]
  
  # Coefficient of beta^4
  v4 <- a1*d2 - a2*d1
  
  # Coefficient of beta^3
  v3 <- a1*e2 - a2*e1 + b1*d2 - b2*d1
  
  # Coefficient of beta^2
  v2 <- a1*f2 - a2*f1 + b1*e2 - b2*e1 + c1*d2 - c2*d1
  
  # Coefficient of beta
  v1 <- b1*f2 - b2*f1 + c1*e2 - c2*e1
  
  # Constant term
  v0 <- c1*f2 - c2*f1
  
  tol <- 1e-8
  
  if (is.na(v0)){
    stop("Constant term v0 is NA; check that coef1/coef2 contain no NA values.")
  }
  
  if (abs(v0)<tol & abs(v1)<tol & abs(v2)<tol & abs(v3)<tol & abs(v4)<tol){
    real_roots <- c()
  } else if(abs(v1)<tol & abs(v2)<tol & abs(v3)<tol &abs(v4)<tol){
    roots <- polyroot(c(v0))
    real_roots <- Re(roots[which(abs(Im(roots))<1e-10)])
  } else if(abs(v2)<tol & abs(v3)<tol &abs(v4)<tol){
    roots <- polyroot(c(v0,v1))
    real_roots <- Re(roots[which(abs(Im(roots))<1e-10)])
  } else if(abs(v3)<tol &abs(v4)<tol){
    roots <- polyroot(c(v0,v1,v2))
    real_roots <- Re(roots[which(abs(Im(roots))<1e-10)])
  } else if(abs(v4)<tol){
    roots <- polyroot(c(v0,v1,v2,v3))
    real_roots <- Re(roots[which(abs(Im(roots))<1e-10)])
  } else {
    roots <- polyroot(c(v0,v1,v2,v3,v4))
    real_roots <- Re(roots[which(abs(Im(roots))<1e-10)])
  }
  
  return(as.vector(real_roots))
}

################################################################################
## Helper Function: Verify an intersection point
################################################################################

#' Verify that a beta value is truly an intersection of two AR functions
#'
#' @param beta The beta value to check.
#' @param coef1 Coefficients for the first AR function.
#' @param coef2 Coefficients for the second AR function.
#' @param tol Numerical tolerance for considering the two values equal.
#' @return A list with `beta`, `AR1`, `AR2`, `difference`, and `is_valid`
#'   (`TRUE` if `beta` is a valid intersection).
#' @noRd
check_intersection <- function(beta, coef1, coef2, tol = 1e-6){

  # Evaluate AR1(beta)
  AR1 <- (coef1[1]*beta^2 + coef1[2]*beta + coef1[3]) / 
    (coef1[4]*beta^2 + coef1[5]*beta + coef1[6])
  
  # Evaluate AR2(beta)
  AR2 <- (coef2[1]*beta^2 + coef2[2]*beta + coef2[3]) / 
    (coef2[4]*beta^2 + coef2[5]*beta + coef2[6])
  
  # Calculate difference
  diff <- abs(AR1 - AR2)
  
  # Check if difference is within tolerance
  is_valid <- diff < tol
  
  return(list(
    beta = beta,
    AR1 = AR1,
    AR2 = AR2,
    difference = diff,
    is_valid = is_valid
  ))
}
