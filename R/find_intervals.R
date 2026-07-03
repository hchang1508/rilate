## find_intervals_03.R
## Finding Confidence Set Intervals
################################################################################
##
## PURPOSE:
## These functions identify intervals where the observed AR statistic is below
## a permutation AR statistic. This is a key step in constructing confidence
## sets using the Anderson-Rubin test inversion method.
##
## WORKFLOW:
## 1. find_intervals() - Main function that finds acceptance regions
##    - Uses AR_intersection() to find where two AR curves cross
##    - Uses compare_AR() to test which curve is lower in each region
##    - Returns intervals where observed AR < permutation AR
##
## 2. compare_AR() - Helper that evaluates and compares AR values at a point
##
## DEPENDENCIES:
## - Requires AR_intersection() from calculate_intersections_02.R
##
################################################################################


################################################################################
## Main Function: find_intervals
################################################################################

#' Find acceptance intervals where the observed AR statistic is below the
#' permutation AR statistic
#'
#' Solves "for which values of the treatment effect `beta` would we accept the
#' null hypothesis?" The returned intervals form the confidence set from an
#' Anderson-Rubin test inversion. Intersection points of the two AR curves
#' partition the beta-axis into regions; each region is tested at an interior
#' point via `compare_AR()`.
#'
#' @param coef1 Vector of 6 coefficients `[a, b, c, d, e, f]` for the observed
#'   AR function, `AR1(beta) = (a*beta^2 + b*beta + c) /`
#'   `(d*beta^2 + e*beta + f)`.
#' @param coef2 Vector of 6 coefficients for the permutation AR function.
#' @return A 2-by-k matrix whose first row holds the left endpoints and second
#'   row holds the right endpoints of the acceptance intervals. Endpoints may be
#'   `-Inf` or `Inf`. An empty matrix means the null is accepted nowhere.
#' @seealso `AR_intersection()` for the intersection points and `compare_AR()`
#'   for the pointwise comparison.
#' @noRd
find_intervals_03 <- function(coef1, coef2) {
  # Find intervals where AR_observed(β) < AR_permutation(β)
  #
  # This function solves a key question: "For which values of the treatment
  # effect β would we ACCEPT the null hypothesis?" The answer forms our
  # confidence set.
  #
  # INPUTS:
  #   coef1 - Vector of 6 coefficients for observed AR function
  #           Format: [a, b, c, d, e, f] where AR(β) = (aβ² + bβ + c)/(dβ² + eβ + f)
  #   coef2 - Vector of 6 coefficients for permutation AR function
  #
  # OUTPUT:
  #   A 2×k matrix where:
  #     Row 1 = left endpoints of acceptance intervals
  #     Row 2 = right endpoints of acceptance intervals
  #   
  #   Example output:
  #        [,1]  [,2]  [,3]
  #   [1,] -Inf   2.5   8.0    <- left endpoints
  #   [2,]  1.0   5.0   Inf    <- right endpoints
  #   
  #   This means: accept null for β ∈ (-∞, 1.0] ∪ [2.5, 5.0] ∪ [8.0, ∞)
  
  # STEP 1: Find where the two AR curves intersect
  # These intersection points divide the β-axis into regions where one curve
  # is consistently above or below the other
  intersections <- AR_intersection(coef1, coef2)
  
  # Initialize empty vectors to store interval endpoints
  left_end <- c()
  right_end <- c()
  
  # STEP 2: Sort intersections in increasing order
  # This ensures we process regions from left to right on the β-axis
  intersections <- sort(intersections)
  
  # STEP 3: Check if there are any intersections
  if (length(intersections) != 0) {
    
    # There are intersections - we need to test each region between them
    # If there are n intersections, there are n+1 regions to check:
    # Region 1: (-∞, intersection[1])
    # Region 2: (intersection[1], intersection[2])
    # ...
    # Region n+1: (intersection[n], ∞)
    
    for (i in 1:(length(intersections) + 1)) {
      
      if (i == 1) {
        # REGION 1: Before the first intersection
        # Test at a point 1 unit to the left of first intersection
        beta <- intersections[i] - 1
        
        # Check if observed AR is below permutation AR at this point
        ind <- compare_AR(beta, coef1, coef2)
        
        # If TRUE, this entire region is in the acceptance set
        if (ind == TRUE) {
          left_end <- c(left_end, -Inf)
          right_end <- c(right_end, intersections[i])
        }
        
      } else if (i == (length(intersections) + 1)) {
        # LAST REGION: After the final intersection
        # Test at a point 1 unit to the right of last intersection
        beta <- intersections[i - 1] + 1
        
        # Check if observed AR is below permutation AR
        ind <- compare_AR(beta, coef1, coef2)
        
        # If TRUE, this entire region is in the acceptance set
        if (ind == TRUE) {
          left_end <- c(left_end, intersections[i - 1])
          right_end <- c(right_end, Inf)
        }
        
      } else {
        # MIDDLE REGIONS: Between two intersections
        # Test at the midpoint between consecutive intersections
        beta <- (intersections[i - 1] + intersections[i]) / 2
        
        # Check if observed AR is below permutation AR
        ind <- compare_AR(beta, coef1, coef2)
        
        # If TRUE, this region is in the acceptance set
        if (ind == TRUE) {
          left_end <- c(left_end, intersections[i - 1])
          right_end <- c(right_end, intersections[i])
        }
      }
    }
    
  } else {
    # NO INTERSECTIONS: The two AR curves never cross
    # This means one curve is always above the other everywhere
    # We just need to test one point to see which is which
    
    # Test at β = 0 (arbitrary choice - any β would work)
    ind <- compare_AR(0, coef1, coef2)
    
    # If observed is below permutation everywhere, accept for all β
    if (ind == TRUE) {
      left_end <- -Inf
      right_end <- Inf
    }
    # If observed is above permutation everywhere, accept nowhere
    # (left_end and right_end stay empty)
  }
  
  # STEP 4: Package results into a matrix
  # Row 1 = left endpoints, Row 2 = right endpoints
  endpoints <- rbind(left_end, right_end)
  
  return(endpoints)
}


################################################################################
## Helper Function: compare_AR
################################################################################

#' Compare two AR statistics at a specific beta value
#'
#' Evaluates both AR functions at `beta` and reports whether the first
#' (observed) lies below the second (permutation). For `abs(beta) > 1e8` the
#' `1/beta` transformation is used to avoid overflow of the `beta^2` terms.
#'
#' @param beta The treatment effect value at which to compare the two curves.
#' @param coef1 Vector of 6 coefficients `[a, b, c, d, e, f]` for the first
#'   (observed) AR function, `AR(beta) = (a*beta^2 + b*beta + c) /`
#'   `(d*beta^2 + e*beta + f)`.
#' @param coef2 Vector of 6 coefficients for the second (permutation) AR
#'   function.
#' @return `TRUE` if `AR1(beta) < AR2(beta)` (observed below permutation),
#'   otherwise `FALSE`.
#' @noRd
compare_AR <- function(beta, coef1, coef2) {
  # Compare two AR statistics at a specific β value
  #
  # INPUTS:
  #   beta  - The treatment effect value where we want to compare
  #   coef1 - Coefficients [a,b,c,d,e,f] for first AR function (observed)
  #   coef2 - Coefficients [a,b,c,d,e,f] for second AR function (permutation)
  #
  # OUTPUT:
  #   TRUE  if AR1(beta) < AR2(beta)  [observed is below permutation]
  #   FALSE if AR1(beta) >= AR2(beta) [observed is at or above permutation]
  #
  # MATHEMATICAL FORM:
  #   AR(β) = (a*β² + b*β + c) / (d*β² + e*β + f)
  #
  # NUMERICAL STABILITY:
  #   For very large |β|, we use the transformation 1/β to avoid overflow
  
  # Check if beta is extremely large (positive or negative)
  if (abs(beta) > 1e+8) {
    
    # USE 1/β TRANSFORMATION for numerical stability
    # When |β| is huge, β² can overflow. Instead, we rewrite:
    # AR(β) = (a*β² + b*β + c) / (d*β² + e*β + f)
    #       = (a + b/β + c/β²) / (d + e/β + f/β²)
    # by dividing numerator and denominator by β²
    
    # Calculate numerator of first AR: a + b/β + c/β²
    num1 <- sum(coef1[1:3] * c(1, 1/beta, 1/beta^2))

    # Calculate denominator of first AR: d + e/β + f/β²
    dem1 <- sum(coef1[4:6] * c(1, 1/beta, 1/beta^2))

    # Calculate numerator of second AR
    num2 <- sum(coef2[1:3] * c(1, 1/beta, 1/beta^2))

    # Calculate denominator of second AR
    dem2 <- sum(coef2[4:6] * c(1, 1/beta, 1/beta^2))
    
    # Compute AR statistics
    AR1 <- num1 / dem1  # Observed AR
    AR2 <- num2 / dem2  # Permutation AR
    
    # Return TRUE if observed is below permutation
    return(AR1 < AR2)
    
  } else {
    
    # STANDARD EVALUATION for reasonable |β| values
    # Directly compute AR(β) = (a*β² + b*β + c) / (d*β² + e*β + f)
    
    # Calculate numerator of first AR: a*β² + b*β + c
    num1 <- sum(coef1[1:3] * c(beta^2, beta, 1))

    # Calculate denominator of first AR: d*β² + e*β + f
    dem1 <- sum(coef1[4:6] * c(beta^2, beta, 1))

    # Calculate numerator of second AR
    num2 <- sum(coef2[1:3] * c(beta^2, beta, 1))

    # Calculate denominator of second AR
    dem2 <- sum(coef2[4:6] * c(beta^2, beta, 1))
    
    # Compute AR statistics
    AR1 <- num1 / dem1  # Observed AR
    AR2 <- num2 / dem2  # Permutation AR
    
    # Return TRUE if observed is below permutation
    return(AR1 < AR2)
  }
}

