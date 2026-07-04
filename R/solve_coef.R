#' Construct the coefficients of a covariate-adjusted AR rational function
#'
#' Builds the 6 coefficients `[a, b, c, d, e, f]` defining the
#' Anderson-Rubin rational function
#' `AR(beta) = (a*beta^2 + b*beta + c) / (d*beta^2 + e*beta + f)`,
#' using covariate-adjusted outcome and compliance regressions. The numerator
#' coefficients (1-3) are formed from the adjusted treatment effects, and the
#' denominator coefficients (4-6) from the regression residual variances.
#'
#' @param data_table Data frame with columns `Y_observed`, `D_observed`,
#'   `assignment`, and any covariates. `assignment` is the (binary) treatment
#'   assignment.
#' @param N1 Number of treated units.
#' @param N0 Number of control units.
#' @return Numeric vector of 6 coefficients `[a, b, c, d, e, f]` for the AR
#'   rational function.
#' @noRd
solve_coef_01 <- function(data_table, N1, N0) {

  # Initialize coefficient vector
  coef = rep(0, 6)

  #### Imputing ty and td ####
  # Remove D_observed for outcome regressions
  data_tableY = data_table[, -which(names(data_table) == 'D_observed')]
  # Remove Y_observed for compliance regressions
  data_tableD = data_table[, -which(names(data_table) == 'Y_observed')]

  # Split data by treatment assignment
  data_tableY1 = data_tableY[which(data_tableY$assignment == 1), ]  # outcome data for treated group
  data_tableY0 = data_tableY[which(data_tableY$assignment == 0), ]  # outcome data for control group

  data_tableD1 = data_tableD[which(data_tableD$assignment == 1), ]  # compliance data for treated group
  data_tableD0 = data_tableD[which(data_tableD$assignment == 0), ]  # compliance data for control group

  # Run regressions
  regy1 = stats::lm(Y_observed ~ . - assignment, data = data_tableY1)  # outcome regression for treated group
  regy0 = stats::lm(Y_observed ~ . - assignment, data = data_tableY0)  # outcome regression for control group
  regd1 = stats::lm(D_observed ~ . - assignment, data = data_tableD1)  # compliance regression for treated group
  regd0 = stats::lm(D_observed ~ . - assignment, data = data_tableD0)  # compliance regression for control group
  
  # Extract adjusted averages from intercepts
  ty1 = regy1$coef['(Intercept)']  # adjusted average outcome for the treated group
  ty0 = regy0$coef['(Intercept)']  # adjusted average outcome for the control group
  td1 = regd1$coef['(Intercept)']  # adjusted average compliance for the treated group
  td0 = regd0$coef['(Intercept)']  # adjusted average compliance for the control group

  # Calculate treatment effects
  ty = ty1 - ty0
  td = td1 - td0
  
  # Calculate numerator coefficients (coef 1-3)
  coef[1] = td^2
  coef[2] = -2 * ty * td
  coef[3] = ty^2

  #### Imputing variances ####
  # Extract residuals from regressions
  resy1 = regy1$residuals
  resy0 = regy0$residuals
  resd1 = regd1$residuals
  resd0 = regd0$residuals

  # Calculate denominator coefficients (coef 4-6)
  coef[4] = sum(resd1^2) / N1^2 + sum(resd0^2) / N0^2
  coef[5] = -2 * sum(resd1 * resy1) / N1^2 + -2 * sum(resd0 * resy0) / N0^2
  coef[6] = sum(resy1^2) / N1^2 + sum(resy0^2) / N0^2

  return(coef)
}
 
