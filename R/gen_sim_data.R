################################################################################
## gen_sim_data.R
## Simulate potential-outcome tables and treatment assignments
##
## PURPOSE:
## Helpers for building the data-generating process used in the package's
## simulations. A "potential outcome" table fixes each unit's (Y(1), Y(0)) and
## its compliance type; an assignment mechanism then draws Z, from which the
## observed (Y, D) are realized.
##
## WORKFLOW:
##   gen_outcome_raw_normal()  -> add_effects_to_raw_outcomes[_additive]()  (potential
##                         outcomes + compliance types)
##   gen_data[_additive]()  wraps those two steps.
##   gen_assignment_CR() / gen_assignment_CR_index()  draw a random assignment.
##   gen_data_onesim()  realizes one observed sample from a fixed table.
##
## COMPLIANCE-TYPE CONVENTION:
## Rows are always ordered always-takers, then never-takers, then compliers,
## encoded by the (d1, d0) columns (treatment taken under Z = 1 and Z = 0):
##   always-taker  (d1 = 1, d0 = 1)
##   never-taker   (d1 = 0, d0 = 0)
##   complier      (d1 = 1, d0 = 0)
################################################################################

#' Generate a raw table of potential outcomes and compliance types
#'
#' Builds one row per unit with placeholder potential outcomes (`y0_raw`
#' equals `y1_raw`, i.e. no treatment effect yet -- "raw") and the compliance
#' indicators `d1`/`d0`. Units are ordered always-takers, never-takers, then
#' compliers.
#'
#' @param Nat Number of always-takers.
#' @param Nnt Number of never-takers.
#' @param Nc Number of compliers.
#' @return A numeric matrix with columns `indices`, `y1_raw`, `y0_raw`, `d1`,
#'   `d0` and `Nat + Nnt + Nc` rows.
#' @noRd
gen_outcome_raw_normal=function(Nat,Nnt,Nc){

    # This function generates a table of potential outcomes and complier status
    # Nat: number of always-takers 
    # Nnt: number of never-takers
    # Nc: number of compliers

    # By saying "raw", we mean that we have not attempted to add any effects to the outcomes yet

    # N: total number of units
    N=Nat+Nnt+Nc
    indices=1:N

    # Generate potential outcomes
    y1=stats::rnorm(N,0,1) #note what y1 equals here does not matter. They are going to be overwritten later in the function add_effects_to_raw_outcomes
    y0=y1 
    d1=c(rep(1,Nat),rep(0,Nnt),rep(1,Nc))
    d0=c(rep(1,Nat),rep(0,Nnt),rep(0,Nc))
    

    # Combine into a table
    output = cbind(indices,y1,y0,d1,d0)

    colnames(output)=c('indices',"y1_raw","y0_raw","d1","d0")
    
    return(output) 

}

#' Add a constant (additive) treatment effect to raw outcomes
#'
#' Sets `Y(1) = Y(0) + tau` for every unit by shifting the `y1` column, giving a
#' homogeneous treatment effect equal to `tau`.
#'
#' @param outcome_raw Raw outcome matrix from `gen_outcome_raw_normal()`.
#' @param tau Constant treatment effect applied to all units.
#' @return The outcome matrix with columns renamed to `indices`, `y1`, `y0`,
#'   `d1`, `d0` and `y1` shifted by `tau`.
#' @noRd
add_effects_to_raw_outcomes_additive=function(outcome_raw,tau){

    # This function adds treatment effects to the raw outcomes
    # outcome_raw: table of raw outcomes and complier status
    # tau1: treatment effect for treated units
    # tau0: treatment effect for control units

    # Not we are assuming that in the row outcomes the order is always-takers, never-takers, compliers
    outcome_raw[,2]=outcome_raw[,2] + tau

    #calibrate the effects to match exactly 
    colnames(outcome_raw)=c('indices',"y1","y0","d1","d0")

    return(outcome_raw) 

}

#' Add heterogeneous treatment effects calibrated to exact group means
#'
#' Draws a random treatment effect for each unit (scaled around its group's
#' target `tau`), adds it to `y1`, then re-centers each compliance group so that
#' the group-average treatment effect `mean(y1 - y0)` equals its target `tau`
#' exactly. This yields within-group heterogeneity while pinning the group means.
#'
#' @param outcome_raw Raw outcome matrix from `gen_outcome_raw_normal()`, ordered
#'   always-takers, never-takers, then compliers.
#' @param tau_at,tau_nt,tau_c Target average treatment effect for always-takers,
#'   never-takers, and compliers respectively.
#' @param N_at,N_nt,N_c Number of units in each of those groups.
#' @param mode Either `"heterogeneous"` (default) to draw unit-level effects
#'   `tau * rnorm(n, 1, 0.1)` around each group's target and re-center the group
#'   mean to `tau`, or `"constant"` to give every unit in a group exactly its
#'   `tau` (a constant principal-stratum effect).
#' @return The outcome matrix with columns renamed to `indices`, `y1`, `y0`,
#'   `d1`, `d0` and `y1` updated so each group's mean effect matches its `tau`.
#' @noRd
add_effects_to_raw_outcomes=function(outcome_raw,tau_at,tau_nt,tau_c,N_at,N_nt,N_c,
                                     mode="heterogeneous"){

    # This function adds treatment effects to the raw outcomes
    # outcome_raw: table of raw outcomes and complier status
    # tau_at: treatment effect for always-takers
    # tau_nt: treatment effect for never-takers
    # tau_c: treatment effect for compliers
    # mode: "constant"     -> every unit gets exactly its group's tau
    #       "heterogeneous"-> unit-level effects tau*rnorm(n,1,0.1), mean-calibrated

    mode = match.arg(mode, c("heterogeneous", "constant"))

    # Not we are assuming that in the row outcomes the order is always-takers, never-takers, compliers

    if (mode == "constant"){
        # Constant principal-stratum effect: no within-group heterogeneity.
        eff_at = rep(tau_at, N_at)
        eff_nt = rep(tau_nt, N_nt)
        eff_c  = rep(tau_c,  N_c)
    } else {
        # Heterogeneous effects scaled around each group's target tau.
        # Extract raw outcomes
        if (tau_at!=0){
            eff_at =  tau_at*stats::rnorm(N_at,1,0.1)  #
            #eff_at =  rep(tau_at,N_at)


        }else {
           #eff_at = rnorm(N_at,0,0.1)
           eff_at =  rep(tau_at,N_at)

        }

        if (tau_nt!=0){
            eff_nt =  tau_nt*stats::rnorm(N_nt,1,0.1)
            #eff_nt = rep(tau_nt,N_nt)

        }else {
            #eff_nt = rnorm(N_nt,0.1)
            eff_nt =rep(tau_nt,N_nt)

        }

        if (tau_c!=0){
            eff_c =  tau_c*stats::rnorm(N_c,1,0.1)
            #eff_c = rep(tau_c,N_c)
            #eff_c =  tau_c
        }else {
            eff_c = stats::rnorm(N_c,0,0.1)
            #eff_c = rep(tau_c,N_c)
        }
    }

    outcome_raw[,2]=outcome_raw[,2] + c(eff_at,eff_nt,eff_c)
    #browser()
    #calibrate the effects to match exactly 
    #For always-takers
    factor_at =  mean(outcome_raw[1:N_at,2]-outcome_raw[1:N_at,3]) -tau_at 
    outcome_raw[1:N_at,2] = outcome_raw[1:N_at,2] - factor_at 

    #For never-takers
    factor_nt = mean(outcome_raw[(N_at+1):(N_at+N_nt),2]-outcome_raw[(N_at+1):(N_at+N_nt),3]) -tau_nt 
    outcome_raw[(N_at+1):(N_at+N_nt),2] = outcome_raw[(N_at+1):(N_at+N_nt),2] -factor_nt 

    #For compliers
    factor_c = mean(outcome_raw[(N_at+N_nt+1):(N_at+N_nt+N_c),2]-outcome_raw[(N_at+N_nt+1):(N_at+N_nt+N_c),3]) - tau_c 
    outcome_raw[(N_at+N_nt+1):(N_at+N_nt+N_c),2] =  outcome_raw[(N_at+N_nt+1):(N_at+N_nt+N_c),2] - factor_c
    

    

    
    colnames(outcome_raw)=c('indices',"y1","y0","d1","d0")

    return(outcome_raw) 

}

#' Draw a completely randomized treatment assignment (matrix form)
#'
#' Randomly selects exactly `N1` of the `N1 + N0` units to be treated, without
#' replacement (complete randomization).
#'
#' @param N1 Number of units assigned to treatment.
#' @param N0 Number of units assigned to control.
#' @return A numeric matrix with columns `indices` and `assignment` (a 0/1
#'   treatment indicator), one row per unit.
#' @export
gen_assignment_CR=function(N1,N0){ # nolint

    # Generate assignment to treatment and control groups: completly randomized
    # N1: number of units in treatment group
    # N0: number of units in control group

    # N: total number of units
    N=N1+N0

    treated=sample(1:N,N1,replace=FALSE)

    # indices of units
    indices=1:N

    # Generate assignment
    assignment=rep(0,N)
    assignment[treated]=1

    # Combine into a table
    output_assignment = cbind(indices,assignment)

    colnames(output_assignment)=c('indices',"assignment")

    return(output_assignment)

}


#' Draw a completely randomized treatment assignment (vector form)
#'
#' Like [gen_assignment_CR()] but returns only the assignment vector. The `i`
#' argument is an unused index slot so the function can be called inside
#' `sapply()`/`replicate()`-style loops over simulation draws.
#'
#' @param N1 Number of units assigned to treatment.
#' @param N0 Number of units assigned to control.
#' @param i Loop/replicate index; ignored, present only for use as an `sapply`
#'   iterator.
#' @return A length-`N1 + N0` numeric vector of 0/1 treatment indicators.
#' @noRd
gen_assignment_CR_index=function(N1,N0,i){ # nolint

  # Generate assignment to treatment and control groups: completly randomized
  # N1: number of units in treatment group
  # N0: number of units in control group
  
  # N: total number of units
  N=N1+N0
  
  treated=sample(1:N,N1,replace=FALSE)
  
  # indices of units
  indices=1:N
  
  # Generate assignment
  assignment=rep(0,N)
  assignment[treated]=1
   
  return(assignment)
  
}
#' Generate a potential-outcome table with heterogeneous effects
#'
#' Convenience wrapper: builds a raw table with `gen_outcome_raw_normal()` and applies
#' group-calibrated heterogeneous effects with `add_effects_to_raw_outcomes()`.
#'
#' @param N_at,N_nt,N_c Number of always-takers, never-takers, and compliers.
#' @param tau_at,tau_nt,tau_c Target average treatment effect for each group.
#' @param mode Passed to `add_effects_to_raw_outcomes()`: `"heterogeneous"`
#'   (default) or `"constant"`.
#' @return A potential-outcome matrix with columns `indices`, `y1`, `y0`, `d1`,
#'   `d0`.
#' @noRd
gen_data=function(N_at,N_nt,N_c,tau_at,tau_nt,tau_c,mode="heterogeneous"){

    #generate an outcome table
    outcome_table=gen_outcome_raw_normal(N_at,N_nt,N_c)

    #generate heterogeneous treatment effects
    outcome_table=add_effects_to_raw_outcomes(outcome_table,tau_at,tau_nt,tau_c,N_at,N_nt,N_c,mode)

    return(outcome_table)

}

#' Generate a potential-outcome table with a constant additive effect
#'
#' Convenience wrapper: builds a raw table with `gen_outcome_raw_normal()` and applies a
#' single homogeneous effect `tau` with `add_effects_to_raw_outcomes_additive()`.
#'
#' @param N_at,N_nt,N_c Number of always-takers, never-takers, and compliers.
#' @param tau Constant treatment effect applied to all units.
#' @return A potential-outcome matrix with columns `indices`, `y1`, `y0`, `d1`,
#'   `d0`.
#' @noRd
gen_data_additive=function(N_at,N_nt,N_c,tau){

    #generate an outcome table
    outcome_table=gen_outcome_raw_normal(N_at,N_nt,N_c)

    #generate heterogeneous treatment effects 
    outcome_table=add_effects_to_raw_outcomes_additive(outcome_table,tau)

    return(outcome_table)

}

#' Realize one observed sample from a fixed potential-outcome table
#'
#' Draws one completely randomized assignment, then reveals the observed
#' treatment `D = d1*Z + d0*(1 - Z)` and outcome `Y = y1*D + y0*(1 - D)` implied
#' by each unit's potential outcomes and compliance type.
#'
#' @param outcome_table Potential-outcome matrix (e.g. from `gen_data()`) with
#'   columns `y1`, `y0`, `d1`, `d0`.
#' @param N1 Number of units assigned to treatment (`Z = 1`).
#' @param N0 Number of units assigned to control (`Z = 0`).
#' @return A numeric matrix with columns `Y_observed`, `D_observed`, and
#'   `assignment`, one row per unit.
#' @export
gen_data_onesim = function(outcome_table,N1,N0){


    assignment=gen_assignment_CR(N1,N0)[,2]

    D = outcome_table[,'d1'] * assignment + outcome_table[,'d0']*(1-assignment)

    Y = outcome_table[,'y1'] * D + outcome_table[,'y0'] * (1-D)

    output = cbind(Y,D,assignment)

    colnames(output)=c('Y_observed',"D_observed","assignment")

    return(output)
}

#' Simulate a randomized experiment with noncompliance in one call
#'
#' Builds a fixed potential-outcome ("true") table for a population split into
#' always-reporters, compliers, and never-reporters, then draws one completely
#' randomized experiment from it. By default treatment effects are constant
#' within each principal stratum (principal-strata-constant-effect); set
#' `mode = "heterogeneous"` for unit-level heterogeneity around each stratum's
#' target.
#'
#' @param N Total number of units.
#' @param N1 Number of units assigned to treatment (`Z = 1`); the remaining
#'   `N - N1` units are controls.
#' @param fractions Numeric vector or list `c(fa, fc, fn)` of population
#'   fractions for always-reporters, compliers, and never-reporters. Must sum to
#'   1. Stratum counts are `round(fa * N)` and `round(fc * N)`, with the
#'   never-reporter count taken as the remainder so the three counts sum to `N`
#'   exactly regardless of rounding.
#' @param taus Numeric vector or list `c(tau_a, tau_c, tau_n)` of
#'   principal-stratum treatment effects for always-reporters, compliers, and
#'   never-reporters.
#' @param mode Either `"constant"` (default) to give every unit in a stratum
#'   exactly its `tau`, or `"heterogeneous"` to draw unit-level effects
#'   `tau * rnorm(n, 1, 0.1)` with the stratum mean calibrated to `tau`.
#' @return A list with two elements: `true_table`, the full potential-outcome
#'   table from `gen_data()`; and `observed`, one realized sample from
#'   `gen_data_onesim()` with columns `Y_observed`, `D_observed`, `assignment`.
#' @export
gen_sim_data = function(N, N1, fractions, taus, mode = "constant") {

    # --- Unpack and validate inputs ------------------------------------------
    mode      = match.arg(mode, c("constant", "heterogeneous"))
    fractions = unlist(fractions)
    taus      = unlist(taus)

    if (length(fractions) != 3) {
        stop("`fractions` must have three elements: (fa, fc, fn).")
    }
    if (length(taus) != 3) {
        stop("`taus` must have three elements: (tau_a, tau_c, tau_n).")
    }
    if (N1 < 0 || N1 > N) {
        stop("`N1` must satisfy 0 <= N1 <= N.")
    }

    # Guard: the three fractions must sum to one (up to floating-point error).
    if (abs(sum(fractions) - 1) > 1e-8) {
        stop("`fractions` (fa, fc, fn) must sum to 1; got sum = ",
             sum(fractions), ".")
    }

    # --- Convert fractions to integer counts that sum to N exactly -----------
    # Round the first two strata; the never-reporter count absorbs the rounding
    # remainder so that N_at + N_c + N_nt == N.
    N_at = round(fractions[1] * N)
    N_c  = round(fractions[2] * N)
    N_nt = N - N_at - N_c

    if (N_nt < 0) {
        stop("Rounding produced a negative never-reporter count; ",
             "check `fractions` and `N`.")
    }

    # --- Build the true table and draw one experiment ------------------------
    # gen_data() expects the order (N_at, N_nt, N_c, tau_at, tau_nt, tau_c).
    true_table = gen_data(N_at, N_nt, N_c, taus[1], taus[3], taus[2], mode)

    observed = gen_data_onesim(true_table, N1, N - N1)

    return(list(true_table = true_table, observed = observed))
}
