################################################################################
## run_AR.R
## User-facing wrapper around the Anderson-Rubin confidence-set algorithms
## (AR_algo1_custom / AR_algo2_custom).
##
## Responsibilities:
##   * accept a tidy data frame with user-named y / d / z columns (+ covariates)
##   * standardize column names to what the algorithms expect
##   * demean covariates
##   * guard against an ill-conditioned (near-collinear) covariate design within
##     either treatment group -- solve_coef_01() fits a separate regression per
##     group, so [1, X] must be well-conditioned in EACH group
##   * generate the permuted-assignment matrix (zsim)
##   * dispatch to Algorithm 1 or 2 -- currently behind a GATE (`execute`), since
##     the algorithms are not yet wired into the package namespace
################################################################################

#' Run the Anderson-Rubin confidence-set procedure (wrapper)
#'
#' Convenience wrapper that prepares data for and dispatches to the
#' Anderson-Rubin permutation algorithms (`AR_algo1_custom` /
#' `AR_algo2_custom`). It standardizes column names, demeans covariates, checks
#' the per-treatment-group covariate design for ill-conditioning, and builds the
#' matrix of permuted assignments.
#'
#' When covariates are supplied and `with_covariate = TRUE` (the default), both a
#' without-covariates and a with-covariates analysis are prepared. Set
#' `with_covariate = FALSE`, or supply no covariates, to prepare only the
#' without-covariates analysis.
#'
#' The actual algorithm call is gated by `execute`. With `execute = FALSE`
#' (default) the function performs all preparation and validation and returns the
#' assembled inputs without running the (not-yet-integrated) algorithms.
#'
#' @param data Data frame containing the outcome, treatment-received, assignment,
#'   and (optionally) covariate columns.
#' @param y,d,z Names of the outcome, treatment-received, and assignment columns
#'   in `data`. Default to `"Y_observed"`, `"D_observed"`, `"assignment"`. Every
#'   other column of `data` is treated as a covariate.
#' @param algorithm Which algorithm to use: `"algo2"` (default, fast) or
#'   `"algo1"`.
#' @param n_rand Number of randomizations (permuted assignments) to draw.
#'   Default `1000`.
#' @param with_covariate If `TRUE` (default) and covariates are present, prepare
#'   both the with- and without-covariates analyses. If `FALSE`, prepare only the
#'   without-covariates analysis.
#' @param alpha Confidence level (default `0.95`).
#' @param tol Numerical tolerance passed to the algorithm (default `1e-8`).
#' @param cond_threshold Condition-number cutoff above which a treatment group's
#'   covariate design `[1, X]` is treated as ill-conditioned (default `1e10`).
#' @param seed Optional integer seed for reproducible `zsim` generation.
#' @param execute If `TRUE`, actually invoke the selected algorithm; if `FALSE`
#'   (default), stop after preparation and return the assembled inputs (the
#'   "gate").
#' @return A list with the prepared inputs and metadata: `algorithm`, `alpha`,
#'   `tol`, `N`, `N1`, `N0`, `runs` (which versions were prepared), `covariates`,
#'   `zsim`, `inputs` (the `without` / `with` data tables), `guard` (per-group
#'   rank and condition numbers), and `results` (`NULL` while gated, otherwise the
#'   algorithm output for each prepared version).
#' @export
run_AR <- function(data,
                   y = "Y_observed",
                   d = "D_observed",
                   z = "assignment",
                   algorithm = c("algo2", "algo1"),
                   n_rand = 1000,
                   with_covariate = TRUE,
                   alpha = 0.95,
                   tol = 1e-8,
                   cond_threshold = 1e10,
                   seed = NULL,
                   execute = FALSE) {

    algorithm = match.arg(algorithm)

    # --- Validate the column mapping -----------------------------------------
    if (!is.data.frame(data)) {
        stop("`data` must be a data frame.")
    }
    key_cols = c(y, d, z)
    if (!all(vapply(key_cols, is.character, logical(1))) ||
        any(lengths(list(y, d, z)) != 1)) {
        stop("`y`, `d`, `z` must each be a single column name.")
    }
    if (length(unique(key_cols)) != 3) {
        stop("`y`, `d`, `z` must name three distinct columns.")
    }
    missing_cols = setdiff(key_cols, names(data))
    if (length(missing_cols) > 0) {
        stop("Column(s) not found in `data`: ",
             paste(missing_cols, collapse = ", "), ".")
    }

    # Covariates are everything that is not y/d/z.
    cov_names = setdiff(names(data), key_cols)

    # --- Standardize to the names the algorithms expect ----------------------
    std = data.frame(
        Y_observed = data[[y]],
        D_observed = data[[d]],
        assignment = data[[z]]
    )
    if (length(cov_names) > 0) {
        std = cbind(std, data[cov_names])
    }

    # --- Validate assignment / treatment columns -----------------------------
    if (!all(std$assignment %in% c(0, 1))) {
        stop("Assignment column `", z, "` must be binary (0/1).")
    }
    if (!all(std$D_observed %in% c(0, 1))) {
        stop("Treatment-received column `", d, "` must be binary (0/1).")
    }

    N  = nrow(std)
    N1 = sum(std$assignment == 1)
    N0 = sum(std$assignment == 0)
    if (N1 == 0 || N0 == 0) {
        stop("Both treatment arms must be non-empty; got N1 = ", N1,
             ", N0 = ", N0, ".")
    }

    # --- Decide which versions to prepare ------------------------------------
    has_cov = length(cov_names) > 0
    run_with = has_cov && isTRUE(with_covariate)
    # The without-covariates analysis is always prepared.
    runs = if (run_with) c("without", "with") else "without"

    # --- Report what the program understood ----------------------------------
    compliance_rate = mean(std$D_observed)
    p_d_treated = mean(std$D_observed[std$assignment == 1])
    p_d_control = mean(std$D_observed[std$assignment == 0])

    cat("========================================\n")
    cat("run_AR(): setup\n")
    cat("========================================\n")
    cat("  Algorithm            :", algorithm, "\n")
    cat("  Sample size N        :", N, "\n")
    cat("  Treated   N1 (Z = 1) :", N1, "\n")
    cat("  Control   N0 (Z = 0) :", N0, "\n")
    cat("  Randomizations       :", n_rand, "\n")
    cat("  Compliance rate      :", round(compliance_rate, 4),
        "(overall D = 1 rate)\n")
    cat("     P(D = 1 | Z = 1)  :", round(p_d_treated, 4), "\n")
    cat("     P(D = 1 | Z = 0)  :", round(p_d_control, 4), "\n")
    cat("     first stage       :", round(p_d_treated - p_d_control, 4),
        "(P(D=1|Z=1) - P(D=1|Z=0))\n")
    if (run_with) {
        cat("  Will return          : WITHOUT- and WITH-covariates results\n")
    } else if (has_cov) {
        cat("  Will return          : WITHOUT-covariates results only",
            "(with_covariate = FALSE)\n")
    } else {
        cat("  Will return          : WITHOUT-covariates results only",
            "(no covariates supplied)\n")
    }
    cat("\n")

    # --- Demean covariates ----------------------------------------------------
    if (has_cov) {
        non_numeric = cov_names[!vapply(std[cov_names], is.numeric,
                                        logical(1))]
        if (length(non_numeric) > 0) {
            stop("Covariate(s) must be numeric to be demeaned/conditioned: ",
                 paste(non_numeric, collapse = ", "), ".")
        }
        for (nm in cov_names) {
            std[[nm]] = std[[nm]] - mean(std[[nm]])
        }
        cat("Demeaned", length(cov_names), "covariate(s): ",
            paste(cov_names, collapse = ", "), "\n\n", sep = "")
    }

    # --- Invertibility guard (only matters for the with-covariates design) ---
    guard = NULL
    if (run_with) {
        guard = check_group_conditioning(std, cov_names, cond_threshold)
    }

    # --- Assemble the data tables for each version ---------------------------
    core_cols = c("Y_observed", "D_observed", "assignment")
    data_without = std[, core_cols, drop = FALSE]
    data_with = if (run_with) std else NULL

    # --- Generate the permuted-assignment matrix (zsim) ----------------------
    if (!is.null(seed)) {
        set.seed(seed)
    }
    zsim = sapply(seq_len(n_rand), gen_assignment_CR_index, N1 = N1, N0 = N0)

    # --- Dispatch / GATE ------------------------------------------------------
    algo_fun_name = switch(algorithm,
                           algo1 = "AR_algo1_custom",
                           algo2 = "AR_algo2_custom")

    results = list(without = NULL, with = NULL)

    if (!isTRUE(execute)) {
        # Gate closed: prepared everything, but do not touch the algorithms.
        for (r in runs) {
            cat("[GATE] prepared ", algo_fun_name, "() for the ",
                toupper(r), "-covariates data; not executed ",
                "(execute = FALSE).\n", sep = "")
        }
        cat("\n")
    } else {
        if (!exists(algo_fun_name, mode = "function")) {
            stop(algo_fun_name, "() not found. Source the algorithm file ",
                 "(e.g. B_AR_", algorithm, ".R) before calling with ",
                 "execute = TRUE.")
        }
        algo_fun = match.fun(algo_fun_name)
        for (r in runs) {
            dt = if (r == "with") data_with else data_without
            results[[r]] = algo_fun(dt, N1 = N1, N0 = N0, zsim = zsim,
                                    tol = tol, alpha = alpha)
        }
    }

    return(list(
        algorithm  = algorithm,
        alpha      = alpha,
        tol        = tol,
        N          = N,
        N1         = N1,
        N0         = N0,
        runs       = runs,
        covariates = cov_names,
        zsim       = zsim,
        inputs     = list(without = data_without, with = data_with),
        guard      = guard,
        results    = results[runs]
    ))
}

#' Check the per-treatment-group covariate design for ill-conditioning
#'
#' `solve_coef_01()` fits a separate regression within each treatment group, so
#' the design `[1, X]` must be full-rank and well-conditioned in BOTH the treated
#' and control groups. This helper computes each group's design rank and
#' condition number and stops if either group's design is rank-deficient or
#' exceeds `cond_threshold`.
#'
#' @param std Standardized data frame with `assignment` and the (already
#'   demeaned) covariate columns.
#' @param cov_names Character vector of covariate column names.
#' @param cond_threshold Condition-number cutoff.
#' @return A named list (`treated`, `control`) of `list(rank, ncol, cond)`
#'   diagnostics; invisibly, since it is called for its side effect of stopping
#'   on ill-conditioning.
#' @noRd
check_group_conditioning <- function(std, cov_names, cond_threshold) {

    groups = list(treated = std$assignment == 1,
                  control = std$assignment == 0)

    diagnostics = list()
    for (g in names(groups)) {
        idx = groups[[g]]
        X = as.matrix(std[idx, cov_names, drop = FALSE])
        design = cbind(`(Intercept)` = 1, X)

        rank = qr(design)$rank
        full = ncol(design)
        cond = if (rank < full) Inf else kappa(design, exact = FALSE)

        diagnostics[[g]] = list(rank = rank, ncol = full, cond = cond)

        if (rank < full) {
            stop("Covariate design [1, X] is rank-deficient in the ", g,
                 " group (rank ", rank, " < ", full, " columns): covariates ",
                 "are collinear. Covariates: ",
                 paste(cov_names, collapse = ", "), ".")
        }
        if (cond > cond_threshold) {
            stop("Covariate design [1, X] is ill-conditioned in the ", g,
                 " group (condition number ", format(cond, digits = 3),
                 " > ", format(cond_threshold, digits = 3),
                 "). Near-collinear covariates: ",
                 paste(cov_names, collapse = ", "), ".")
        }
    }

    invisible(diagnostics)
}
