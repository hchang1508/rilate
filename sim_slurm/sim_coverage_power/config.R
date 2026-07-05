################################################################################
## config.R
## Shared configuration + helper functions for the coverage/power simulation.
##
## The design mirrors C_sim_gam_tau.R's gamma x tau grid, translated into the
## rilate package's principal-stratum DGP (assess_coverage / gen_sim_data):
##
##   * gamma (instrument strength)  ->  complier fraction fc
##       gamma = 0   -> fc = 0    -> fractions = c(0.5,  0,   0.5)   (weak IV)
##       gamma = 0.5 -> fc = 0.5  -> fractions = c(0.25, 0.5, 0.25)  (strong IV)
##       gamma = 1   -> fc = 1    -> fractions = c(0,    1,   0)     (strong IV)
##     (non-complier mass split evenly between always-/never-takers)
##   * tau (treatment effect / LATE) -> tau_c = taus[2], taus = c(0, tau, 0)
##
## Every parameter cell is run with BOTH Algorithm 1 and Algorithm 2 on the
## *same* observed data and the *same* permutation matrix, so their confidence
## sets are directly comparable via Hausdorff distance.
################################################################################

## ---- Fixed design constants -------------------------------------------------
SIM_CFG <- list(
    N            = 100,      # total units
    N1           = 50,       # units assigned to treatment (Z = 1)
    nsim         = 10000,    # simulations per parameter cell
    chunks       = 1000,     # array chunks per cell (10 sims each -> 9000 tasks)
    n_rand       = 1000,     # permuted assignments per run_rilate() call
    n_rand_hi    = 10000,    # high-resolution permutation count (algo2 only)
    run_hi       = FALSE,    # run the algo2 @ n_rand_hi case? (OFF for now -- it
                             # costs ~1 h/sim on a cluster core; enable per-run
                             # with env RUN_HI=1 once the 1000 sweep is validated)
    n_cov        = 2,        # # pure-noise N(0,1) covariates (no predictive power)
    alpha        = 0.95,     # confidence level (sig level = 1 - alpha = 0.05)
    mode         = "constant",
    flag_tol     = 1e-4,     # Hausdorff distance above which algo1/algo2 differ
    # RNG namespaces (kept far apart so streams never collide)
    seed_true    = 10000L,   # per-cell fixed "true table" seed base
    seed_obs     = 200000L,  # per-simulation observed-draw seed base
    seed_cov     = 400000L,  # per-simulation covariate-draw seed base
    seed_zsim    = 900000L   # per-simulation permutation-matrix seed base
)
SIM_CFG$nsim_per_chunk <- SIM_CFG$nsim / SIM_CFG$chunks   # 10

## ---- The 9-cell gamma x tau grid --------------------------------------------
## Ordered so combo_id = (gamma_idx - 1) * 3 + tau_idx, gamma outer / tau inner.
FC_LEVELS  <- c(0, 0.5, 1)     # gamma -> complier fraction
TAU_LEVELS <- c(0, 0.5, 1)     # tau_c (true LATE)

## Even split of the non-complier mass into (always, never)-takers.
fractions_for_fc <- function(fc) {
    rest <- (1 - fc) / 2
    c(fa = rest, fc = fc, fn = rest)
}

build_grid <- function() {
    rows <- list()
    k <- 0L
    for (gi in seq_along(FC_LEVELS)) {
        for (ti in seq_along(TAU_LEVELS)) {
            k <- k + 1L
            fc  <- FC_LEVELS[gi]
            tau <- TAU_LEVELS[ti]
            fr  <- fractions_for_fc(fc)
            rows[[k]] <- data.frame(
                combo_id = k,
                gamma    = fc,          # instrument strength (== fc here)
                fc       = fc,
                fa       = fr["fa"],
                fn       = fr["fn"],
                tau      = tau,         # true LATE (tau_c)
                stringsAsFactors = FALSE
            )
        }
    }
    grid <- do.call(rbind, rows)
    rownames(grid) <- NULL
    grid
}

SIM_GRID <- build_grid()

################################################################################
## Confidence-set helpers (self-contained; do not rely on package internals)
################################################################################

## Coerce a confidence set to a clean list of length-2 numeric intervals.
.as_interval_list <- function(cs) {
    if (is.null(cs) || length(cs) == 0) return(list())
    if (is.numeric(cs) && length(cs) == 2) return(list(as.numeric(cs)))
    lapply(cs, function(iv) as.numeric(iv[1:2]))
}

## Is `point` covered by the confidence set (union of intervals)?
point_in_cs <- function(point, cs, tol = 1e-8) {
    ivs <- .as_interval_list(cs)
    if (length(ivs) == 0) return(FALSE)
    for (iv in ivs) {
        if (point >= iv[1] - tol && point <= iv[2] + tol) return(TRUE)
    }
    FALSE
}

## Does the confidence set contain an unbounded interval?
cs_is_infinite <- function(cs) {
    ivs <- .as_interval_list(cs)
    if (length(ivs) == 0) return(FALSE)
    any(vapply(ivs, function(iv) !is.finite(iv[1]) || !is.finite(iv[2]),
               logical(1)))
}

## Total (finite) width of the confidence set; Inf if any interval unbounded.
cs_width <- function(cs) {
    ivs <- .as_interval_list(cs)
    if (length(ivs) == 0) return(0)
    sum(vapply(ivs, function(iv) iv[2] - iv[1], numeric(1)))
}

## Sort + merge overlapping/adjacent intervals into a canonical list.
.merge_intervals <- function(ivs, tol = 1e-12) {
    if (length(ivs) == 0) return(list())
    m <- do.call(rbind, ivs)
    m <- m[order(m[, 1]), , drop = FALSE]
    out <- list()
    cur <- m[1, ]
    if (nrow(m) >= 2) {
        for (i in 2:nrow(m)) {
            if (m[i, 1] <= cur[2] + tol) {
                cur[2] <- max(cur[2], m[i, 2])
            } else {
                out[[length(out) + 1]] <- cur
                cur <- m[i, ]
            }
        }
    }
    out[[length(out) + 1]] <- cur
    out
}

## Are two interval lists identical (matching finite endpoints / infinities)?
.sets_equal <- function(A, B, tol = 1e-9) {
    if (length(A) != length(B)) return(FALSE)
    for (i in seq_along(A)) {
        for (j in 1:2) {
            a <- A[[i]][j]; b <- B[[i]][j]
            if (is.finite(a) && is.finite(b)) {
                if (abs(a - b) > tol) return(FALSE)
            } else if (a != b) {         # both infinite: signs must match
                return(FALSE)
            }
        }
    }
    TRUE
}

## Distance from a finite point x to the union of (finite) intervals.
.point_to_set <- function(x, ivs) {
    min(vapply(ivs, function(iv) {
        if (x < iv[1]) iv[1] - x else if (x > iv[2]) x - iv[2] else 0
    }, numeric(1)))
}

## Directed Hausdorff distance sup_{a in A} d(a, B), A/B finite interval lists.
## d(.,B) is 1-Lipschitz, piecewise linear; its max over A is attained at an
## endpoint of A or at the midpoint of a gap between two B intervals that lies
## inside A. Evaluating those candidate points is exact.
.directed_hausdorff <- function(A, B) {
    cand <- unlist(lapply(A, function(iv) c(iv[1], iv[2])))
    if (length(B) >= 2) {
        Bs <- B[order(vapply(B, `[`, numeric(1), 1))]
        for (i in seq_len(length(Bs) - 1)) {
            mid <- (Bs[[i]][2] + Bs[[i + 1]][1]) / 2
            inside_A <- any(vapply(A, function(iv) mid >= iv[1] & mid <= iv[2],
                                   logical(1)))
            if (inside_A) cand <- c(cand, mid)
        }
    }
    max(vapply(cand, .point_to_set, numeric(1), ivs = B))
}

## Hausdorff distance between two confidence sets (unions of intervals).
## Infinite sets: 0 if identical, otherwise Inf (they cannot be within any
## finite tolerance of each other).
hausdorff_cs <- function(csA, csB) {
    A <- .merge_intervals(.as_interval_list(csA))
    B <- .merge_intervals(.as_interval_list(csB))
    if (length(A) == 0 && length(B) == 0) return(0)
    if (length(A) == 0 || length(B) == 0) return(Inf)
    a_inf <- any(!is.finite(unlist(A)))
    b_inf <- any(!is.finite(unlist(B)))
    if (a_inf || b_inf) {
        return(if (.sets_equal(A, B)) 0 else Inf)
    }
    max(.directed_hausdorff(A, B), .directed_hausdorff(B, A))
}
