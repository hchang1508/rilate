RI for LATE (under construction)
================

**Author:** Arya Gadage and Haoge Chang<br/> **Last updated:** July 04,
2026

# Introduction

This package implements the randomization-based inferential procedure
for the Local Average Treatment Effect (LATE) developed in
[\[2\]](https://academic.oup.com/biomet/article/113/2/asag010/8487895)
and [\[4\]](https://academic.oup.com/jrsssa/article/168/1/109/7084141).

# Overview

In randomized experiments, participants do not always take the treatment
they were assigned to — the treatment they actually receive can differ
from their assignment. This gap is called **noncompliance**. When it is
present, a natural target of inference is the **Local Average Treatment
Effect (LATE)** [\[3\]](https://www.jstor.org/stable/2951620): the
average causal effect of treatment among **compliers**, the
subpopulation whose treatment status changes as a result of treatment
assignment.

To state this precisely, consider a sample of individuals
$i = 1, \ldots, n$, and for each let

- $Y_i(1), Y_i(0)$ — the potential outcomes under treatment and control;
- $D_i(1), D_i(0)$ — the treatment actually taken when assigned to
  treatment and to control.

A complier is an individual with $D_i(1) > D_i(0)$: they take the
treatment when assigned to it, but not otherwise. We use
$\mathcal{C} = \lbrace i : D_i(1) > D_i(0) \rbrace$ to denote the set of
compliers. The LATE is the average treatment effect over this group,

$$
\tau_{\text{LATE}} = \frac{1}{n_c} \sum_{i \in \mathcal{C}} \left( Y_i(1) - Y_i(0) \right)
$$

where $n_c$ is the number of compliers. This package provides a
randomization-based inferential procedure for $\tau_{\text{LATE}}$ with
two complementary guarantees:

- **Tests and p-values** are finite-sample valid for testing the sharp
  null $Y_i(1) = Y_i(0)$ for all $i \in \mathcal{C}$, and asymptotically
  valid for testing $\tau_{\text{LATE}} =0$ (under treatment-effect
  heterogeneity).
- **Confidence sets** are finite-sample valid under treatment-effect
  additive homogeneity ($Y_i(1) = Y_i(0) + \tau$ for all
  $i \in \mathcal{C}$, with $\tau$ not depending on $i$), and
  asymptotically valid for the LATE under heterogeneity. The procedures
  are robust against low compliance rates.

The package supports covariate adjustment via linear regression, which
can improve precision and guard against chance covariate imbalance.

## Complete Randomization Design (Important)

The package currently supports completely randomized designs with a
binary treatment. The remaining statistical regularity conditions for
asymptotic guarantees are stated in
[\[2\]](https://academic.oup.com/biomet/article/113/2/asag010/8487895).

# Note

1.  The package includes two algorithms that implement the same
    inferential procedure but differ in speed. Algorithm 2 is faster and
    is used by default; Algorithm 1 is retained for comparison and
    validation.
2.  The authors used Claude to migrate the original code from
    [\[2\]](https://academic.oup.com/biomet/article/113/2/asag010/8487895)
    into this repository and to help design the test cases.

# Demonstration

## Install packages

`rilate` is not yet on CRAN. Install the development version from
[GitHub](https://github.com/hchang1508/rilate) with:

``` r
# install.packages("remote")
remotes::install_github("hchang1508/rilate")
```

Then load it with `library(rilate)`, as in the example below.

## run_rilate function

The main entry point is `run_rilate()`. The `y`, `d`, and `z` arguments
name the outcome, treatment-received, and assignment columns of `data`
(and default to `"Y_observed"`, `"D_observed"`, `"assignment"`). Passing
`x` runs *both* an unadjusted analysis and a covariate-adjusted one.
Each returns a confidence set for the LATE and a randomization p-value
for the sharp null and weak null hypotheses. `run_rilate()` also prints
a human-readable setup report to the console. The number of simulated
randomization draw samples is controlled by the `n_rand` argument and
defaults to 1000.

``` r
library(rilate)

#####################################################################################
########generate random data ########################################################
#####################################################################################
set.seed(7)
sim <- gen_sim_data(
  N         = 800,
  N1        = 400,                  # 400 assigned to treatment, 400 to control
  fractions = c(0.20, 0.50, 0.30),  # always-takers, compliers, never-takers
  taus      = c(0, 1, 0),           # principal-stratum effects; complier LATE = 1
  mode      = "constant"
)
dat <- as.data.frame(sim$observed)  # columns: Y_observed, D_observed, assignment


#####################################################################################
########run run_rilate ##############################################################
#####################################################################################
# run_rilate() prints a setup report and a per-step progress trace to the
# console; pass verbose = FALSE to run silently and keep the returned object.
res_nocov <- run_rilate(dat, y = "Y_observed", d = "D_observed",
                        z = "assignment", n_rand = 1000, verbose = FALSE)

res_nocov$results$without_covariates$confidence_set   # 95% confidence set for the LATE
#> [[1]]
#> [1] 0.8491505 1.5757713
res_nocov$results$without_covariates$p_value          # randomization p-value (H0: LATE = 0)
#> [1] 0.000999001
```

`run_rilate()` returns a list with the prepared inputs and metadata
(`algorithm`, `alpha`, `tol`, `N`, `N1`, `N0`, `runs`, `covariates`,
`zsim`, `inputs`, `guard`) plus `results`. Each element of `results`
(`without_covariates` and, when applicable, `with_covariates`) holds a
`confidence_set` (a list of `[lower, upper]` intervals) and a `p_value`
(the randomization p-value for the null $\tau_{\text{LATE}} = 0$).

## run_rilate() with covariates

We can use pretreatment covariates to improve the precision of our
analysis.

``` r
set.seed(7)
n <- nrow(dat)
dat$x1 <- rnorm(n)
dat$x2 <- rnorm(n)
dat$Y_observed <- dat$Y_observed + 2 * dat$x1 - 1.5 * dat$x2

head(dat)
#>     Y_observed D_observed assignment         x1         x2
#> 1  6.668696762          1          1  2.2872472  0.1286965
#> 2 -2.765516317          1          0 -1.1967717 -0.5498658
#> 3 -0.003811554          1          0 -0.6942925 -1.3860440
#> 4 -1.024435857          1          0 -0.4122930 -0.1416287
#> 5 -1.468333867          1          1 -0.9706733 -0.9624574
#> 6 -4.048748998          1          1 -0.9472799  0.8046061

# Passing `x` runs both an unadjusted and a covariate-adjusted analysis.
res <- run_rilate(dat, y = "Y_observed", d = "D_observed", z = "assignment",
                  x = c("x1", "x2"), n_rand = 1000, seed = 1, verbose = FALSE)
```

The result carries a `without_covariates` and a `with_covariates` entry,
each with a `confidence_set` and a `p_value`:

``` r
fmt_cs   <- function(cs) paste(vapply(cs, function(iv)
  sprintf("[%.3f, %.3f]", iv[1], iv[2]), character(1)), collapse = " U ")
cs_width <- function(cs) sum(vapply(cs, function(iv) diff(iv), numeric(1)))

data.frame(
  analysis = c("without covariates", "with covariates"),
  p_value  = round(c(res$results$without_covariates$p_value, res$results$with_covariates$p_value), 4),
  CS_95    = c(fmt_cs(res$results$without_covariates$confidence_set),
               fmt_cs(res$results$with_covariates$confidence_set)),
  width    = round(c(cs_width(res$results$without_covariates$confidence_set),
                     cs_width(res$results$with_covariates$confidence_set)), 3),
  row.names = NULL
)
#>             analysis p_value          CS_95 width
#> 1 without covariates   0.022 [0.198, 2.422] 2.224
#> 2    with covariates   0.001 [1.070, 1.382] 0.312
```

Adjusting for the two prognostic covariates shrinks the confidence set
roughly sevenfold (from about `[0.20, 2.42]`, width `2.22`, to
`[1.07, 1.38]`, width `0.31`) and sharpens the p-value for
$\tau_{\text{LATE}} = 0$ from about `0.022` to `0.001`. This is the
precision gain covariate adjustment buys when the covariates are
predictive for the outcome.

## run_rilate() parameter references

Below are its inputs, split into those you must supply and those that
have sensible defaults.

### Required inputs

| Argument | Type | Description |
|----|----|----|
| `data` | data frame | Contains the outcome, treatment-received, assignment, and (optionally) covariate columns. |

The `y`, `d`, and `z` arguments below name columns of `data`. They are
technically optional because they default to `"Y_observed"`,
`"D_observed"`, and `"assignment"`, but the analysis is undefined
without valid outcome, treatment, and assignment columns — so treat them
as required unless your columns already use those default names.

| Argument | Default | Description |
|----|----|----|
| `y` | `"Y_observed"` | Name of the outcome column in `data`. |
| `d` | `"D_observed"` | Name of the treatment-*received* column in `data`. Must be binary (0/1). |
| `z` | `"assignment"` | Name of the treatment-*assignment* column in `data`. Must be binary (0/1). |

### Optional inputs

| Argument | Default | Description |
|----|----|----|
| `x` | `NULL` | Character vector naming covariate columns to adjust for. When `NULL`, no covariates are used and the covariate-adjusted analysis is skipped. Named columns must exist in `data`, must be numeric, and must not overlap `y`/`d`/`z`. |
| `algorithm` | `"algo2"` | Which algorithm to run: `"algo2"` (faster, default) or `"algo1"` (retained for validation). Both implement the same procedure. |
| `n_rand` | `1000` | Number of randomizations (permuted assignments) to draw. Larger values give more precise p-values at higher computational cost. |
| `with_covariate` | `TRUE` | When `TRUE` and covariates are supplied, prepare *both* an unadjusted and a covariate-adjusted analysis. When `FALSE`, prepare only the unadjusted analysis. |
| `alpha` | `0.95` | Confidence level for the returned confidence set. |
| `tol` | `1e-8` | Numerical tolerance passed to the algorithm. |
| `cond_threshold` | `1e10` | Condition-number cutoff above which a treatment group’s covariate design `[1, X]` is treated as ill-conditioned (triggers an error). Only relevant when covariates are used. |
| `seed` | `NULL` | Optional integer seed for reproducible generation of the permuted-assignment matrix (`zsim`). |

# Assess Coverage

# Compare Two algorithms (time )

# Miscellaneous

## gen_sim_data() references

`gen_sim_data()` simulates a completely randomized experiment with
one-sided noncompliance in a single call. It first builds a
potential-outcome (“true”) table for a population split into three
principal strata — always-takers, compliers, and never-takers — and then
draws one completely randomized assignment from it, revealing the
observed treatment received `D = d1 * Z + d0 * (1 - Z)` and outcome
`Y = y1 * D + y0 * (1 - D)`. It is the data generator used throughout
these examples and in the package’s simulation tests.

The stratum counts are derived from `fractions`: `N_at = round(fa * N)`,
`N_c = round(fc * N)`, and the never-taker count `N_nt = N - N_at - N_c`
absorbs any rounding remainder so the three counts always sum to `N`
exactly. By default (`mode = "constant"`) every unit in a stratum
receives exactly that stratum’s treatment effect;
`mode = "heterogeneous"` instead draws unit-level effects
`tau * rnorm(n, 1, 0.1)` and re-centers each stratum’s mean back to its
target `tau`.

### Inputs

| Argument | Type / Default | Description |
|----|----|----|
| `N` | integer | Total number of units. |
| `N1` | integer | Number of units assigned to treatment (`Z = 1`); the remaining `N - N1` units are controls. Must satisfy `0 <= N1 <= N`. |
| `fractions` | `c(fa, fc, fn)` | Population fractions for always-takers, compliers, and never-takers. Must sum to 1 (up to floating-point error). |
| `taus` | `c(tau_a, tau_c, tau_n)` | Principal-stratum treatment effects for always-takers, compliers, and never-takers. `tau_c` is the complier LATE. |
| `mode` | `"constant"` | `"constant"` gives every unit in a stratum exactly its `tau`; `"heterogeneous"` draws unit-level effects around each stratum’s target and calibrates the stratum mean to `tau`. |

### Value

A list with two elements:

- `true_table` — the full potential-outcome table (columns `indices`,
  `y1`, `y0`, `d1`, `d0`), with rows ordered always-takers,
  never-takers, then compliers.
- `observed` — one realized sample (columns `Y_observed`, `D_observed`,
  `assignment`), one row per unit. This is the data frame you pass to
  `run_rilate()`.

``` r
set.seed(7)
sim <- gen_sim_data(
  N         = 800,
  N1        = 400,                  # 400 assigned to treatment, 400 to control
  fractions = c(0.20, 0.50, 0.30),  # always-takers, compliers, never-takers
  taus      = c(0, 1, 0),           # complier LATE = 1; no effect for others
  mode      = "constant"
)
head(as.data.frame(sim$observed))
#>   Y_observed D_observed assignment
#> 1  2.2872472          1          1
#> 2 -1.1967717          1          0
#> 3 -0.6942925          1          0
#> 4 -0.4122930          1          0
#> 5 -0.9706733          1          1
#> 6 -0.9472799          1          1
```

### How the baseline outcome `y0` is generated

The baseline (control) potential outcome $Y_i(0)$ is drawn once, up
front, and is never modified afterward. Each unit gets an i.i.d.
standard-normal draw,

$$
Y_i(0) \sim \mathcal{N}(0, 1),
$$

independently of the unit’s compliance type (always-taker, complier, or
never-taker). At this “raw” stage the treated outcome is initialized to
the same value, $Y_i(1) = Y_i(0)$, so there is no treatment effect yet.

The treatment effect is then added **only to the `y1` column**; `y0` is
left untouched:

- **`mode = "constant"`** — every unit in a stratum gets exactly its
  stratum effect, so $Y_i(1) = Y_i(0) + \tau_{s(i)}$ where $s(i)$ is the
  unit’s stratum.
- **`mode = "heterogeneous"`** — each unit draws
  $Y_i(1) = Y_i(0) + \tau_{s(i)}
  \cdot \varepsilon_i$ with $\varepsilon_i \sim \mathcal{N}(1, 0.1)$,
  and the stratum is then re-centered so that
  $\operatorname{mean}(y1 - y0) = \tau_{s(i)}$ exactly within each
  stratum.

Two consequences worth noting:

1.  The treatment effect enters **entirely through `y1`**. The baseline
    `y0` is homogeneous noise that carries no information about
    compliance type — it is identically distributed across all three
    strata.
2.  Because `y0` has no covariate structure of its own, the covariate
    example above must *add* prognostic covariates on top of
    `Y_observed` for covariate adjustment to have anything to explain.
    If you need stratum-specific baselines, that would have to be added
    to `gen_outcome_raw_normal()`; it is not modeled by default.

# References

\[1\] Anderson, T. W., & Rubin, H. (1949). Estimation of the parameters
of a single equation in a complete system of stochastic equations.
*Annals of Mathematical Statistics*, 20(1), 46-63.

\[2\] Aronow, P. M., Chang, H., & Lopatto, P. (2026).
Randomization-Based Confidence Sets for the Local Average Treatment
Effect. *Biometrika*, 113(2).

\[3\] Imbens, G. W., & Angrist, J. D. (1994). Identification and
Estimation of Local Average Treatment Effects. Econometrica, 62(2),
467-475.

\[4\] Imbens, G. W., & Rosenbaum, P. R. (2005). Robust, accurate
confidence intervals with a weak instrument: quarter of birth and
education. *Journal of the Royal Statistical Society: Series A*, 168(1),
109-126.
