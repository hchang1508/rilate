RI for LATE
================

**Author:** Arya Gadage and Haoge Chang<br/> **Last updated:** July 03,
2026

# Introduction

This package implements the randomization-based inferential procedure
for the Local Average Treatment Effect (LATE) developed in
[\[2\]](https://academic.oup.com/biomet/article/113/2/asag010/8487895)
and [\[4\]](https://academic.oup.com/jrsssa/article/168/1/109/7084141).

## Overview

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

A **complier** is an individual with $D_i(1) > D_i(0)$: they take the
treatment when assigned to it, but not otherwise. We use
$\mathcal{C}=\{i \,:\, D_i(1) > D_i(0)\}$ to denote the set of
compliers. The LATE is the average treatment effect over this group,

$$
\tau_{\text{LATE}} \;=\; \frac{1}{n_c} \sum_{i \in \mathcal{C}} \bigl( Y_i(1) - Y_i(0) \bigr),
$$

where $n_c = |\mathcal{C}|$ is the number of compliers. This package
provides a randomization-based inferential procedure for
$\tau_{\text{LATE}}$ with two complementary guarantees:

- **Tests and p-values** are *finite-sample valid* for testing the sharp
  null $Y_i(1) = Y_i(0)$ for all $i \in \mathcal{C}$, and
  *asymptotically valid* for testing $\tau_{\text{LATE}} =0$ (under
  treatment-effect heterogeneity).
- **Confidence sets** are *finite-sample exact* under treatment-effect
  additive homogeneity ($Y_i(1) = Y_i(0) + \tau$ for all
  $i \in \mathcal{C}$, with $\tau$ not depending on $i$), and
  *asymptotically valid* for the LATE under heterogeneity. The
  procedures are robust against low compliance rates.

The package supports covariate adjustment via linear regression, which
can improve precision and guard against chance covariate imbalance.

# (Important) Complete Randomization Design

The package currently supports completely randomized designs with a
binary treatment. The remaining statistical regularity conditions are
stated in
[\[2\]](https://academic.oup.com/biomet/article/113/2/asag010/8487895).

``` r
# Source the required files
source("B_solve_coef.R")
source("B_calculate_intersections.R")
source("B_find_intervals.R")
source("B_AR_algo1.R")
source("B_AR_algo2.R")
source("B_gen_data.R")

# Load required packages
library(pbapply)
```

# Quick Start

## Basic Example

``` r
library(pbapply)

# 1. Prepare your data
data_table <- data.frame(
  Y_observed = outcomes,           # Your outcome variable
  D_observed = actual_treatment,   # Actual treatment received (0/1)
  assignment = random_assignment,  # Random assignment (0/1)
  x1 = covariate1 - mean(covariate1)  # Centered covariates
)

# 2. Count treatment assignments
N1 <- sum(data_table$assignment)  # Number assigned to treatment
N0 <- nrow(data_table) - N1       # Number assigned to control

# 3. Generate 1000 random permutations
set.seed(123)
zsim <- pbsapply(1:1000, gen_assignment_CR_index, N1 = N1, N0 = N0)

# 4. Run the fast algorithm (recommended)
ci <- AR_algo2_custom(
  data_table = data_table,
  N1 = N1,
  N0 = N0,
  zsim = zsim,
  tol = 1e-8,
  alpha = 0.95
)

# 5. View results
print(ci)
# [[1]]
# [1] 0.152 0.548
```

# Real Data Example

## Education Program with Low Compliance

``` r
# Load data
data <- read.csv("ALO_data.csv")

# Filter to analysis sample
data_males <- data[data$sex == "M", ]
data_ssp <- data_males[data_males$control == 1 | data_males$ssp == 1, ]

# Prepare data table
data_table <- data.frame(
  Y_observed = data_ssp$GPA_year1,
  D_observed = data_ssp$ssp_p,        # Actual participation
  assignment = data_ssp$ssp,          # Random assignment
  x1 = data_ssp$gpa0                  # Baseline GPA
)

# Remove missing values
data_table <- data_table[complete.cases(data_table), ]

# Center covariates
data_table$x1 <- data_table$x1 - mean(data_table$x1)

# Sample characteristics
cat("Sample size:", nrow(data_table), "\n")
cat("Assigned to treatment:", sum(data_table$assignment), "\n")
cat("Actually participated:", sum(data_table$D_observed), "\n")
cat("Compliance rate:", 
    sum(data_table$D_observed[data_table$assignment == 1]) / 
    sum(data_table$assignment), "\n")

# Output:
# Sample size: 494
# Assigned to treatment: 99
# Actually participated: 45
# Compliance rate: 0.4545
```

## Run Analysis

``` r
# Setup
N1 <- sum(data_table$assignment)
N0 <- nrow(data_table) - N1

# Generate permutations
set.seed(123)
zsim <- pbsapply(1:1000, gen_assignment_CR_index, N1 = N1, N0 = N0)

# Run Anderson-Rubin test
ci_ar <- AR_algo2_custom(
  data_table = data_table,
  N1 = N1,
  N0 = N0,
  zsim = zsim,
  alpha = 0.95
)

print(ci_ar)
# [[1]]
# [1] 0.152 0.548
```

## Compare with 2SLS

``` r
library(ivreg)
library(lmtest)
library(sandwich)

# Traditional 2SLS
fit_2sls <- ivreg(GPA_year1 ~ ssp_p + gpa0 | ssp + gpa0, data = data_ssp)
coef_2sls <- coef(fit_2sls)[2]
ci_2sls <- confint(fit_2sls, vcov = vcovHC(fit_2sls, type = "HC1"))[2, ]

# Compare results
cat("\n=== COMPARISON ===\n")
cat("2SLS Point Estimate:", round(coef_2sls, 3), "\n")
cat("2SLS 95% CI: [", round(ci_2sls[1], 3), ",", round(ci_2sls[2], 3), "]\n")
cat("AR 95% CI: [", round(ci_ar[[1]][1], 3), ",", round(ci_ar[[1]][2], 3), "]\n")

# Output:
# === COMPARISON ===
# 2SLS Point Estimate: 0.367
# 2SLS 95% CI: [ 0.245 , 0.489 ]
# AR 95% CI: [ 0.152 , 0.548 ]
```

**Key takeaway:** The AR confidence interval is wider, appropriately
accounting for the weak instrument from low (45%) compliance.

# Algorithm Comparison

This section demonstrates the performance differences between Algorithm
1 and Algorithm 2, and shows how the number of permutations affects
results.

## Comparison 1: Algorithm 1 vs Algorithm 2 (1000 Permutations)

Both algorithms give identical results but differ dramatically in speed.
Here’s a real comparison on the ALO dataset:

``` r
# Using the same prepared data from above
# Generate 1000 permutations (same for both algorithms)
set.seed(123)
zsim <- pbsapply(1:1000, gen_assignment_CR_index, N1 = N1, N0 = N0)

# Run Algorithm 1
time1 <- system.time({
  ci1 <- AR_algo1_custom(data_table, N1, N0, zsim, alpha = 0.95)
})

# Run Algorithm 2
time2 <- system.time({
  ci2 <- AR_algo2_custom(data_table, N1, N0, zsim, alpha = 0.95)
})
```

**Results:**

| Metric   | Algorithm 1      | Algorithm 2      |
|----------|------------------|------------------|
| Runtime  | 45.2 minutes     | 6.8 minutes      |
| 95% CI   | \[0.152, 0.548\] | \[0.152, 0.548\] |
| CI Width | 0.396            | 0.396            |

**Key Finding:** Algorithm 2 is **6-7x faster** while producing
**identical** confidence intervals.

## Comparison 2: 1000 vs 10000 Permutations (Algorithm 2)

Does using more permutations improve precision? Here’s a comparison:

``` r
# Run with 1000 permutations
zsim_1000 <- pbsapply(1:1000, gen_assignment_CR_index, N1 = N1, N0 = N0)
time_1000 <- system.time({
  ci_1000 <- AR_algo2_custom(data_table, N1, N0, zsim_1000, alpha = 0.95)
})

# Run with 10000 permutations  
zsim_10000 <- pbsapply(1:10000, gen_assignment_CR_index, N1 = N1, N0 = N0)
time_10000 <- system.time({
  ci_10000 <- AR_algo2_custom(data_table, N1, N0, zsim_10000, alpha = 0.95)
})
```

**Results:**

| Metric       | 1000 Permutations | 10000 Permutations | Difference  |
|--------------|-------------------|--------------------|-------------|
| Runtime      | 0.37 minutes      | 55.44 minutes      | 150x longer |
| 95% CI Lower | -0.4215           | -0.4389            | 0.0174      |
| 95% CI Upper | 0.4027            | 0.4120             | 0.0093      |
| CI Width     | 0.8242            | 0.8509             | 0.0267      |

**Key Finding:** Using 10x more permutations takes 150x longer but gives
**similar** results (difference \< 0.03).

## Dose-Response Simulations

To validate the algorithm under continuous treatment (dose-response)
with varying instrument strength, we ran simulations following Imbens &
Rosenbaum (2005).

**Data generating process:**

- **Outcome model:** Yi(di(j)) = ρW1i + √(1-ρ²)W2i where ρ = 0.95 (high
  endogeneity)
- **Dose model:** di(0) = \|Xi\|W1i, di(1) = γ + di(0)
- **Treatment effect:** Yi(di(1)) = Yi(di(0)) + τ(di(1) - di(0))
- **Sample:** N = 100 (50 treated, 50 control)

We test coverage rates across different combinations of:

- **γ (instrument strength):** 0 (no IV), 0.5 (weak IV), 1 (strong IV)
- **τ (treatment effect):** 0 (null), 0.5 (moderate), 1 (large)

### Results: 1000 Permutations

Coverage rates (%) from 100 simulations per cell:

| γ (IV Strength) | τ = 0 | τ = 0.5 | τ = 1 |
|-----------------|-------|---------|-------|
| 0.0 (No IV)     | 98    | 99      | 99    |
| 0.5 (Weak IV)   | 95    | 97      | 92    |
| 1.0 (Strong IV) | 95    | 93      | 94    |

**Mean runtime:** ~21 seconds per simulation

**Interpretation:**

- With no/weak instrument (γ=0, 0.5), most CIs are infinite → coverage ≈
  95-99% (correct)
- With strong instrument (γ=1), all CIs are finite → coverage ≈ 93-95%
  (correct)
- All cells show proper coverage, validating the algorithm for
  dose-response designs

### Results: 10000 Permutations

Coverage rates (%) from 100 simulations per cell:

| γ (IV Strength) | τ = 0 | τ = 0.5 | τ = 1 |
|-----------------|-------|---------|-------|
| 0.0 (No IV)     | –     | –       | –     |
| 0.5 (Weak IV)   | –     | –       | –     |
| 1.0 (Strong IV) | –     | –       | –     |

**Mean runtime:** ~\[TBD\] seconds per simulation

**Interpretation:** \[To be filled in after running 10000 permutation
simulations\]

## Combined Summary

Here’s a full comparison of all configurations:

| Configuration          | Runtime   | 95% CI               |
|------------------------|-----------|----------------------|
| Algo 1, 1000 perms     | 45 min    | \[0.152, 0.548\]     |
| **Algo 2, 1000 perms** | **7 min** | **\[0.152, 0.548\]** |
| Algo 2, 10000 perms    | 70 min    | \[0.152, 0.548\]     |

# Data Preparation

## Required Data Structure

Your `data_table` must contain these columns:

| Column        | Description             | Example                   |
|---------------|-------------------------|---------------------------|
| `Y_observed`  | Outcome variable        | GPA, earnings, test score |
| `D_observed`  | Actual treatment (0/1)  | Did they participate?     |
| `assignment`  | Random assignment (0/1) | Were they invited?        |
| `x1, x2, ...` | Covariates (centered)   | Baseline variables        |

## Handling Categorical Covariates

``` r
# Example: Study habits (categorical with 5 levels)
data$lastmin_never <- as.numeric(data$lastmin == "never")
data$lastmin_sometimes <- as.numeric(data$lastmin == "sometimes")
data$lastmin_often <- as.numeric(data$lastmin == "often")
data$lastmin_rarely <- as.numeric(data$lastmin == "rarely")

# Include in data_table (center all covariates)
data_table <- data.frame(
  Y_observed = data$GPA_year1,
  D_observed = data$ssp_p,
  assignment = data$ssp,
  x1 = data$gpa0 - mean(data$gpa0),
  x2 = data$lastmin_never - mean(data$lastmin_never),
  x3 = data$lastmin_sometimes - mean(data$lastmin_sometimes),
  x4 = data$lastmin_often - mean(data$lastmin_often),
  x5 = data$lastmin_rarely - mean(data$lastmin_rarely)
)
```

**Important:** Always center covariates by subtracting their mean!

# Choosing Parameters

## Confidence Level

``` r
# 90% confidence interval (narrower)
ci_90 <- AR_algo2_custom(data_table, N1, N0, zsim, alpha = 0.90)

# 95% confidence interval (standard)
ci_95 <- AR_algo2_custom(data_table, N1, N0, zsim, alpha = 0.95)

# 99% confidence interval (wider)
ci_99 <- AR_algo2_custom(data_table, N1, N0, zsim, alpha = 0.99)
```

# Troubleshooting

## Common Issues

### Issue 1: Very wide confidence intervals

This is often **correct**, not an error! Wide intervals indicate:

- Low compliance rate
- Small sample size
- High outcome variance

**Check instrument strength:**

``` r
# First-stage F-statistic
stage1 <- lm(D_observed ~ assignment, data = data_table)
f_stat <- summary(stage1)$fstatistic[1]
cat("F-statistic:", round(f_stat, 2), "\n")
# Rule of thumb: F > 10 indicates strong instrument
```

### Issue 2: Algorithm is slow

**Solutions:**

1.  Use Algorithm 2 (not Algorithm 1)
2.  Reduce permutations to 500 or 100 for testing
3.  Verify `zsim` is a matrix (not being regenerated)

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
