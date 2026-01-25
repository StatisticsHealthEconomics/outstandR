# Calculate individual-level patient data statistics

Computes mean and variance statistics for individual-level patient data
using various approaches, including Matching-Adjusted Indirect
Comparison (MAIC), Simulated Treatment Comparison (STC), and
G-computation via Maximum Likelihood Estimation (MLE) or Bayesian
inference.

## Usage

``` r
IPD_stats(strategy, ipd, ald, scale, ...)

# Default S3 method
IPD_stats(...)

# S3 method for class 'mim'
IPD_stats(strategy, ipd, ald, scale, ...)

# S3 method for class 'stc'
IPD_stats(strategy, ipd, ald, scale, ...)

# S3 method for class 'maic'
IPD_stats(strategy, ipd, ald, scale, ...)

# S3 method for class 'gcomp_ml'
IPD_stats(strategy, ipd, ald, scale, ...)

# S3 method for class 'gcomp_bayes'
IPD_stats(strategy, ipd, ald, scale, ...)
```

## Arguments

- strategy:

  A list corresponding to different approaches

- ipd:

  Individual-level data

- ald:

  Aggregate-level data

- scale:

  A scaling parameter for the effect calculation.

- ...:

  Additional arguments

## Value

A list containing:

- mean:

  Estimated mean treatment effect.

- var:

  Estimated variance of the treatment effect.

## Multiple imputation marginalisation

Using Stan, compute marginal relative treatment effect for *A* vs *C*
for each MCMC sample by transforming from probability to linear
predictor scale. Approximate by using imputation and combining estimates
using Rubin's rules, in contrast to `IPD_stats.gcomp_bayes()`.

## Simulated treatment comparison statistics

IPD from the *AC* trial are used to fit a regression model describing
the observed outcomes \\y\\ in terms of the relevant baseline
characteristics \\x\\ and the treatment variable \\z\\.

## Matching-adjusted indirect comparison statistics

Marginal *A* vs *C* treatment effect estimates using bootstrapping
sampling.

## G-computation maximum likelihood statistics

Compute a non-parametric bootstrap with default \\R=1000\\ resamples.

## G-computation Bayesian statistics

Using Stan, compute marginal log-odds ratio for *A* vs *C* for each MCMC
sample by transforming from probability to linear predictor scale.

## Examples

``` r
if (FALSE) { # \dontrun{
strategy <- strategy_maic()
ipd <- data.frame(id = 1:100, treatment = sample(c("A", "C"), 100, replace = TRUE), outcome = rnorm(100))
ald <- data.frame(treatment = c("A", "C"), mean = c(0.2, 0.1), var = c(0.05, 0.03))
IPD_stats(strategy, ipd, ald, scale = "log")
} # }
```
