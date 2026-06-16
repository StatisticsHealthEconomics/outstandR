# G-computation using Stan

Calculate draws of binary responses from posterior predictive
distribution from the Bayesian G-computation method using Hamiltonian
Monte Carlo.

## Usage

``` r
calc_gcomp_bayes(strategy, analysis_params, ...)
```

## Arguments

- strategy:

  A list specifying the model strategy, including:

  formula

  :   A linear regression `formula` object.

  family

  :   A `family` object specifying the distribution and link function
      (e.g., `binomial`).

  iter

  :   Number of iterations for the MCMC sampling.

  warmup

  :   Number of warmup iterations for the MCMC sampling.

  chains

  :   Number of MCMC chains.

- ipd:

  Individual-level patient data. Dataframe with one row per patient with
  outcome, treatment and covariate columns.

- ald:

  Aggregate-level data. Single row matrix with summary statistics for
  each covariate and treatment outcomes. The format is 'mean.*' and
  'sd.*' for covariates and 'y.*.sum', 'y.*.bar', 'y.\*.sd' for
  treatments B and C. We assume a common distribution for each treatment
  arm.

## Value

A list of \\y^\*\_A\\ and \\y^\*\_C\\ posterior predictions:

- `` `0` ``:

  Posterior means for reference treatment group "C".

- `` `1` ``:

  Posterior means for comparator treatment group "A".

## Examples

``` r
if (FALSE) { # \dontrun{
strategy <- list(
  formula = y ~ trt + age,
  family = binomial(),
  iter = 2000,
  warmup = 500,
  chains = 4
)
ipd <- data.frame(trt = c("A", "C"),
                  y = c(1, 0),
                  age = c(30, 40))
ald <- data.frame()
calc_gcomp_bayes(strategy, ipd, ald)
} # }
```
