# Strategy class and subclasses

The `strategy` class is a virtual class that defines the statistical
approach for population adjustment in indirect treatment comparisons
These objects are constructors that validate hyperparameters and
encapsulate modelling settings before execution by
[`outstandR()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/outstandR.md)

## Details

Objects of class `strategy` have a common structure but carry different
subclasses to trigger specific S3 method dispatch

- General fields:

  Shared by all strategies:

  - `formula`: The linear regression formula for the outcome model

  - `family`: A base R family object specifying the distribution and
    link

  - `trt_var`: The name of the treatment variable.

- maic subclass:

  Additional fields for Matching-Adjusted Indirect Comparison:

  - `n_boot`: Number of bootstrap resamples for variance estimation.

- stc subclass:

  Additional fields for Simulated Treatment Comparison:

  - `N`: Synthetic sample size for the target population.

- gcomp_ml subclass:

  Additional fields for Maximum Likelihood G-computation:

  - `rho`: Named square matrix of covariate correlations.

  - `marginal_distns`: Names of the marginal distributions for
    covariates.

  - `marginal_params`: Parameters for the marginal distributions.

  - `N`: Synthetic sample size for the pseudo-population.

  - `n_boot`: Number of bootstrap resamples.

- gcomp_bayes subclass:

  Additional fields for Bayesian G-computation:

  - `rho`, `marginal_distns`, `marginal_params`, `N`: Same as
    `gcomp_ml`.

  - `...`: Additional arguments passed to the Stan engine via
    [`rstanarm::stan_glm()`](https://mc-stan.org/rstanarm/reference/stan_glm.html).

- mim subclass:

  Additional fields for Multiple Imputation Marginalization:

  - `rho`: Correlation matrix.

  - `N`: Number of iterations/simulated individuals.
