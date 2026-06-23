<div id="main" class="col-md-9" role="main">

# Bayesian G-computation using Stan

<div class="ref-description section level2">

Calculate draws of binary responses from posterior predictive
distribution from the Bayesian G-computation method using Hamiltonian
Monte Carlo.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
calc_gcomp_bayes(strategy, analysis_params, ...)
```

</div>

</div>

<div class="section level2">

## Arguments

-   strategy:

    A list specifying the model strategy, including:

    -   `outcome_model`: A linear regression `formula` object.

    -   `family`: A `family` object specifying the distribution and link
        function (e.g., `binomial`).

    -   `iter`: Number of iterations for the MCMC sampling.

    -   `warmup`: Number of warmup iterations for the MCMC sampling.

    -   `chains`: Number of MCMC chains.

-   analysis_params:

    List of analysis parameters. Must contain `ipd` and `ald`.

-   ...:

    Additional arguments passed to `rstanarm::stan_glm()`.

</div>

<div class="section level2">

## Value

A list containing:

-   `means`: A list containing:

    -   Posterior means for comparator treatment group.

    -   Posterior means for reference treatment group.

-   `model`: A list containing the `fit` object (from `stan_glm`),
    `rho`, `N`, and `stan_args`.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
strategy <- list(
  outcome_model = y ~ trt:X1,
  family = binomial(),
  rho = NA,
  N = 1000L,
  marginal_distns = NA,
  marginal_params = NA,
  trt_var = "trt",
  iter = 2000,
  warmup = 500,
  chains = 4)

ipd <- data.frame(
   trt = sample(c("A", "C"), size = 100, replace = TRUE),
   X1 = rnorm(100, 1, 1),
   y = sample(c(1,0), size = 100, prob = c(0.7, 0.3), replace = TRUE))

ald <- data.frame(
  trt = c(NA, NA, "B", "C", "B", "C"),
  variable = c("X1", "X1", "y", "y", NA, NA),
  statistic = c("mean", "sd", "sum", "sum", "N", "N"),
  value = c(0.5, 0.1, 10, 12, 20, 25))

res <- 
  calc_gcomp_bayes(
    strategy,
    analysis_params = list(
      ipd = ipd, ald = ald, 
      ref_trt = "C",
      ipd_comp = "A"))

str(res, max.level = 2, list.len = 3, vec.len = 2)
#> List of 3
#>  $ means          :List of 2
#>   ..$ A: num [1:2000] 0.626 0.69 ...
#>   ..$ C: num [1:2000] 0.61 0.667 ...
#>  $ point_estimates:List of 2
#>   ..$ A: num 0.603
#>   ..$ C: num 0.586
#>  $ model          :List of 4
#>   ..$ fit      :List of 28
#>   .. ..- attr(*, "class")= chr [1:3] "stanreg" "glm" ...
#>   ..$ rho      : logi NA
#>   ..$ N        : int 1000
#>   .. [list output truncated]
```

</div>

</div>

</div>
