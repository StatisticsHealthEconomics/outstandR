<div id="main" class="col-md-9" role="main">

# Simulate Aggregate-Level Data Pseudo-Population

<div class="ref-description section level2">

Generates a synthetic cohort using a normal copula based on
aggregate-level data.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
simulate_ALD_pseudo_pop(
  formula,
  ipd = NULL,
  ald = NULL,
  trt_var,
  rho = NA,
  N = 1000,
  marginal_distns = NA,
  marginal_params = NA,
  seed = NULL,
  verbose = FALSE
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   formula:

    Linear regression `formula` object. Prognostic factors (PF) are main
    effects and effect modifiers (EM) are interactions with the
    treatment variable, e.g., y \~ X1 + trt + trt:X2. For covariates as
    both PF and EM use `*` syntax.

-   ipd:

    Individual-level patient data. Dataframe with one row per patient
    with outcome, treatment and covariate columns.

-   ald:

    Aggregate-level data. Long format summary statistics for each
    covariate and treatment outcomes. We assume a common distribution
    for each treatment arm.

-   rho:

    A named square matrix of covariate correlations or single value;
    default NA takes from IPD.

-   N:

    Sample size for the synthetic cohort. Default is 1000.

-   marginal_distns:

    Marginal distributions names; vector default NA. Available
    distributions are given in stats::Distributions. See
    `copula::Mvdc()` for details

-   marginal_params:

    Marginal distributions parameters; named list of lists, default NA.
    See `copula::Mvdc()` for details.

-   seed:

    Random seed

-   verbose:

    Default `FALSE`

</div>

<div class="section level2">

## Value

A data frame representing the synthetic pseudo-population. It contains
`N` rows (one for each simulated individual) and columns for every
covariate specified in `marginal_distns` of `formula`.

</div>

<div class="section level2">

## Examples

<div class="sourceCode">

``` r
if (FALSE) { # \dontrun{

## Example 1: Simulating data with explicitly defined marginals and
## a provided correlation matrix (rho)

my_formula <- y ~ x1 + x2*trt + trt

# Define marginal distributions (Normal for age, Lognormal for BMI)
dists <- c(x1 = "norm", x2 = "norm")

# Define parameters matching the chosen distributions
params <- list(
  x1 = list(mean = 60, sd = 8),
  x2 = list(mean = 3, sd = 0.1)
)

# Define a 2x2 correlation matrix
corr_mat <- matrix(c(1, 0.25,
                     0.25, 1), nrow = 2,
                   dimnames = list(c("x1", "x2"),
                                   c("x1", "x2")))

# Generate the synthetic cohort
sim_cohort_marginals_rho <- simulate_ALD_pseudo_pop(
  formula = my_formula,
  ald = mock_ald,
  trt_var = "trt",
  rho = corr_mat,
  N = 100,
  marginal_distns = dists,
  marginal_params = params
)

head(sim_cohort_marginals_rho)

## Example 2: Estimating the correlation matrix (rho) from provided IPD

# Create some mock Individual Patient Data (IPD)
mock_ipd <- data.frame(
  x1 = rnorm(200, 60, 8),
  x2 = rlnorm(200, 3.3, 0.1)
)

mock_ald <- data.frame(
  variable = c("x1", "x1", "x2", "x2"),
  statistic = c("mean", "sd", "mean", "sd"),
  value = c(1000, 500, 0.7, 0.1),
  trt = NA
)

# Generate the synthetic cohort using IPD for the correlation structure
sim_cohort <- simulate_ALD_pseudo_pop(
  formula = my_formula,
  ipd = mock_ipd,
  ald = mock_ald,
  trt_var = "trt",
  rho = NA,
  N = 100
)

head(sim_cohort)

# Generate the synthetic cohort using IPD for the correlation structure
# with marginals
sim_cohort_marginals <- simulate_ALD_pseudo_pop(
  formula = my_formula,
  ipd = mock_ipd,
  ald = mock_ald,
  trt_var = "trt",
  N = 100,
  marginal_distns = dists,
  marginal_params = params
)

head(sim_cohort_marginals)
} # }
```

</div>

</div>

</div>
