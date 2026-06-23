<div id="main" class="col-md-9" role="main">

# G-computation maximum likelihood mean outcomes by arm

<div class="ref-description section level2">

G-computation maximum likelihood mean outcomes by arm

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
gcomp_ml_means(
  outcome_model,
  family,
  ipd,
  ald,
  trt_var,
  rho = NA,
  N = 1000,
  ref_trt,
  comp_trt,
  marginal_distns = NA,
  marginal_params = NA
)
```

</div>

</div>

<div class="section level2">

## Arguments

-   family:

    A 'family' object specifying the distribution and link function
    (e.g., 'binomial'). See stats::family() for more details.

-   ipd:

    Individual-level patient data. Dataframe with one row per patient
    with outcome, treatment and covariate columns.

-   ald:

    Aggregate-level data. Long format summary statistics for each
    covariate and treatment outcomes. We assume a common distribution
    for each treatment arm.

-   rho:

    A named square matrix of covariate correlations; default NA.

-   N:

    Synthetic sample size for g-computation

-   marginal_distns, marginal_params:

    Marginal distributions and parameters

</div>

<div class="section level2">

## Value

A list containing:

-   `stats`: Named vector of marginal mean probabilities

-   `model`: The fitted glm object

</div>

<div class="section level2">

## See also

<div class="dont-index">

`strategy_gcomp_ml()` `gcomp_ml.boot()`

</div>

</div>

</div>
