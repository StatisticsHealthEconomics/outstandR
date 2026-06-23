<div id="main" class="col-md-9" role="main">

# Simulating Target Cohorts with Custom Copulas

<div class="section level2">

## Introduction

In G-computation and MIM, `outstandR` simulates a synthetic
pseudo-population reflecting the ALD trial using a Gaussian copula. This
approach decouples the dependence structure (correlation matrix *Σ*)
from the marginal distributions of covariates.

The simulation follows three steps: 1. Draw samples from a multivariate
normal: *Z*_(*i*) ∼ *M**V**N*_(*K*)(0,*Σ*) 2. Apply the Probability
Integral Transform: *U*_(*i**k*) = *Φ*(*Z*_(*i**k*)) 3. Map to target
scales using inverse CDFs:
*X*_(*i**k*) = *F*_(*k*)⁻¹(*U*_(*i**k*);*θ*_(*k*))

First, let’s load the package and prepare the mixed continuous/binary
outcome data.

<div id="cb1" class="sourceCode">

``` r
library(outstandR)
data(AC_IPD_contY_mixedX) 
data(BC_ALD_contY_mixedX) 

# Define mixed outcome regression formula
lin_form_mixedX <- as.formula("y ~ X1 + X2 + X3 + trt + trt:(X1 + X2 + X4)")
```

</div>

------------------------------------------------------------------------

</div>

<div class="section level2">

## Specifying Marginal Distributions

By default, the package assumes a multivariate normal distribution.
However, users can define specific marginal distributions using the
`marginal_distns` argument. `outstandR` automatically retrieves reported
means and SDs from the ALD to parameterize these via the method of
moments.

Supported marginal distributions include: \* `"gamma"` (Gamma) \*
`"binom"` (Binomial) \* `"norm"` (Normal)

<div id="cb2" class="sourceCode">

``` r
custom_distns <- c(X1 = "gamma", X2 = "binom")

outstandR_custom <- outstandR(
  ipd_trial = AC_IPD_contY_mixedX,
  ald_trial = BC_ALD_contY_mixedX,
  strategy = strategy_gcomp_ml(
    formula = list(outcome_model = lin_form_mixedX),
    family = gaussian(link = "identity"),
    marginal_distns = custom_distns,
    N = 1000),
  seed = 12345
)

print(outstandR_custom)
#> Object of class 'outstandR' 
#> ITC algorithm: GCOMP_ML 
#> Model: gaussian 
#> Scale: mean_difference 
#> Common treatment: C 
#> Individual patient data study: A vs C 
#> Aggregate level data study: B vs C 
#> Confidence interval level: 0.95 
#> 
#> Contrasts:
#> 
#> # A tibble: 3 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl>      <dbl>      <dbl>
#> 1 AB           -0.293    0.0930     -0.891      0.305
#> 2 AC           -1.22     0.0619     -1.71      -0.732
#> 3 BC           -0.926    0.0311     -1.27      -0.581
#> 
#> Absolute:
#> 
#> # A tibble: 3 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl>      <dbl>      <dbl>
#> 1 A            -0.507    0.0176     -0.767     -0.247
#> 2 B            -0.214    0.0804     -0.769      0.342
#> 3 C             0.713    0.0493      0.277      1.15
```

</div>

------------------------------------------------------------------------

</div>

<div class="section level2">

## Manual Parameter Specification

For sensitivity analyses, or when the method of moments conversion isn’t
natively supported, you can manually override target parameters using
`marginal_params`.

For instance, to simulate an older cohort (mean = 70, SD = 10.5) for
covariate `X1` using a Gamma distribution:

<div id="cb3" class="sourceCode">

``` r
target_mean <- 70
target_sd <- 10.5

sens_params <- list(
  X1 = list(shape = (target_mean/target_sd)^2,
            rate = target_mean/(target_sd^2))
)

outstandR_sens <- outstandR(
  ipd_trial = AC_IPD_contY_mixedX,
  ald_trial = BC_ALD_contY_mixedX,
  strategy = strategy_gcomp_ml(
    formula = list(outcome_model = lin_form_mixedX),
    family = gaussian(link = "identity"),
    marginal_distns = c(X1 = "gamma"),
    marginal_params = sens_params,
    N = 1000),
  seed = 12345
)

print(outstandR_sens)
#> Object of class 'outstandR' 
#> ITC algorithm: GCOMP_ML 
#> Model: gaussian 
#> Scale: mean_difference 
#> Common treatment: C 
#> Individual patient data study: A vs C 
#> Aggregate level data study: B vs C 
#> Confidence interval level: 0.95 
#> 
#> Contrasts:
#> 
#> # A tibble: 3 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl>      <dbl>      <dbl>
#> 1 AB          -22.8   1073.         -87.1      41.4  
#> 2 AC          -23.8   1073.         -88.0      40.4  
#> 3 BC           -0.926    0.0311      -1.27     -0.581
#> 
#> Absolute:
#> 
#> # A tibble: 3 × 5
#>   Treatments Estimate Std.Error lower.0.95 upper.0.95
#>   <chr>         <dbl>     <dbl>      <dbl>      <dbl>
#> 1 A              62.1      260.       30.5       93.8
#> 2 B              85.0      853.       27.7      142. 
#> 3 C              85.9      853.       28.7      143.
```

</div>

</div>

</div>
