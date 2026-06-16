# Prepare Covariate Distributions

Resolves missing distributions and parameters by looking at the ALD.
Allows for partial specification (e.g., user specifies one variable,
function auto-detects the rest).

## Usage

``` r
prepare_covariate_distns(
  formula,
  ald,
  trt_var,
  marginal_distns,
  marginal_params,
  verbose = FALSE
)
```
