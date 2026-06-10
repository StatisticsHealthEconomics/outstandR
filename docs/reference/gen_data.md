# Generate simulated datasets of IPD covariates and binary outcome for a trial

Generate simulated datasets of IPD covariates and binary outcome for a
trial

## Usage

``` r
gen_data(
  N,
  b_trt,
  b_X,
  b_EM,
  b_0,
  meanX_EM,
  sdX_EM,
  meanX,
  sdX,
  corX,
  allocation,
  family = binomial("logit")
)
```

## Arguments

- N:

  Total number of patients

- b_trt:

  `b` coefficient for active treatment vs. common comparator

- b_X:

  `b` coefficients for each prognostic variable `X`

- b_EM:

  `b` coefficients effect modifiers

- b_0:

  Intercept coefficient

- meanX:

  Mean of each normally-distributed covariate `X`

- sdX:

  Standard deviation of each covariate `X`

- corX:

  Covariate correlation coefficient of `X`

- allocation:

  Allocation to active treatment as proportion of total; 0 to 1

- event_rate:

  Event rate

## Value

Data frame of `X`, `trt` and `y`

## Examples

``` r

if (FALSE) { # \dontrun{
x <- gen_data(
 N = 100,
 b_trt = log(0.17),
 b_X = -log(0.5),
 b_EM = -log(0.67),
 b_0 = -0.62,
 meanX = 0.6,
 sdX = 0.4,
 event_rate = 0.35, 
 corX = 0.2,
 allocation = 2/3) 

head(x)
} # }
```
