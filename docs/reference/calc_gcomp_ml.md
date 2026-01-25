# G-computation Maximum Likelihood Bootstrap

Computes the mean difference in treatment effects using bootstrap
resampling.

## Usage

``` r
calc_gcomp_ml(strategy, analysis_params)
```

## Arguments

- strategy:

  A list specifying the model strategy, including:

  - `R`: Number of bootstrap replications.

  - `formula`: A linear regression `formula` object.

  - `family`: A `family` object specifying the distribution and link
    function (e.g., `binomial`).

  - `N`: Synthetic sample size for g-computation.

- analysis_params:

  List of analysis parameters.

## Value

A list containing:

- `means`: A list containing:

  - `A`: Bootstrap estimates for comparator treatment group "A".

  - `C`: Bootstrap estimates for reference treatment group "C".

- `model`: A list containing the `fit` object, `rho`, and `N`.

## Examples

``` r
strategy <- list(
  formula = y ~ trt:X1,
  family = binomial(),
  rho = NA,
  N = 1000L,
  n_boot = 100L,
  marginal_distns = NA,
  marginal_params = NA,
  trt_var = "trt")

ipd <- data.frame(trt = sample(c("A", "C"), size = 100, replace = TRUE),
                  X1 = rnorm(100, 1, 1),
                  y = sample(c(1,0), size = 100, prob = c(0.7,0.3), replace = TRUE))

ald <- data.frame(trt = c(NA, NA, "B", "C", "B", "C"),
                  variable = c("X1", "X1", "y", "y", NA, NA),
                  statistic = c("mean", "sd", "sum", "sum", "N", "N"),
                  value = c(0.5, 0.1, 10, 12, 20, 25))

calc_gcomp_ml(
  strategy,
  analysis_params = 
    list(ipd = ipd, ald = ald, 
         ref_trt = "C", 
         ipd_comp = "A"))
#> $means
#> $means$A
#>   [1] 0.7537442 0.7804083 0.7803788 0.8176517 0.6261085 0.7942802 0.7799355
#>   [8] 0.7591019 0.6518243 0.7312572 0.7541696 0.7059253 0.7811411 0.8268642
#>  [15] 0.6940952 0.8258627 0.7403065 0.6394265 0.7602504 0.8284318 0.7329183
#>  [22] 0.6752033 0.7460324 0.8146492 0.7374148 0.6383276 0.8390803 0.7564710
#>  [29] 0.7378948 0.7886655 0.6853042 0.6716206 0.7932023 0.6966904 0.6595554
#>  [36] 0.7932463 0.7112485 0.6219996 0.7542252 0.8355904 0.6754437 0.6716659
#>  [43] 0.6411188 0.7887678 0.7710465 0.7059794 0.7424005 0.7184483 0.7249245
#>  [50] 0.7017129 0.6366139 0.7383667 0.5994948 0.6743979 0.6765540 0.7350509
#>  [57] 0.7290709 0.7415034 0.6902279 0.7392571 0.7830743 0.7346571 0.5640236
#>  [64] 0.7367141 0.6906191 0.7084894 0.8354787 0.6757976 0.7226717 0.7427268
#>  [71] 0.7026815 0.6695894 0.7676483 0.7057674 0.6270349 0.7571799 0.7454205
#>  [78] 0.7119828 0.6421126 0.7369034 0.6467641 0.7610800 0.6804569 0.6840998
#>  [85] 0.6859370 0.7995461 0.8335337 0.6359600 0.7665221 0.8142967 0.6147780
#>  [92] 0.7877127 0.6851999 0.6996645 0.6995421 0.7771037 0.6883149 0.7256025
#>  [99] 0.6940002 0.6910828
#> 
#> $means$C
#>   [1] 0.7341149 0.7348943 0.7775082 0.7798873 0.6098913 0.6659045 0.7434782
#>   [8] 0.7479637 0.5933621 0.6887711 0.7349731 0.6764418 0.7076330 0.8007199
#>  [15] 0.6331876 0.7668043 0.7324601 0.6506011 0.7039040 0.7182162 0.7227844
#>  [22] 0.6894363 0.6942953 0.8117161 0.7185000 0.6587124 0.8214528 0.7160370
#>  [29] 0.7111210 0.7293131 0.6978250 0.6663719 0.7957795 0.6669556 0.6180348
#>  [36] 0.7828916 0.6944620 0.6809854 0.7553352 0.8052713 0.6937876 0.6060482
#>  [43] 0.6881049 0.7252061 0.7622770 0.6958937 0.7389865 0.7819663 0.6823700
#>  [50] 0.6877324 0.6417942 0.6919524 0.6025808 0.6760991 0.6483074 0.7345234
#>  [57] 0.7308785 0.5999942 0.6087972 0.7532254 0.8210469 0.7336323 0.5651038
#>  [64] 0.7172907 0.7437540 0.6492722 0.7409792 0.6741772 0.6881777 0.7484329
#>  [71] 0.6400894 0.7261319 0.7232688 0.6698483 0.6241309 0.7600538 0.7707850
#>  [78] 0.6990852 0.7179687 0.7009608 0.6804860 0.7713426 0.6878557 0.6818330
#>  [85] 0.6924526 0.7734790 0.8151975 0.6606551 0.7645366 0.7548032 0.6118537
#>  [92] 0.7292321 0.6947536 0.6100965 0.6881375 0.8108806 0.6865149 0.6980783
#>  [99] 0.6774447 0.5844859
#> 
#> 
#> $model
#> $model$fit
#> 
#> Call:  glm(formula = formula, family = family, data = ipd)
#> 
#> Coefficients:
#> (Intercept)      trtA:X1      trtC:X1  
#>     0.85293      0.16062     -0.03278  
#> 
#> Degrees of Freedom: 99 Total (i.e. Null);  97 Residual
#> Null Deviance:       120.4 
#> Residual Deviance: 120   AIC: 126
#> 
#> $model$rho
#> [1] NA
#> 
#> $model$N
#> [1] 1000
#> 
#> 
         
```
