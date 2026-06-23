<div id="main" class="col-md-9" role="main">

# G-computation Maximum Likelihood Bootstrap

<div class="ref-description section level2">

Computes the mean difference in treatment effects using bootstrap
resampling.

</div>

<div class="section level2">

## Usage

<div class="sourceCode">

``` r
calc_gcomp_ml(strategy, analysis_params)
```

</div>

</div>

<div class="section level2">

## Arguments

-   strategy:

    A list specifying the model strategy, including:

    -   `R`: Number of bootstrap replications.

    -   `outcome_model`: A linear regression `formula` object for the
        outcome model.

    -   `family`: A `family` object specifying the distribution and link
        function (e.g., `binomial`).

    -   `N`: Synthetic sample size for g-computation.

-   analysis_params:

    List of analysis parameters.

</div>

<div class="section level2">

## Value

A list containing:

-   `means`: A list containing:

    -   `A`: Bootstrap estimates for comparator treatment group "A".

    -   `C`: Bootstrap estimates for reference treatment group "C".

-   `model`: A list containing the `fit` object, `rho`, and `N`.

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
#>   [1] 0.7524448 0.7208018 0.6857046 0.7773186 0.7379123 0.7953863 0.7416379
#>   [8] 0.7333357 0.6740201 0.7496222 0.7004516 0.7565848 0.7328825 0.7406903
#>  [15] 0.7488645 0.7340921 0.6933722 0.7124750 0.8259316 0.8126080 0.7030263
#>  [22] 0.7026071 0.7626186 0.7341363 0.7593515 0.7545032 0.8159394 0.7371818
#>  [29] 0.7284083 0.6787501 0.6071921 0.7056624 0.7314274 0.6486667 0.7948652
#>  [36] 0.8041567 0.6859274 0.7155094 0.7249215 0.7912577 0.6099682 0.6655606
#>  [43] 0.8584396 0.6814494 0.7191635 0.7445806 0.6603892 0.7538379 0.7519644
#>  [50] 0.7192660 0.7830961 0.8151415 0.6831753 0.7150227 0.8449769 0.7668341
#>  [57] 0.7483105 0.7879836 0.7614439 0.8144559 0.6243036 0.6539474 0.7320960
#>  [64] 0.7168957 0.7408405 0.7626611 0.7756132 0.7237764 0.7889745 0.7579104
#>  [71] 0.7375602 0.7413274 0.6781683 0.6714261 0.7565422 0.8332126 0.6637191
#>  [78] 0.5921426 0.7549736 0.7339375 0.7560676 0.7910103 0.7243667 0.7994549
#>  [85] 0.8220225 0.8462224 0.7428240 0.5864449 0.7003216 0.7453676 0.8061360
#>  [92] 0.8159053 0.6917901 0.6791657 0.7111803 0.7263968 0.7310195 0.7909051
#>  [99] 0.7227581 0.7323826
#> 
#> $means$C
#>   [1] 0.7309018 0.7420960 0.6819819 0.6978632 0.7191187 0.7537383 0.7023496
#>   [8] 0.7228270 0.6564365 0.7003893 0.6392001 0.6893771 0.6915394 0.7598735
#>  [15] 0.6611513 0.7114534 0.6330323 0.6893841 0.8136764 0.8007598 0.7033007
#>  [22] 0.7148683 0.7531519 0.7223997 0.7248275 0.7375199 0.7290035 0.7066567
#>  [29] 0.7381150 0.7441146 0.5489844 0.6454383 0.7261025 0.6269568 0.7808563
#>  [36] 0.7845551 0.6790214 0.7108215 0.7857677 0.7722908 0.5920379 0.6711567
#>  [43] 0.7735766 0.7082050 0.7105959 0.6968303 0.6837809 0.7379012 0.7804497
#>  [50] 0.7257683 0.7395591 0.7641128 0.5839501 0.7344368 0.7883609 0.7415819
#>  [57] 0.7914461 0.7409043 0.7648980 0.7197088 0.6058059 0.6712301 0.7169338
#>  [64] 0.6840646 0.7066697 0.6927932 0.7584608 0.7285138 0.7396870 0.7564756
#>  [71] 0.7302828 0.6995002 0.6820383 0.6638061 0.7293550 0.8527352 0.6752132
#>  [78] 0.6727563 0.7111311 0.6847701 0.7357006 0.7539906 0.7102805 0.7856127
#>  [85] 0.8156423 0.8199455 0.7346005 0.6136072 0.7113708 0.7010819 0.7942727
#>  [92] 0.7816347 0.6858805 0.6637144 0.7716140 0.7212641 0.7290398 0.7679028
#>  [99] 0.6818292 0.7198914
#> 
#> 
#> $point_estimates
#> $point_estimates$A
#>         1 
#> 0.7321823 
#> 
#> $point_estimates$C
#>         0 
#> 0.7163028 
#> 
#> 
#> $model
#> $model$fit
#> 
#> Call:  glm(formula = outcome_model, family = family, data = ipd)
#> 
#> Coefficients:
#> (Intercept)      trtA:X1      trtC:X1  
#>     0.93459      0.14321     -0.01688  
#> 
#> Degrees of Freedom: 99 Total (i.e. Null);  97 Residual
#> Null Deviance:       116.7 
#> Residual Deviance: 116.3     AIC: 122.3
#> 
#> $model$rho
#> [1] NA
#> 
#> $model$N
#> [1] 1000
#> 
#> $model$n_boot
#> [1] 100
#> 
#> 
         
```

</div>

</div>

</div>
