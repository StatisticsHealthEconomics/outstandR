# Calculate individual-level patient data statistics

Computes mean and variance statistics for individual-level patient data
using various approaches, including Matching-Adjusted Indirect
Comparison (MAIC), Simulated Treatment Comparison (STC), and
G-computation via Maximum Likelihood Estimation (MLE) or Bayesian
inference.

## Usage

``` r
calc_IPD_stats(strategy, analysis_params, ...)

# Default S3 method
calc_IPD_stats(...)

# S3 method for class 'stc'
calc_IPD_stats(strategy, analysis_params, var_method = NULL, ...)

# S3 method for class 'maic'
calc_IPD_stats(strategy, analysis_params, var_method = NULL, ...)

# S3 method for class 'gcomp_ml'
calc_IPD_stats(strategy, analysis_params, var_method = NULL, ...)

# S3 method for class 'gcomp_bayes'
calc_IPD_stats(strategy, analysis_params, var_method = NULL, ...)

# S3 method for class 'mim'
calc_IPD_stats(strategy, analysis_params, var_method = NULL, ...)
```

## Arguments

- strategy:

  A list corresponding to different modelling approaches

- analysis_params:

  A list containing:

  - `ipd`: Individual-level patient data (data frame)

  - `ald`: Aggregate-level trial data (data frame)

  - `ref_trt`: Treatment label for the reference arm (common; e.g., "C")

  - `ipd_comp`: Treatment label for the comparator arm in the IPD (e.g.,
    "A")

  - `scale`: Scaling parameter ("log_odds", "risk_difference",
    "log_relative_risk")

- ...:

  Additional arguments

- var_method:

  A string specifying the variance estimation method, either "sample"
  (default) or "sandwich".

## Value

A list containing:

- `contrasts`: A list with elements `mean` and `var`.

- `absolute`: A list with elements `mean` and `var`.

## Simulated treatment comparison statistics

IPD for reference "C" and comparator "A" trial arms are used to fit a
regression model describing the observed outcomes \\y\\ in terms of the
relevant baseline characteristics \\x\\ and the treatment variable
\\z\\.

## Matching-adjusted indirect comparison statistics

Marginal IPD comparator treatment "A" vs reference treatment "C"
treatment effect estimates using bootstrapping sampling.

## G-computation maximum likelihood statistics

Compute a non-parametric bootstrap with default \\R=1000\\ resamples.

## G-computation Bayesian statistics

Using Stan, compute marginal relative effects for IPD comparator "A" vs
reference "C" treatment arms for each MCMC sample by transforming from
probability to linear predictor scale.

## Multiple imputation marginalisation

Using Stan, compute marginal relative treatment effect for IPD
comparator "A" vs reference "C" arms for each MCMC sample by
transforming from probability to linear predictor scale. Approximate by
using imputation and combining estimates using Rubin's rules.

## Examples

``` r
strategy <- strategy_maic(formula = as.formula(y~trt:X1), family = binomial())
ipd <- data.frame(trt = sample(c("A", "C"), size = 100, replace = TRUE),
                  X1 = rnorm(100, 1, 1),
                  y = sample(c(1,0), size = 100, prob = c(0.7,0.3), replace = TRUE))

ald <- data.frame(trt = c(NA, "B", "C", "B", "C"),
                  variable = c("X1", "y", "y", NA, NA),
                  statistic = c("mean", "sum", "sum", "N", "N"),
                  value = c(0.5, 10, 12, 20, 25))

calc_IPD_stats(strategy,
  analysis_params = list(ipd = ipd, ald = ald, scale = "log_odds"))
#> $contrasts
#> $contrasts$mean
#> [1] 0.04072276
#> 
#> $contrasts$var
#> [1] 0.2451418
#> 
#> 
#> $absolute
#> $absolute$mean
#>         A         C 
#> 0.7025638 0.6951590 
#> 
#> $absolute$var
#>           A           C 
#> 0.005650160 0.005183217 
#> 
#> 
#> $method_name
#> [1] "MAIC"
#> 
#> $model
#> $model$weights
#>   weights1   weights2   weights3   weights4   weights5   weights6   weights7 
#>  0.8644187  0.7471808  1.0941498  0.5137571  1.6671230  1.0104253  0.7662525 
#>   weights8   weights9  weights10  weights11  weights12  weights13  weights14 
#>  1.0924254  0.7196236  0.8194893  1.1114083  0.7223890  1.2919443  0.8837867 
#>  weights15  weights16  weights17  weights18  weights19  weights20  weights21 
#>  0.8534521  1.0946767  0.3083721  0.7779886  0.7720410  0.5818904  1.2813310 
#>  weights22  weights23  weights24  weights25  weights26  weights27  weights28 
#>  0.4048974  0.7521981  1.2995935  1.0658084  0.8159727  1.5211228  0.5130418 
#>  weights29  weights30  weights31  weights32  weights33  weights34  weights35 
#>  0.8678377  0.9656656  0.7451850  1.2546478  1.4252054  2.1407440  0.3410411 
#>  weights36  weights37  weights38  weights39  weights40  weights41  weights42 
#>  0.3434997  1.0843389  1.7071231  0.3736615  1.0078974  1.2780704  0.8389136 
#>  weights43  weights44  weights45  weights46  weights47  weights48  weights49 
#>  0.9167890  0.6745588  1.0534299  0.4039551  0.3274506  0.7164308  1.8885636 
#>  weights50  weights51  weights52  weights53  weights54  weights55  weights56 
#>  1.2057152  1.0639350  1.5633846  0.9396529  3.0833536  0.9397011  0.4908616 
#>  weights57  weights58  weights59  weights60  weights61  weights62  weights63 
#>  1.1619506  0.5912057  0.7738476  0.3515805  0.8976565  0.5447694  0.4509756 
#>  weights64  weights65  weights66  weights67  weights68  weights69  weights70 
#>  0.7542419  1.0597763  0.9247797  1.8938653  0.6103745  0.8746193  0.5307043 
#>  weights71  weights72  weights73  weights74  weights75  weights76  weights77 
#>  0.4544261  0.6079698  0.8191396  1.3460228  0.7095652  0.5928967  0.4365664 
#>  weights78  weights79  weights80  weights81  weights82  weights83  weights84 
#>  1.4783556  0.5532154  0.2758676  0.5065433  0.3944689  1.0978242  0.3406655 
#>  weights85  weights86  weights87  weights88  weights89  weights90  weights91 
#>  0.9603423  0.8679057  0.6058239  1.2072252  0.4784969  0.5336189  0.7415894 
#>  weights92  weights93  weights94  weights95  weights96  weights97  weights98 
#>  0.7947574  1.2495765  0.7779616  0.7617820  0.4187021  0.3652338  1.6012872 
#>  weights99 weights100 
#>  0.2849569  1.3295498 
#> 
#> $model$ESS
#>     ESS 
#> 79.4969 
#> 
#> 
  
```
