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
using imputation and combining estimates using pooling.

## Examples

``` r
strategy <- strategy_maic(
  formula = list(outcome_model = y ~ trt,  # default 
                 balance_model = ~ X1),
  family = binomial())
#> Treatment is guessed as: trt

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
#> [1] 0.0404359
#> 
#> $contrasts$var
#> [1] 0.2441335
#> 
#> 
#> $absolute
#> $absolute$mean
#>         A         C 
#> 0.5667870 0.5562863 
#> 
#> $absolute$var
#>           A           C 
#> 0.006112309 0.007795801 
#> 
#> 
#> $method_name
#> [1] "MAIC"
#> 
#> $model
#> $model$weights
#>   weights1   weights2   weights3   weights4   weights5   weights6   weights7 
#>  0.5097557  0.9862214  0.1767722  1.0189126  0.4083001  0.3004525  0.5227465 
#>   weights8   weights9  weights10  weights11  weights12  weights13  weights14 
#>  0.2650783  0.5314778  0.8255845  0.4200238  1.9246909  1.6187657  0.4651878 
#>  weights15  weights16  weights17  weights18  weights19  weights20  weights21 
#>  0.4238390  0.4253863  2.4215912  2.1405183  2.2258621  1.3013718  1.3575855 
#>  weights22  weights23  weights24  weights25  weights26  weights27  weights28 
#>  0.2828358  1.1742708  1.9678484  0.4041438  0.4765942  1.0933591  0.3810760 
#>  weights29  weights30  weights31  weights32  weights33  weights34  weights35 
#>  0.5973394  1.4076416  0.8244193  0.3731263  0.7741472  0.9083386  0.1294497 
#>  weights36  weights37  weights38  weights39  weights40  weights41  weights42 
#>  0.5622625  0.5166102  1.6685593  0.2084445  0.4586885  0.9201768  0.3863122 
#>  weights43  weights44  weights45  weights46  weights47  weights48  weights49 
#>  0.5076362  0.5611625  0.4537915  0.3643097  0.3339151  0.3729876  1.1400251 
#>  weights50  weights51  weights52  weights53  weights54  weights55  weights56 
#>  0.8037654  0.7215501  1.6619502  0.1326304  0.6915303  1.1697328  1.9035924 
#>  weights57  weights58  weights59  weights60  weights61  weights62  weights63 
#>  1.1017251  0.3539011  1.4816379  0.5736033  1.3541201  0.7266439  1.5106265 
#>  weights64  weights65  weights66  weights67  weights68  weights69  weights70 
#>  0.4651451  0.3815356  0.9272789  0.6387655  0.2490713  0.4498635  0.8699813 
#>  weights71  weights72  weights73  weights74  weights75  weights76  weights77 
#>  0.4311347  0.6371576  0.5781595  1.2766576  0.8184790  1.6135860  1.6471549 
#>  weights78  weights79  weights80  weights81  weights82  weights83  weights84 
#>  0.8892978  0.2861711  0.8390482  1.1599676  0.7896686  0.6933948  1.1559473 
#>  weights85  weights86  weights87  weights88  weights89  weights90  weights91 
#>  0.8915050  1.0140840  1.2954799  0.6523693  0.4526894  0.7622167  1.0480156 
#>  weights92  weights93  weights94  weights95  weights96  weights97  weights98 
#>  0.5038243  0.1916380  0.1984988  0.1956611  2.0358831  1.1200094  0.5659745 
#>  weights99 weights100 
#>  0.6973233  0.8650603 
#> 
#> $model$ESS
#>      ESS 
#> 71.57563 
#> 
#> 
  
```
