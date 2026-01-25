# G-computation Maximum Likelihood Log-Odds Ratio

Marginal *A* vs *C* log-odds ratio (mean difference in expected
log-odds) estimated by transforming from probability to linear predictor
scale.

## Usage

``` r
gcomp_ml_log_odds_ratio(formula, ipd, ald)
```

## Arguments

- formula:

  Linear regression `formula` object

- ipd:

  Individual-level data

- ald:

  Aggregate-level data

## Value

Difference of log-odds

## Details

\\\log(\hat{\mu}\_A/(1 - \hat{\mu}\_A)) - \log(\hat{\mu}\_C/(1 -
\hat{\mu}\_C))\\

## See also

[`strategy_gcomp_ml()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md),
[`gcomp_ml.boot()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/gcomp_ml.boot.md)
