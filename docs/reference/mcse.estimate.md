# Monte Carlo SE of any continuous performance metric

Monte Carlo SE of any continuous performance metric

## Usage

``` r
mcse.estimate(pm)
```

## Arguments

- pm:

  pm

## Value

\\sqrt(sum((pm - pm_mean)^2)/(nsim\*(nsim-1)))\\
