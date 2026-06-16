# Calculate the difference between treatments using all evidence

This is the main, top-level wrapper for `{outstandR}`. Methods taken
from (Remiro‐Azócar et al. 2022) .

## Usage

``` r
outstandR(AC.IPD, BC.ALD, strategy, CI = 0.95, ...)
```

## Arguments

- AC.IPD:

  Individual-level patient data. Suppose between studies *A* and *C*.

- BC.ALD:

  Aggregate-level data. Suppose between studies *B* and *C*.

- strategy:

  Computation strategy function. These can be
  [`strategy_maic()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md),
  [`strategy_stc()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md),
  [`strategy_gcomp_ml()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)
  and
  [`strategy_gcomp_bayes()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/strategy.md)

- CI:

  Confidence interval; between 0,1

- ...:

  Additional arguments

## Value

List of length 3 of statistics as a `outstandR` class object. Containing
statistics between each pair of treatments. These are the mean
contrasts, variances and confidence intervals, respectively.

## References

Remiro‐Azócar A, Heath A, Baio G (2022). “Parametric G‐computation for
compatible indirect treatment comparisons with limited individual
patient data.” *Res. Synth. Methods*, 1--31. ISSN 1759-2879,
[doi:10.1002/jrsm.1565](https://doi.org/10.1002/jrsm.1565) , 2108.12208.

## Examples

``` r
data(AC_IPD)  # AC patient-level data
data(BC_ALD)  # BC aggregate-level data

lin_form <- as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")

# matching-adjusted indirect comparison
outstandR_maic <- outstandR(AC_IPD, BC_ALD, strategy = strategy_maic(formula = lin_form))

# simulated treatment comparison
outstandR_stc <- outstandR(AC_IPD, BC_ALD, strategy = strategy_stc(lin_form))

# G-computation with maximum likelihood
# outstandR_gcomp_ml <- outstandR(AC_IPD, BC_ALD, strategy = strategy_gcomp_ml(lin_form))

# G-computation with Bayesian inference
outstandR_gcomp_bayes <- outstandR(AC_IPD, BC_ALD, strategy = strategy_gcomp_bayes(lin_form))
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 1).
#> Chain 1: 
#> Chain 1: Gradient evaluation took 7.1e-05 seconds
#> Chain 1: 1000 transitions using 10 leapfrog steps per transition would take 0.71 seconds.
#> Chain 1: Adjust your expectations accordingly!
#> Chain 1: 
#> Chain 1: 
#> Chain 1: Iteration:    1 / 4000 [  0%]  (Warmup)
#> Chain 1: Iteration:  400 / 4000 [ 10%]  (Warmup)
#> Chain 1: Iteration:  800 / 4000 [ 20%]  (Warmup)
#> Chain 1: Iteration: 1200 / 4000 [ 30%]  (Warmup)
#> Chain 1: Iteration: 1600 / 4000 [ 40%]  (Warmup)
#> Chain 1: Iteration: 2000 / 4000 [ 50%]  (Warmup)
#> Chain 1: Iteration: 2001 / 4000 [ 50%]  (Sampling)
#> Chain 1: Iteration: 2400 / 4000 [ 60%]  (Sampling)
#> Chain 1: Iteration: 2800 / 4000 [ 70%]  (Sampling)
#> Chain 1: Iteration: 3200 / 4000 [ 80%]  (Sampling)
#> Chain 1: Iteration: 3600 / 4000 [ 90%]  (Sampling)
#> Chain 1: Iteration: 4000 / 4000 [100%]  (Sampling)
#> Chain 1: 
#> Chain 1:  Elapsed Time: 0.71 seconds (Warm-up)
#> Chain 1:                0.733 seconds (Sampling)
#> Chain 1:                1.443 seconds (Total)
#> Chain 1: 
#> 
#> SAMPLING FOR MODEL 'bernoulli' NOW (CHAIN 2).
#> Chain 2: 
#> Chain 2: Gradient evaluation took 2.2e-05 seconds
#> Chain 2: 1000 transitions using 10 leapfrog steps per transition would take 0.22 seconds.
#> Chain 2: Adjust your expectations accordingly!
#> Chain 2: 
#> Chain 2: 
#> Chain 2: Iteration:    1 / 4000 [  0%]  (Warmup)
#> Chain 2: Iteration:  400 / 4000 [ 10%]  (Warmup)
#> Chain 2: Iteration:  800 / 4000 [ 20%]  (Warmup)
#> Chain 2: Iteration: 1200 / 4000 [ 30%]  (Warmup)
#> Chain 2: Iteration: 1600 / 4000 [ 40%]  (Warmup)
#> Chain 2: Iteration: 2000 / 4000 [ 50%]  (Warmup)
#> Chain 2: Iteration: 2001 / 4000 [ 50%]  (Sampling)
#> Chain 2: Iteration: 2400 / 4000 [ 60%]  (Sampling)
#> Chain 2: Iteration: 2800 / 4000 [ 70%]  (Sampling)
#> Chain 2: Iteration: 3200 / 4000 [ 80%]  (Sampling)
#> Chain 2: Iteration: 3600 / 4000 [ 90%]  (Sampling)
#> Chain 2: Iteration: 4000 / 4000 [100%]  (Sampling)
#> Chain 2: 
#> Chain 2:  Elapsed Time: 0.73 seconds (Warm-up)
#> Chain 2:                0.817 seconds (Sampling)
#> Chain 2:                1.547 seconds (Total)
#> Chain 2: 
```
