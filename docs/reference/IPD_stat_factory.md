# Factory function for creating calc_IPD_stats methods

Creates a method for computing IPD mean and variance statistics based on
the supplied function.

## Usage

``` r
IPD_stat_factory(ipd_fun)
```

## Arguments

- ipd_fun:

  A function that computes mean and variance statistics for
  individual-level patient data.

## Value

A function that computes mean and variance statistics for a given
strategy.
