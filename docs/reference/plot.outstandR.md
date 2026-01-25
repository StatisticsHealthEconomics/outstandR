# Default Plot Method for outstandR Objects

Default Plot Method for outstandR Objects

## Usage

``` r
# S3 method for class 'outstandR'
plot(x, ..., type = c("both", "contrasts", "absolute"), labels = NULL)
```

## Arguments

- x:

  An object of class 'outstandR' or a list of 'outstandR' objects.

- ...:

  Additional 'outstandR' objects for comparison.

- type:

  Character, one of "both" (default), "contrasts", or "absolute".

- labels:

  Optional character vector of names for the models.

## Value

A
[`ggplot2::ggplot()`](https://ggplot2.tidyverse.org/reference/ggplot.html)
object representing the forest plot of the results.
