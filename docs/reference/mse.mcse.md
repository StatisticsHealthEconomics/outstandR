# Monte Carlo SE of MSE estimate

Monte Carlo SE of MSE estimate

## Usage

``` r
mse.mcse(theta.hat, theta)
```

## Arguments

- theta.hat:

  Theta hat

- theta:

  Theta

## Value

\\sqrt(sum((tmp - mse.est)^2)/(nsim\*(nsim-1)))\\
