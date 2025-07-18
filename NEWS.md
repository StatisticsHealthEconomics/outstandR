# outstandR (development version) 0.009


* Original aggregate data format from Remiro-Azocar study changed to a long format

* Original covariate names from Remiro-Azocar study generalised and user-defined

* Input data can now contain any combination of binary and continuous covariates
  * For MAIC (4fb1e21)
  
* Treatment names taken from data, rather than hard coded as "A", "B", "C" (5119596)

* `simulate_ALD_pseudo_pop()` can now take user-provided marginal arguments for a target distribution with
    `marginal_distns`, `marginal_params`. These are defined optionally in
    `strategy_gcomp_ml()` and `strategy_gcomp_stan()`

* Included correlation of ALD simulations as optional function argument (d48d0ab)

* `print` method available (983cb2f)

* `outstandR()` call returns absolute values, as well as contrasts (2f1fbd7)

* `outstandR()` `family` argument now takes binary (`binomial`), continuous (`gaussian`) and count data (`poisson`)

* Separate vignettes written for binary, continuous and count outcome data (fc10a68)

* `pkgdown` site built and published on GitHub (4df16f2)

* `outstandR()` `scale` argument gives the ability to select outcome scale from log-odds, risk difference, relative risk

* Moved covariate generation in to separate [`simcovariates`](https://github.com/n8thangreen/simcovariates) package
