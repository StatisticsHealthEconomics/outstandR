# outstandR (development version) 1.0.1

### New features and improvements

* Changed `formula` argument of `strategy`s so that they now take a list containing `outcome_model` and `balance_model`, 
  depending on the type of model. This makes the MAIC strategy more explicit and future-proofs for doubly-robust methods. ()

* Added a `verbose` argument to the main `outstandR()` function (defaulting to `TRUE`). 
  Users can now silence console output when running simulations or repetitive tasks.
  Added proactive warnings for computationally expensive operations to prevent users from thinking the session has hung. The package will now gently warn you if:
    * The number of bootstrap iterations (`n_boot`) for MAIC is very high.
    * The combination of bootstrap iterations and synthetic cohort size ($N \times R$) for Maximum Likelihood G-computation is exceptionally large.
  
### Deprecated and defunct

* `strategy_stc()` has been deprecated and will be removed in a future release. 
   Simulated Treatment Comparison (STC) is being phased out in favor of G-computation approaches (`strategy_gcomp_ml()`), 
   which offer better statistical properties for these types of adjustments. ()

# outstandR 1.0.0

`outstandR` 1.0.0 marks a major milestone, transforming the original code from the 
Remiro-Azocar study into a generalized, user-friendly R package.

### Breaking changes & major updates

* The aggregate-level data (ALD) input has been redesigned to use 
  a standard long format, replacing the legacy format used in the original study.
  
* Treatment labels and covariate names are no longer hard-coded.
  `outstandR` now dynamically extracts treatment arms (e.g., rather than "A", "B", "C") 
  and general covariates directly from your input data (5119596).

* Covariate generation functions have been separated and moved 
  into their own standalone package, [`simcovariates`](https://www.google.com/search?q=%5Bhttps://github.com/n8thangreen/simcovariates%5D(https://github.com/n8thangreen/simcovariates)).

### Expanded modelling capabilities

* The `family` argument in `outstandR()` has been expanded to 
  support binary (`binomial`), continuous (`gaussian`), and count data (`poisson`).

* Input data can now contain any combination of binary and 
  continuous covariates, specifically upgraded for MAIC (4fb1e21).

* `simulate_ALD_pseudo_pop()` now accepts 
  user-provided marginal arguments for a target distribution. These can be defined optionally 
  via `marginal_distns` and `marginal_params` in `strategy_gcomp_ml()` and 
  `strategy_gcomp_bayes()`.

* Users can now explicitly select the outcome scale (e.g., log-odds,
  risk difference, relative risk) using the new `scale` argument in `outstandR()`.

* Added support for optional correlation structures in ALD 
  simulations (d48d0ab).

### Output & documentation improvements

* The main `outstandR()` function now returns absolute values 
  in addition to treatment contrasts (2f1fbd7).

* Added a dedicated `print()` method for `outstandR` objects to cleanly
  display analysis results (983cb2f).

* Comprehensive, separate vignettes demonstrating workflows for
  binary, continuous, and count outcome data (fc10a68).

* A complete package documentation website has been built using `pkgdown` 
  and is now available via GitHub Pages (4df16f2).


