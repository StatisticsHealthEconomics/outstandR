# outstandR (development version) 1.1.0

### New features

-   **Explicit Model Specifications**: The package architecture has been
    fundamentally updated to explicitly separate outcome prediction from
    population weighting. The `formula` argument across all
    `strategy_*()` functions now accepts a list containing
    `outcome_model` and `balance_model`. This conceptual split makes the
    statistical intent of methods like MAIC much clearer (i.e., defining
    exactly which covariates to balance vs. which to regress on) and
    future-proofs the package for upcoming doubly-robust methods that
    require both models simultaneously.

    ```r
    # New recommended approach 
    strategy_maic(
      formula = list(
        outcome_model = y ~ trt + age + sex, 
        balance_model = ~ age + sex
      )
    )

    # Legacy approach (still supported) 
    strategy_maic(formula = y ~ trt + age + sex)
    ```

-   To support this new architecture, `strategy_maic()` introduces
    intelligent defaults for missing formula components:

    -   If a traditional single formula is passed (legacy support), it
        automatically infers the `balance_model` by stripping the
        response and treatment variables.
    
    -   If `balance_model` is omitted from a list, it acts as a delayed
        function defaulting to a linear sum of all covariates found in
        the aggregate-level data (ALD).
    
    -   If `outcome_model` is omitted, it gracefully defaults to
        `y ~ trt`.

-   `outstandR()` gains a new `verbose` argument (defaulting to `TRUE`),
    allowing users to silence console output during simulations or
    repetitive tasks. It also introduces proactive warnings for
    computationally expensive operations (e.g., very high `n_boot` for
    MAIC, or large $N \times R$ for ML G-computation) to prevent users
    from thinking the session has hung.

### Minor improvements and fixes

-   `guess_treatment_name()` (used under the hood by all strategies when
    `trt_var` is `NULL`) is now smarter. It correctly identifies the
    treatment variable by checking for duplicated variables in
    interaction terms, or by falling back to the last main effect term.

-   `check_balance_formula()` has been added to provide strict
    validation for balance models. It ensures the formula is one-sided
    (no response variable) and throws an informative warning if the
    treatment variable is mistakenly included.

-   `check_formula()` now explicitly intercepts `Surv` objects in
    outcome models. It throws a clear, informative error stating that
    survival data support is officially scheduled for a later version.

### Deprecated and defunct

-   `strategy_stc()` is now deprecated and will be removed in a future
    release. Simulated Treatment Comparison (STC) is being phased out in
    favor of G-computation approaches (e.g., `strategy_gcomp_ml()`),
    which offer better statistical properties for population adjustment.

# outstandR 1.0.0

`outstandR` 1.0.0 marks a major milestone, transforming the original
code from the Remiro-Azocar study into a generalized, user-friendly R
package.

### Breaking changes & major updates

-   The aggregate-level data (ALD) input has been redesigned to use a
    standard long format, replacing the legacy format used in the
    original study.

-   Treatment labels and covariate names are no longer hard-coded.
    `outstandR` now dynamically extracts treatment arms (e.g., rather
    than "A", "B", "C") and general covariates directly from your input
    data (5119596).

-   Covariate generation functions have been separated and moved into
    their own standalone package,
    [`simcovariates`](https://www.google.com/search?q=%5Bhttps://github.com/n8thangreen/simcovariates%5D(https://github.com/n8thangreen/simcovariates)).

### Expanded modelling capabilities

-   The `family` argument in `outstandR()` has been expanded to support
    binary (`binomial`), continuous (`gaussian`), and count data
    (`poisson`).

-   Input data can now contain any combination of binary and continuous
    covariates, specifically upgraded for MAIC (4fb1e21).

-   `simulate_ALD_pseudo_pop()` now accepts user-provided marginal
    arguments for a target distribution. These can be defined optionally
    via `marginal_distns` and `marginal_params` in `strategy_gcomp_ml()`
    and `strategy_gcomp_bayes()`.

-   Users can now explicitly select the outcome scale (e.g., log-odds,
    risk difference, relative risk) using the new `scale` argument in
    `outstandR()`.

-   Added support for optional correlation structures in ALD simulations
    (d48d0ab).

### Output & documentation improvements

-   The main `outstandR()` function now returns absolute values in
    addition to treatment contrasts (2f1fbd7).

-   Added a dedicated `print()` method for `outstandR` objects to
    cleanly display analysis results (983cb2f).

-   Comprehensive, separate vignettes demonstrating workflows for
    binary, continuous, and count outcome data (fc10a68).

-   A complete package documentation website has been built using
    `pkgdown` and is now available via GitHub Pages (4df16f2).
