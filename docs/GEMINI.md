# outstandR: Outcome Regression Standardisation

`outstandR` is an R package for population-adjusted indirect treatment
comparison (ITC). It implements various methods to adjust for
cross-trial differences in effect modifiers when individual patient data
(IPD) is available for one trial and only aggregate-level data (ALD) is
available for another.

## Core Architecture

The package follows a “Strategy” pattern where users define an
adjustment method using `strategy_*` functions and pass it to the main
[`outstandR()`](https://StatisticsHealthEconomics.github.io/outstandR/reference/outstandR.md)
wrapper.

### Key Methods

- **MAIC (Matching-Adjusted Indirect Comparison):** Non-parametric
  likelihood reweighting.
- **G-computation (ML & Bayesian):** Model-based standardisation (STC is
  deprecated in favour of G-computation).
- **MIM (Multiple Imputation Marginalization):** Parametric
  G-computation using multiple imputation.

### Main Entry Points

- `outstandR(ipd_trial, ald_trial, strategy, ...)`: The primary
  user-facing function.
- `strategy_maic(formula, ...)`: Define a MAIC strategy.
- `strategy_gcomp_ml(formula, ...)`: Define a Maximum Likelihood
  G-computation strategy.
- `strategy_gcomp_bayes(formula, ...)`: Define a Bayesian G-computation
  strategy.
- `strategy_mim(formula, ...)`: Define a Multiple Imputation
  Marginalization strategy.

## Development Workflows

### Standard Commands

- **Documenting:**
  [`devtools::document()`](https://devtools.r-lib.org/reference/document.html)
  (updates `NAMESPACE` and `man/`).
- **Testing:**
  [`devtools::test()`](https://devtools.r-lib.org/reference/test.html)
  or
  [`testthat::test_local()`](https://testthat.r-lib.org/reference/test_package.html).
- **Checking:**
  [`devtools::check()`](https://devtools.r-lib.org/reference/check.html)
  (standard R package quality check).
- **Building Site:**
  [`pkgdown::build_site()`](https://pkgdown.r-lib.org/reference/build_site.html)
  (updates the documentation website in `docs/`).

### Project Structure

- `R/`: Source code.
  - `outstandR.R`: Main wrapper.
  - `strategy_.R`: Strategy class definitions and logic.
  - `maic.R`, `gcomp_ml.R`, `mim.R`: Method-specific implementations.
- `tests/testthat/`: Unit tests organized by feature/method.
- `vignettes/`: Detailed usage examples and methodology documentation.
- `man/`: Auto-generated documentation (do not edit manually).
- `data/`: Built-in datasets for examples and testing.

## Engineering Standards

- **Consistency:** Follow standard R package conventions.
- **Documentation:** Use `roxygen2` with Markdown support
  (`Roxygen: list(markdown = TRUE)` in `DESCRIPTION`).
- **Testing:** New features or bug fixes MUST include `testthat` tests.
- **S3 Classes:** Strategies and the main output use S3 classes
  (`outstandR`, `strategy`, `maic`, etc.).
- **User Feedback:** Use `cli` and `glue` for informative console
  output.
- **Dependencies:** Prefer `dplyr`, `tidyr`, and `rlang` for data
  manipulation and tidy evaluation.

## Data Formats

### IPD (Individual Patient Data)

A data frame in long format containing: - A treatment column (e.g.,
`trt`). - An outcome column (e.g., `y`). - Covariates used in the
adjustment.

### ALD (Aggregate-Level Data)

A data frame with specific columns: - `variable`: Covariate name. -
`statistic`: “mean”, “sd”, “prop”, “sum”, or “N”. - `value`: Numeric
value. - `trt`: Treatment label (or `NA` if common across arms).
