
#' Simulate Aggregate-Level Data Pseudo-Population
#'
#' Generates a synthetic cohort using a normal copula based on aggregate-level data.
#'
#' @eval reg_args(include_formula = TRUE, include_family = FALSE)
#' @eval study_data_args(include_ipd = TRUE, include_ald = TRUE)
#' @param rho A named square matrix of covariate correlations or single value; default NA takes from IPD.
#' @param N Sample size for the synthetic cohort. Default is 1000.
#' @param marginal_distns Marginal distributions names; vector default NA.
#'    Available distributions are given in stats::Distributions. See [copula::Mvdc()] for details
#' @param marginal_params Marginal distributions parameters;
#'    named list of lists, default NA. See [copula::Mvdc()] for details.
#' @param seed Random seed
#' @param verbose Default `FALSE`
#' 
#' @return A data frame representing the synthetic pseudo-population.
#'         It contains \code{N} rows (one for each simulated individual) and 
#'         columns for every covariate specified in \code{marginal_distns} of \code{formula}.
#' @importFrom copula normalCopula mvdc
#' 
#' @keywords internal
#' @examples
#' \dontrun{
#' 
#' ## Example 1: Simulating data with explicitly defined marginals and
#' ## a provided correlation matrix (rho)
#' 
#' my_formula <- y ~ x1 + x2*trt + trt
#' 
#' # Define marginal distributions (Normal for age, Lognormal for BMI)
#' dists <- c(x1 = "norm", x2 = "norm")
#' 
#' # Define parameters matching the chosen distributions
#' params <- list(
#'   x1 = list(mean = 60, sd = 8),
#'   x2 = list(mean = 3, sd = 0.1)
#' )
#' 
#' # Define a 2x2 correlation matrix
#' corr_mat <- matrix(c(1, 0.25,
#'                      0.25, 1), nrow = 2,
#'                    dimnames = list(c("x1", "x2"),
#'                                    c("x1", "x2")))
#' 
#' # Generate the synthetic cohort
#' sim_cohort_marginals_rho <- simulate_ALD_pseudo_pop(
#'   formula = my_formula,
#'   ald = mock_ald,
#'   trt_var = "trt",
#'   rho = corr_mat,
#'   N = 100,
#'   marginal_distns = dists,
#'   marginal_params = params
#' )
#' 
#' head(sim_cohort_marginals_rho)
#' 
#' ## Example 2: Estimating the correlation matrix (rho) from provided IPD
#' 
#' # Create some mock Individual Patient Data (IPD)
#' mock_ipd <- data.frame(
#'   x1 = rnorm(200, 60, 8),
#'   x2 = rlnorm(200, 3.3, 0.1)
#' )
#' 
#' mock_ald <- data.frame(
#'   variable = c("x1", "x1", "x2", "x2"),
#'   statistic = c("mean", "sd", "mean", "sd"),
#'   value = c(1000, 500, 0.7, 0.1),
#'   trt = NA
#' )
#' 
#' # Generate the synthetic cohort using IPD for the correlation structure
#' sim_cohort <- simulate_ALD_pseudo_pop(
#'   formula = my_formula,
#'   ipd = mock_ipd,
#'   ald = mock_ald,
#'   trt_var = "trt",
#'   rho = NA,
#'   N = 100
#' )
#' 
#' head(sim_cohort)
#' 
#' # Generate the synthetic cohort using IPD for the correlation structure
#' # with marginals
#' sim_cohort_marginals <- simulate_ALD_pseudo_pop(
#'   formula = my_formula,
#'   ipd = mock_ipd,
#'   ald = mock_ald,
#'   trt_var = "trt",
#'   N = 100,
#'   marginal_distns = dists,
#'   marginal_params = params
#' )
#' 
#' head(sim_cohort_marginals)
#' }
#' 
simulate_ALD_pseudo_pop <- function(formula,
                                    ipd = NULL, ald = NULL,
                                    trt_var,
                                    rho = NA,
                                    N = 1000,
                                    marginal_distns = NA,
                                    marginal_params = NA,
                                    seed = NULL,
                                    verbose = FALSE) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  
  if (verbose) {
    cli::cli_alert_warning(c(
      "Gaussian copula methods may not properly capture the joint distribution ",
      "for certain marginal combinations (e.g., highly skewed variables, ",
      "extreme proportions, or strict bounds). It is highly recommended to ",
      "check the summary statistics of the simulated pseudo-population."
    ))
  }
  
  resolved <- 
    prepare_covariate_distns(
      formula, ald, trt_var, 
      marginal_distns, marginal_params, 
      verbose)
  
  marginal_distns <- resolved$distns
  marginal_params <- resolved$params
  
  covariate_names <- names(marginal_distns)
  n_covariates <- length(covariate_names)
  
  # --- Standard Simulation Logic ---
  
  # handle no covariates (intercept-only or treatment-only models)
  if (n_covariates == 0) {
    return(data.frame(row.names = seq_len(length.out = N)))
  }
  
  # don't require copula for single covariate
  if (n_covariates == 1) {
    # dynamically call appropriate random number generator
    rng_fun <- get(paste0("r", marginal_distns[1]))
    sim_vals <- do.call(rng_fun, c(list(n=N), marginal_params[[1]]))
    
    x_star <- matrix(sim_vals, ncol = 1,
                     dimnames = list(NULL, covariate_names))
    # return(as.data.frame(setNames(list(sim_vals), covariate_names)))  ##TODO
    return(as.data.frame(x_star))
  }
  
  # prepare correlation matrix
  if (!is.matrix(rho)) {
    if (is.na(rho)) {
      
      if (is.null(ipd)) {
        stop("'rho' must be provided when 'ipd' is not available.", call. = FALSE)
      }
      
      rho <- cor(ipd[, covariate_names], use = "pairwise.complete.obs")
    } else {
      rho <- matrix(rho, n_covariates, n_covariates,
                    dimnames = list(covariate_names, covariate_names))
      diag(rho) <- 1
    }
  }
  
  rho <- rho[covariate_names, covariate_names]
  
  # define Copula and Multivariate Distribution
  cor_params <- rho[lower.tri(rho, diag = FALSE)]
  
  cop <- copula::normalCopula(param = cor_params,
                              dim = n_covariates,
                              dispstr = "un")
  
  mvd <- copula::mvdc(copula = cop,
                      margins = marginal_distns,
                      paramMargins = marginal_params)
  
  # simulate data
  x_star <- as.data.frame(copula::rMvdc(n = N, mvd))
  colnames(x_star) <- covariate_names
  
  ##TOOD: calculate mean and sd of x_star and compare with target ALD means
  
  return(x_star)
}
