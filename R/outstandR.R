
#' @title Calculate the difference between treatments using all evidence
#' 
#' @description
#' This is the main, top-level wrapper for `{outstandR}`.
#' Methods taken from
#' \insertCite{RemiroAzocar2022}{outstandR}.
#' 
#' @param ipd_trial Individual-level patient data. For example, suppose between studies _A_ and _C_.
#'   In a long format and must contain a treatment column and outcome column consistent with the formula object.
#'   The labels in the treatment are used internally so there must be a common treatment with the aggregate-level data trial.
#' @param ald_trial Aggregate-level data. For example, suppose between studies _B_ and _C_. The column names are
#'  - `variable`: Covariate name. In the case of treatment arm sample size this is `NA`,
#'  - `statistic`: Summary statistic name from "mean", standard deviation "sd", probability "prop", or "sum",
#'  - `value`: Numerical value of summary statistic,
#'  - `trt`: Treatment label. Because we assume a common covariate distribution between treatment arms this is `NA`.
#' @param strategy Computation strategy function. These can be
#'    `strategy_maic()`, `strategy_stc()`, `strategy_gcomp_ml()` and `strategy_gcomp_bayes()`.
#' @param ref_trt Reference / common / anchoring treatment name.
#' @param CI Confidence interval level; between 0,1 with default 0.95.
#' @param scale Relative treatment effect scale. If `NULL`, the scale is automatically determined from the model.
#'   Choose from "log-odds", "log_relative_risk", "risk_difference", "delta_z", "mean_difference", "rate_difference" depending on the data type.
#' @param var_method Variance estimation method.
#' @param seed Random seed.
#' @param verbose Logical. If `TRUE`, prints progress messages and warnings.
#' @param ... Additional arguments. Currently, can pass named arguments to `rstanarm::stan_glm()` via `strategy_gcomp_bayes()`.
#' 
#' @return List of length 11 of statistics as a `outstandR` class object.
#'   Containing statistics between each pair of treatments.
#'   These are the mean, variances and confidence intervals,
#'   for contrasts and absolute values.
#'   
#' @importFrom Rdpack reprompt
#' @importFrom stats update
#' @importFrom withr local_seed
#' 
#' @seealso [strategy_maic()] [strategy_stc()] [strategy_gcomp_ml()] [strategy_gcomp_bayes()]
#' 
#' @references
#' \insertRef{RemiroAzocar2022}{outstandR}
#' 
#' @export
#' @examples
#' data(AC_IPD_binY_contX)  # A vs C individual patient-level data
#' data(BC_ALD_binY_contX)  # B vs C aggregate-level data
#' 
#' # linear formula
#' lin_form <- as.formula("y ~ PF_cont_1 + PF_cont_2 + trt*EM_cont_1 + trt*EM_cont_2")
#'                                 
#' # sampling values of additional arguments picked for speed
#' # select appropriate to specific analysis
#' 
#' # matching-adjusted indirect comparison
#' outstandR_maic <- outstandR(
#'   AC_IPD_binY_contX, BC_ALD_binY_contX,
#'   strategy = strategy_maic(formula = lin_form, n_boot = 100))
#' 
#' # simulated treatment comparison
#' outstandR_stc <- outstandR(
#'   AC_IPD_binY_contX, BC_ALD_binY_contX,
#'   strategy = strategy_stc(lin_form))
#' 
#' \donttest{
#' # G-computation with maximum likelihood
#' outstandR_gcomp_ml <- outstandR(
#'   AC_IPD_binY_contX, BC_ALD_binY_contX,
#'   strategy = strategy_gcomp_ml(lin_form, n_boot = 100, N =100))
#' 
#' # G-computation with Bayesian inference
#' outstandR_gcomp_bayes <- outstandR(
#'   AC_IPD_binY_contX, BC_ALD_binY_contX,
#'   strategy = strategy_gcomp_bayes(lin_form),
#'   chains = 1, iter = 1000, warmup = 20)
#' 
#' # Multiple imputation marginalization
#' outstandR_mim <- outstandR(
#'   AC_IPD_binY_contX, BC_ALD_binY_contX,
#'   strategy = strategy_mim(lin_form,
#'                           N = 100), # size of pseudo-population
#'   chains = 1, iter = 1000, warmup = 20)
#' }
#' 
outstandR <- function(ipd_trial, ald_trial, strategy,
                      ref_trt = NA,
                      CI = 0.95, 
                      scale = NULL, 
                      var_method = NULL,
                      seed = NULL,
                      verbose = TRUE,
                      ...) {
  if (!is.null(seed)) {
    withr::local_seed(seed)
  }
  
  cl <- match.call()
  
  validate_outstandr(ipd_trial, ald_trial, strategy, CI, scale)

  if (verbose) {
    cli::cli_h1("Starting outstandR Analysis")
    cli::cli_alert_info("Strategy: {.strong {class(strategy)[1]}}")
  }
  
  trt_var <- strategy$trt_var
  
  # evaluate the balance model if it was delayed as a function
  if (is.function(strategy$balance_model)) {
    strategy$balance_model <- strategy$balance_model(ald_trial)
    
    if (!is.null(strategy$balance_model)) {
      # Log it for the user and run the validation we skipped earlier
      cli::cli_alert_info("Auto-generated balance model: {.var {deparse(strategy$balance_model)}}")
      check_balance_formula(strategy$balance_model, trt_var)
    }
  }
  
  # combine outcome_model and balance_model into a single formula
  combined_formula <- strategy$outcome_model
  
  if (!is.null(strategy$balance_model)) {
    balance_terms <- attr(terms(strategy$balance_model), "term.labels")
    
    if (length(balance_terms) > 0) {
      # add balance terms to the right-hand side of outcome formula
      combined_formula <- stats::update(
        combined_formula, 
        paste(". ~ . +", paste(balance_terms, collapse = " + "))
      )
    }
  }
  
  ipd <- prep_ipd(combined_formula, ipd_trial)
  ald <- prep_ald(combined_formula, ald_trial, trt_var = trt_var)
  
  ref_trt <- get_ref_trt(ref_trt, trt_var, ipd_trial, ald_trial)
  
  # treatment names for each study
  ipd_comp <- get_comparator(ipd, ref_trt, trt_var)
  ald_comp <- get_comparator(ald, ref_trt, trt_var)
  
  if (is.null(scale)) scale <- get_treatment_effect(strategy$family$link)
  
  analysis_params <- list(
    ipd = ipd, 
    ald = ald,
    scale = scale,
    trt_var = trt_var,
    ref_trt = ref_trt,
    ipd_comp = ipd_comp,
    ald_comp = ald_comp,
    var_method = var_method,
    verbose = verbose
  )
  
  analysis_params <- add_seed(strategy, analysis_params, seed)
  
  ipd_stats <- calc_IPD_stats(strategy, analysis_params, ...) 
  ald_stats <- calc_ALD_stats(strategy, analysis_params) 
  
  stats <- result_stats(ipd_stats, ald_stats, CI)
  
  structure(
    .Data = list(
      results = stats,
      call = cl,
      outcome_model = strategy$outcome_model, 
      balance_model = strategy$balance_model,
      formula = combined_formula,
      CI = CI,
      ref_trt = ref_trt,
      ipd_comp = ipd_comp,
      ald_comp = ald_comp,
      scale = scale,
      var_method = var_method,
      family = strategy$family$family,
      model = c(method_name = ipd_stats$method_name,
                ipd_stats$model)),
    class = c("outstandR", class(stats)))
}

# only pass random seed to specific strategies ---

add_seed <- function(strategy, fn_args, seed) {
  UseMethod("add_seed")
}

#' @export
add_seed.default <- function(strategy, fn_args, seed) {
  return(fn_args)
}

#' @export
add_seed.gcomp_bayes <- function(strategy, fn_args, seed) {
  if (!is.null(seed)) {
    fn_args$seed <- seed
  }
  return(fn_args)
}
