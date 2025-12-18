
#' Calculate simulated treatment comparison statistics
#'
#' @param strategy An object of class `strategy` created by functions such as
#'   [strategy_maic()], [strategy_stc()], or [strategy_mim()].
#'   Contains modelling details like the formula and family.
#' @param analysis_params List of analysis parameters. Must contain `ipd`
#'   (individual patient data).
#' @param ... Additional arguments.
#'
#' @return A list containing:
#' * `means`: A list containing:
#'     * `A`: Mean for comparator treatment group "A".
#'     * `C`: Mean for reference treatment group "C".
#' * `model`: The fitted [stats::glm()] object.
#'
#' @importFrom stats glm coef
#' @keywords internal
calc_stc <- function(strategy, analysis_params, ...) {
  args <- list(...)
  
  n_boot <- if (!is.null(args$N)) args$N else 1000
  
  # single fit
  run_stc_once <- function(data) {
    # centre covariates within this specific bootstrap sample
    centre_vars <- get_eff_mod_names(strategy$formula)
    data[, centre_vars] <- scale(data[, centre_vars], scale = FALSE)
    
    fit <- glm(formula = strategy$formula,
               family = strategy$family,
               data = data)
    
    # extract means
    coef_fit <- coef(fit)
    treat_coef_name <- grep(pattern = paste0("^", strategy$trt_var, "[^:]*$"),
                            names(coef_fit), value = TRUE)
    
    linkinv <- strategy$family$linkinv
    
    list(
      A = as.numeric(linkinv(coef_fit[1] + coef_fit[treat_coef_name])),
      C = as.numeric(linkinv(coef_fit[1])),
      fit = fit
    )
  }
  
  main_res <- run_stc_once(analysis_params$ipd)
  
  boot_results <- replicate(n_boot, simplify = FALSE, {
    # resample IPD with replacement
    boot_idx <- sample(nrow(analysis_params$ipd), replace = TRUE)
    boot_data <- analysis_params$ipd[boot_idx, ]
    
    run_stc_once(boot_data)
  })
  
  # collate bootstrap samples into vectors
  boot_A <- sapply(boot_results, function(x) x$A)
  boot_C <- sapply(boot_results, function(x) x$C)
  
  list(
    means = list(
      A = boot_A,
      C = boot_C),
    point_estimates = list(
      A = main_res$A, 
      C = main_res$C),  ##TODO: should we use these instead of `means`?
    model = list(
      fit = main_res$fit,
      N = n_boot)
  )
}
