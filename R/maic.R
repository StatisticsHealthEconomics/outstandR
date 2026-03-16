
#' Estimate MAIC weights 
#' 
#' Matching-adjusted indirect comparison weights.
#' Method is taken from
#' \insertCite{Signorovitch2010}{outstandR}.
#' 
#' @param X_EM Centred \eqn{S=1} effect modifiers IPD covariates; matrix or data frame
#' @return Estimated weights for each individual; vector
#' 
#' @importFrom stats optim
#' @references
#' \insertRef{Signorovitch2010}{outstandR}
#' @keywords internal
#' 
maic_weights <- function(X_EM) {
  X_EM <- as.matrix(X_EM)
  
  N <- nrow(X_EM)    # number of individuals
  K <- ncol(X_EM)    # number of covariates
  
  init <- rep(1, K)  # arbitrary starting point for optimizer
  
  ##TODO: what about scaling X_EM?
  ##      because large values return error
  
  # find betas
  Q.min <- optim(fn = Q, X = X_EM, par = init, method = "BFGS")
  
  # check for convergence issues
  if (Q.min$convergence != 0) {
    warning(
      paste("optim did not converge (code:", Q.min$convergence, "). Message:", Q.min$message))
  }
  
  # finite solution is the logistic regression parameters
  hat_beta <- Q.min$par
  
  # Calculate the log weights
  # log(w_i) = sum(beta_k * X_ik)
  log.hat_w <- as.vector(X_EM %*% hat_beta)
  
  exp(log.hat_w)
}

#' Objective function to minimize for standard method of moments MAIC
#'
#' @param beta Beta coefficient to find
#' @param X Covariate value matrix, centred
#' @return Numeric value
#' @keywords internal
#' 
Q <- function(beta, X) {
  sum(exp(X %*% beta))
}

#' MAIC bootstrap sample
#' 
#' Matching-adjusted indirect comparison bootstrap sampling.
#' 
#' @eval study_data_args(include_ipd = TRUE, include_ald = TRUE)
#' @param indices Vector of indices, same length as original,
#'   which define the bootstrap sample
#' @eval reg_args(include_formula = TRUE, include_family = TRUE)
#' @param hat_w MAIC weights; default `NULL` which calls [maic_weights()]
#' 
#' @return Vector of fitted probabilities for treatments _A_ and _C_
#' @importFrom glue glue
#' @seealso [calc_IPD_stats.maic()]
#' 
#' @keywords internal
#' 
maic.boot <- function(ipd, indices = 1:nrow(ipd),
                      outcome_model, balance_model,
                      family, ald,
                      trt_var,
                      moments = 1,
                      int = FALSE,
                      hat_w = NULL) {
  
  dat <- ipd[indices, ]  # bootstrap sample
  n_ipd <- length(indices)
  n_trts <- length(unique(dat[[trt_var]]))
  
  # ensure bootstrap sample contains more than one treatment level
  if (n_trts < 2) {
    warning("Bootstrap sample contains less than two treatment levels. Returning NA.")
    return(c(pC = NA, pA = NA, rep(NA, n_ipd), ESS = NA))
  }
  
  base_vars <- all.vars(balance_model)
  
  if (length(base_vars) > 0) {
    term_list <- base_vars
    
    # Add interactions if requested
    if (isTRUE(int) && length(base_vars) > 1) {
      term_list <- c(term_list, paste0("(", paste(base_vars, collapse = " + "), ")^2"))
    }
    
    # Add squared terms if moments == 2
    if (moments == 2) {
      term_list <- c(term_list, paste0("I(", base_vars, "^2)"))
    }
    
    # Create the expanded formula and generate the IPD matrix
    expanded_formula <- as.formula(paste("~", paste(term_list, collapse = " + ")))
    
    # model.matrix computes the squares and interactions directly from the data
    ipd_matrix <- model.matrix(expanded_formula, data = dat)[, -1, drop = FALSE]
    balance_var_names <- colnames(ipd_matrix)
    
    X_EM_prepared <- matrix(NA, nrow = n_ipd, ncol = length(balance_var_names))
    colnames(X_EM_prepared) <- balance_var_names
    
    # determine covariate types and centre
    for (em_name in balance_var_names) {
      
      ipd_col <- ipd_matrix[, em_name]
      
      is_second_moment <- grepl("^I\\(.*\\^2\\)$", em_name)
      
      # catch squared terms and calculate target using SD ---
      
      if (is_second_moment) {
        
        # Extract the original variable name (e.g., "PF_cont_1" from "I(PF_cont_1^2)")
        base_name <- sub("^I\\((.*)\\^2\\)$", "\\1", em_name)
        
        # Get mean and sd of the base variable from the ALD
        base_mean <- ald |>
          dplyr::filter(.data$variable == base_name, .data$statistic == "mean") |>
          dplyr::pull(.data$value)
        
        base_sd <- ald |>
          dplyr::filter(.data$variable == base_name, .data$statistic == "sd") |>
          dplyr::pull(.data$value)
        
        if (length(base_mean) == 0 || length(base_sd) == 0) {
          stop(paste("Both 'mean' and 'sd' must be in the ALD to balance the variance of:", base_name), call. = FALSE)
        }
        
        # Calculate the ALD target: E[X^2] = SD^2 + Mean^2
        ald_squared_target <- (base_sd^2) + (base_mean^2)
        
        # Centre the IPD squared term
        X_EM_prepared[, em_name] <- ipd_col - ald_squared_target
        
        # Scale to improve optimizer performance
        col_sd <- sd(X_EM_prepared[, em_name])
        
        if (col_sd > 0) {
          X_EM_prepared[, em_name] <- X_EM_prepared[, em_name] / col_sd
        }
        
        next # Move to the next term in the loop
      }
      
      # standard logic for base variables and interactions ---
      
      ald_mean_val <- ald |>
        dplyr::filter(.data$variable == em_name, 
                      .data$statistic == "mean") |>
        dplyr::pull(.data$value)
      
      ald_prop_val <- ald |>
        dplyr::filter(.data$variable == em_name, 
                      .data$statistic == "prop") |>
        dplyr::pull(.data$value)
      
      if (length(ald_mean_val) > 0) {
        # continuous if 'mean' in ALD
        X_EM_prepared[, em_name] <- ipd_col - ald_mean_val
        
        col_sd <- sd(X_EM_prepared[, em_name])
        if (col_sd > 0) {
          X_EM_prepared[, em_name] <- X_EM_prepared[, em_name] / col_sd
        }
      } else if (length(ald_prop_val) > 0) {
        # binary if 'prop' in ALD
        X_EM_prepared[, em_name] <- ipd_col - ald_prop_val
      } else {
        stop(paste("Target statistic not found in ALD for covariate:", em_name), call. = FALSE)
      }
    }
    
    # calculate MAIC weights if not provided
    if (is.null(hat_w)) {
      hat_w <- maic_weights(X_EM = X_EM_prepared)
    }
  } else {
    hat_w <- rep(1, n_ipd)
  }
  
  # Calculate Effective Sample Size (ESS)
  ESS <- sum(hat_w)^2 / sum(hat_w^2)
  
  # so can use non-integer weights
  if (family$family == "binomial") {
    family <- quasibinomial()
  }
  
  # fit weighted regression model
  fit <- glm(formula = outcome_model,
             family = family,
             weights = hat_w / mean(hat_w),
             data = cbind(dat, hat_w = hat_w))
  
  # probabilities using inverse link
  linkinv <- family$linkinv
  
  # extract model coefficients
  coef_fit <- coef(fit)
  
  # index of the treatment variable
  trt_coef_name <- grep(trt_var, names(coef_fit), value = TRUE)
  
  # Handle potential failure to find exact match (e.g. factors)
  if (length(trt_coef_name) == 0) {
    stop("Could not find treatment coefficient in fitted model.", call. = FALSE)
  }
  
  # extract specific coefficients
  intercept <- coef_fit["(Intercept)"]
  trt_effect <- coef_fit[trt_coef_name]
  
  # calculate probabilities (assuming linear additivity on link scale)
  pC <- unname(linkinv(intercept))
  pA <- unname(linkinv(intercept + trt_effect))
  
  c(pC = pC, 
    pA = pA, 
    weights = unname(hat_w), 
    ESS = ESS)
}


#' Calculate MAIC
#' 
#' @param strategy An object of class `strategy` created by functions such as 
#'   [strategy_maic()], [strategy_stc()], or [strategy_mim()]. 
#'   Contains modelling details like the formula and family.
#' @param analysis_params List of analysis parameters. Must contain `ipd`
#'   (individual patient data) and `ald` (aggregated lead data).
#'
#' @return A list containing:
#' * `means`: A list containing:
#'     * `A`: Bootstrap estimates for comparator treatment group "A".
#'     * `C`: Bootstrap estimates for reference treatment group "C".
#' * `model`: A list containing model diagnostics derived from the original data:
#'     * `weights`: Vector of calculated weights for the patients in `ipd`.
#'     * `ESS`: The Effective Sample Size.
#'
#' @keywords internal
#' @importFrom boot boot
#' 
calc_maic <- function(strategy,
                      analysis_params) {
  
  verbose <- isTRUE(analysis_params$verbose)
  
  if (verbose) {
    cli::cli_h2("MAIC Execution")
    cli::cli_alert_info("Calculating weights using method of moments...")
  }
  
  calc_moments <- if (!is.null(strategy$moments)) strategy$moments else 1
  calc_int <- if (!is.null(strategy$int)) strategy$int else FALSE
  
  args_list <- 
    list(R = strategy$n_boot,
         balance_model = strategy$balance_model,
         outcome_model = strategy$outcome_model,
         family = strategy$family,
         trt_var = strategy$trt_var,
         data = analysis_params$ipd,
         ald = analysis_params$ald,
         moments = calc_moments,
         int = calc_int)
  
  if (verbose) {
    cli::cli_alert_info("Starting Bootstrap with {.val {strategy$n_boot}} replicates.")
    if (strategy$n_boot > 1000) {
      cli::cli_alert_warning("High iteration count detected. This may take some time.")
    }
  }
  
  maic_boot <- do.call(boot::boot, c(statistic = maic.boot, args_list))
  
  # boot return vector is: [pC (1), pA (2), weights (3 to N+2), ESS (N+3)]
  N_ipd <- nrow(analysis_params$ipd)
  
  idx_pC <- 1
  idx_pA <- 2
  idx_weights_start <- 3
  idx_weights_end <- 2 + N_ipd
  idx_ESS <- 2 + N_ipd + 1
  
  # 3. Extract Results
  # For MEANS: Use 't' (bootstrap replicates) to estimate uncertainty
  means_A_boot <- maic_boot$t[, idx_pA]
  means_C_boot <- maic_boot$t[, idx_pC]
  
  # For MODEL DIAGNOSTICS (Weights & ESS): Use original data
  # This gives the weights for the actual patients in the IPD
  # not shuffled bootstrap rows
  weights_orig <- maic_boot$t0[idx_weights_start:idx_weights_end]
  ESS_orig     <- maic_boot$t0[idx_ESS]
  
  list(
    means = list(
      A = means_A_boot,
      C = means_C_boot),
    model = list(
      weights = weights_orig,
      ESS = ESS_orig)
  )
}
