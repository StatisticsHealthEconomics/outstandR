
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
#' @param data Individual-level patient data (data frame).
#' @param indices Vector of indices, same length as original,
#'   which define the bootstrap sample.
#' @param balance_matrix Pre-computed balance matrix.
#' @param outcome_x_matrix Pre-computed outcome design matrix.
#' @param outcome_y Pre-computed outcome vector.
#' @param ald_targets Vector of ALD targets.
#' @param scaling_factors Vector of scaling factors.
#' @param trt_var Treatment variable name.
#' @param family A 'family' object specifying the distribution and link function.
#' @param hat_w MAIC weights; default `NULL` which calls [maic_weights()].
#' @param ipd Backwards compatibility IPD data (optional).
#' @param outcome_model Backwards compatibility outcome model formula (optional).
#' @param balance_model Backwards compatibility balance model formula (optional).
#' @param ald Backwards compatibility ALD data (optional).
#' @param moments Backwards compatibility moments (default 1).
#' @param int Backwards compatibility interactions flag (default FALSE).
#' 
#' @return Vector of fitted probabilities for treatments _A_ and _C_
#' @importFrom glue glue
#' @importFrom stats glm.fit
#' @seealso [calc_IPD_stats.maic()]
#' 
#' @keywords internal
#' 
maic.boot <- function(data, indices, 
                      balance_matrix, outcome_x_matrix, outcome_y, 
                      ald_targets, scaling_factors, 
                      trt_var, family, hat_w = NULL,
                      # backwards compatibility arguments
                      ipd = NULL, outcome_model = NULL, balance_model = NULL,
                      ald = NULL, moments = 1, int = FALSE) {
  
  # Backwards compatibility check
  is_old_signature <- !is.null(ipd) || 
                      !is.null(outcome_model) || 
                      !is.null(balance_model) || 
                      !is.null(ald) || 
                      (inherits(balance_matrix, "formula") || inherits(balance_matrix, "call"))
  
  if (is_old_signature) {
    old_ipd <- if (!is.null(ipd)) ipd else if (missing(data)) NULL else data
    old_indices <- if (missing(indices)) NULL else indices
    if (is.null(old_indices) && !is.null(old_ipd)) {
      old_indices <- 1:nrow(old_ipd)
    }
    old_outcome_model <- if (!is.null(outcome_model)) outcome_model else if (missing(balance_matrix)) NULL else balance_matrix
    old_balance_model <- if (!is.null(balance_model)) balance_model else if (missing(outcome_x_matrix)) NULL else outcome_x_matrix
    old_family <- if (!is.null(family)) family else if (missing(family)) {
      if (!missing(outcome_y) && inherits(outcome_y, "family")) outcome_y else NULL
    } else family
    old_ald <- if (!is.null(ald)) ald else if (missing(ald_targets)) NULL else ald_targets
    old_trt_var <- if (!is.null(trt_var)) trt_var else if (missing(trt_var)) {
      if (!missing(scaling_factors) && is.character(scaling_factors)) scaling_factors else NULL
    } else trt_var
    old_moments <- if (!is.null(moments)) moments else 1
    old_int <- if (!is.null(int)) int else FALSE
    
    data <- old_ipd
    indices <- old_indices
    family <- old_family
    trt_var <- old_trt_var
    
    calc_moments <- if (!is.null(old_moments)) old_moments else 1
    calc_int <- if (!is.null(old_int)) old_int else FALSE
    
    base_vars <- all.vars(old_balance_model)
    term_list <- base_vars
    
    if (isTRUE(calc_int) && length(base_vars) > 1) {
      term_list <- c(term_list, paste0("(", paste(base_vars, collapse = " + "), ")^2"))
    }
    if (calc_moments == 2) {
      term_list <- c(term_list, paste0("I(", base_vars, "^2)"))
    }
    
    if (length(term_list) == 0) {
      balance_matrix <- matrix(nrow = nrow(data), ncol = 0)
    } else {
      expanded_formula <- as.formula(paste("~", paste(term_list, collapse = " + ")))
      balance_matrix <- model.matrix(expanded_formula, data = data)[, -1, drop = FALSE]
    }
    
    outcome_x_matrix <- model.matrix(old_outcome_model, data = data)
    outcome_y <- data[[all.vars(old_outcome_model)[1]]]
    
    balance_var_names <- colnames(balance_matrix)
    ald_targets <- setNames(numeric(length(balance_var_names)), balance_var_names)
    scaling_factors <- setNames(rep(1, length(balance_var_names)), balance_var_names)
    
    for (em_name in balance_var_names) {
      is_second_moment <- grepl("^I\\(.*\\^2\\)$", em_name)
      if (is_second_moment) {
        base_name <- sub("^I\\((.*)\\^2\\)$", "\\1", em_name)
        base_mean <- old_ald$value[old_ald$variable == base_name & old_ald$statistic == "mean"]
        base_sd   <- old_ald$value[old_ald$variable == base_name & old_ald$statistic == "sd"]
        
        if (length(base_mean) == 0 || length(base_sd) == 0) {
          stop(paste("Both 'mean' and 'sd' must be in the ALD to balance the variance of:", base_name), call. = FALSE)
        }
        ald_targets[em_name] <- (base_sd^2) + (base_mean^2)
        centred_col <- balance_matrix[, em_name] - ald_targets[em_name]
        if (sd(centred_col) > 0) scaling_factors[em_name] <- sd(centred_col)
      } else {
        ald_mean_val <- old_ald$value[old_ald$variable == em_name & old_ald$statistic == "mean"]
        ald_prop_val <- old_ald$value[old_ald$variable == em_name & old_ald$statistic == "prop"]
        
        if (length(ald_mean_val) > 0) {
          ald_targets[em_name] <- ald_mean_val
          centred_col <- balance_matrix[, em_name] - ald_targets[em_name]
          if (sd(centred_col) > 0) scaling_factors[em_name] <- sd(centred_col)
        } else if (length(ald_prop_val) > 0) {
          ald_targets[em_name] <- ald_prop_val
        } else {
          stop(paste("Target statistic not found in ALD for covariate:", em_name), call. = FALSE)
        }
      }
    }
  }
  
  n_ipd <- length(indices)
  n_trts <- length(unique(data[[trt_var]][indices]))
  
  # ensure bootstrap sample contains more than one treatment level
  if (n_trts < 2) {
    warning("Bootstrap sample contains less than two treatment levels. Returning NA.")
    return(c(pC = NA, pA = NA, rep(NA, n_ipd), ESS = NA))
  }
  
  if (ncol(balance_matrix) > 0) {
    # 1. Subset the pre-made matrix
    X_EM_boot <- balance_matrix[indices, , drop = FALSE]
    
    # Centering and scaling in one step
    X_EM_prepared <- scale(X_EM_boot, center = ald_targets, scale = scaling_factors)
    
    # clean matrix
    attr(X_EM_prepared, "scaled:center") <- NULL
    attr(X_EM_prepared, "scaled:scale") <- NULL
    
    # 4. Find weights
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
  fit <- stats::glm.fit(
    x = outcome_x_matrix[indices, , drop = FALSE], 
    y = outcome_y[indices], 
    weights = hat_w / mean(hat_w), 
    family = family)
  
  
  # extract model coefficients
  coef_fit <- fit$coefficients  #coef_fit <- coef(fit)
  
  # index of the treatment variable
  trt_coef_name <- grep(trt_var, names(coef_fit), value = TRUE)
  
  # Handle potential failure to find exact match (e.g. factors)
  if (length(trt_coef_name) == 0) {
    stop("Could not find treatment coefficient in fitted model.", call. = FALSE)
  }
  
  # extract specific coefficients
  intercept <- coef_fit["(Intercept)"]
  trt_effect <- coef_fit[trt_coef_name]
  
  # probabilities using inverse link
  linkinv <- family$linkinv
  
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
  
  # extract parameters
  ipd <- analysis_params$ipd
  ald <- analysis_params$ald
  
  n_trts <- length(unique(ipd[[strategy$trt_var]]))
  if (n_trts < 2) {
    warning("Bootstrap sample contains less than two treatment levels. Returning NA.", call. = FALSE)
    warning("Bootstrap sample contains less than two treatment levels. Returning NA.", call. = FALSE)
    return(list(
      means = list(
        A = rep(NA, strategy$n_boot),
        C = rep(NA, strategy$n_boot)
      ),
      model = list(
        weights = rep(NA, nrow(ipd)),
        ESS = NA
      )
    ))
  }
  
  calc_moments <- if (!is.null(strategy$moments)) strategy$moments else 1
  calc_int <- if (!is.null(strategy$int)) strategy$int else FALSE
  
  # --- 1. BUILD THE BALANCE MATRIX ONCE ---
  
  base_vars <- all.vars(strategy$balance_model)
  term_list <- base_vars
  
  # Add interactions if requested
  if (isTRUE(calc_int) && length(base_vars) > 1) {
    term_list <- c(term_list, paste0("(", paste(base_vars, collapse = " + "), ")^2"))
  }
  
  # Add squared terms if moments == 2
  if (calc_moments == 2) {
    term_list <- c(term_list, paste0("I(", base_vars, "^2)"))
  }
  
  # Create the expanded formula and generate the IPD matrix
  if (length(term_list) == 0) {
    balance_matrix <- matrix(nrow = nrow(ipd), ncol = 0)
  } else {
    expanded_formula <- as.formula(paste("~", paste(term_list, collapse = " + ")))
    balance_matrix <- model.matrix(expanded_formula, data = ipd)[, -1, drop = FALSE]
  }
  
  # --- 2. BUILD THE OUTCOME MATRICES ONCE ---
  
  # For glm.fit, we need the numeric X matrix and Y vector
  outcome_x_matrix <- model.matrix(strategy$outcome_model, data = ipd)
  # Assumes the first variable in the formula is the outcome (y)
  outcome_y <- ipd[[all.vars(strategy$outcome_model)[1]]]
  
  # --- 3. PRE-CALCULATE TARGETS AND SCALES ---
  
  balance_var_names <- colnames(balance_matrix)
  ald_targets <- setNames(numeric(length(balance_var_names)), balance_var_names)
  scaling_factors <- setNames(rep(1, length(balance_var_names)), balance_var_names)
  
  # 3. Do the expensive dplyr/string logic ONCE
  for (em_name in balance_var_names) {
    
    is_second_moment <- grepl("^I\\(.*\\^2\\)$", em_name)
    
    if (is_second_moment) {
      base_name <- sub("^I\\((.*)\\^2\\)$", "\\1", em_name)
      
      base_mean <- ald$value[ald$variable == base_name & ald$statistic == "mean"]
      base_sd   <- ald$value[ald$variable == base_name & ald$statistic == "sd"]
      
      if (length(base_mean) == 0 || length(base_sd) == 0) {
        stop(paste("Both 'mean' and 'sd' must be in the ALD to balance the variance of:", base_name), call. = FALSE)
      }
      
      ald_targets[em_name] <- (base_sd^2) + (base_mean^2)
      
      # Pre-calculate standard deviation for scaling based on original data
      centred_col <- balance_matrix[, em_name] - ald_targets[em_name]
      if (sd(centred_col) > 0) scaling_factors[em_name] <- sd(centred_col)
      
    } else {
      # Standard logic
      ald_mean_val <- ald$value[ald$variable == em_name & ald$statistic == "mean"]
      ald_prop_val <- ald$value[ald$variable == em_name & ald$statistic == "prop"]
      
      if (length(ald_mean_val) > 0) {
        ald_targets[em_name] <- ald_mean_val
        centred_col <- balance_matrix[, em_name] - ald_targets[em_name]
        
        if (sd(centred_col) > 0) scaling_factors[em_name] <- sd(centred_col)
        
      } else if (length(ald_prop_val) > 0) {
        ald_targets[em_name] <- ald_prop_val
        # Proportions usually aren't scaled, so scale stays 1
      } else {
        stop(paste("Target statistic not found in ALD for covariate:", em_name), call. = FALSE)
      }
    }
  }
  
  # --- 4. RUN THE BOOTSTRAP ---
  
  args_list <- list(
    data = ipd,                    # Required by boot, even if we use matrices
    statistic = maic.boot,
    R = strategy$n_boot,
    # Pass our pre-computed matrices and vectors
    balance_matrix = balance_matrix,
    outcome_x_matrix = outcome_x_matrix,
    outcome_y = outcome_y,
    ald_targets = ald_targets,    
    scaling_factors = scaling_factors,
    family = strategy$family,
    trt_var = strategy$trt_var
  )
  
  if (verbose) {
    cli::cli_alert_info("Starting Bootstrap with {.val {strategy$n_boot}} replicates.")
  }
  
  maic_boot <- do.call(boot::boot, args_list)
  
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
