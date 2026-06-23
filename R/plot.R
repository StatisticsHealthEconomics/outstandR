#' Default Plot Method for outstandR Objects
#' 
#' @param x An object of class 'outstandR' or a list of 'outstandR' objects.
#' @param ... Additional 'outstandR' objects for comparison.
#' @param type Character, one of "both" (default), "contrasts", or "absolute".
#' @param labels Optional character vector of names for the models.
#' @param include_naive Logical. Should naive (unadjusted) estimates be included in the plot? Default is TRUE.
#' @import dplyr ggplot2 purrr
#' @export
#' @return A [ggplot2::ggplot()] object representing the forest plot of the results.
#'
plot.outstandR <- function(x, ...,
                           type = c("both", "contrasts", "absolute"),
                           labels = NULL,
                           include_naive = TRUE) {
  type <- match.arg(type)
  
  # Combine x with any additional objects passed via ...
  models <- list(x, ...)
  
  # Ensure all objects are of class 'outstandR'
  if (!all(sapply(models, function(m) inherits(m, "outstandR")))) {
    stop("All objects passed to plot() must be of class 'outstandR'.")
  }
  
  # Ensure all objects are on the same scale
  scales <- sapply(models, function(m) {
    if (is.null(m$scale)) "Unknown" else m$scale
  })
  if (length(unique(scales)) > 1) {
    stop("All objects passed to plot() must be on the same scale. Found scales: ",
         paste(unique(scales), collapse = ", "), ".")
  }
  
  # Internal helper to convert nested lists into tidy data frame
  tidy_outstandr <- function(res_list) {
    if (is.null(res_list) || all(is.na(unlist(res_list$means)))) return(NULL)
    
    data.frame(
      Treatments = names(res_list$means),
      Estimate   = unname(unlist(res_list$means)),
      lower.0.95 = sapply(res_list$CI, `[`, 1),
      upper.0.95 = sapply(res_list$CI, `[`, 2),
      stringsAsFactors = FALSE
    )
  }
  
  plot_df <- map_df(seq_along(models), function(i) {
    m <- models[[i]]
    m_name <- if(!is.null(labels)) labels[i] else m$model$method_name
    
    out_list <- list()

    if (type %in% c("both", "contrasts")) {
      cont <- tidy_outstandr(m$results$contrasts)
      if (!is.null(cont)) out_list$cont <- mutate(cont, Type = "Relative Contrasts")
    }
    
    # Process Absolute (Will be NULL for MIM)
    if (type %in% c("both", "absolute")) {
      abs_val <- tidy_outstandr(m$results$absolute)
      if (!is.null(abs_val)) out_list$abs <- mutate(abs_val, Type = "Absolute Estimates")
    }
    
    bind_rows(out_list) %>% mutate(Model = m_name)
  })
  
  # ---------------------------------------------------------------------------
  # Inject Naive Estimates
  # ---------------------------------------------------------------------------
  if (include_naive && !is.null(models[[1]]$naive)) {
    m_naive <- models[[1]]$naive
    naive_list <- list()
    
    # Check if naive is structured identically to standard results (nested)
    if ("contrasts" %in% names(m_naive) || "absolute" %in% names(m_naive)) {
      if (type %in% c("both", "contrasts") && !is.null(m_naive$contrasts)) {
        cont_n <- tidy_outstandr(m_naive$contrasts)
        if (!is.null(cont_n)) naive_list$cont <- mutate(cont_n, Type = "Relative Contrasts", Model = "Naive (Unadjusted)")
      }
      if (type %in% c("both", "absolute") && !is.null(m_naive$absolute)) {
        abs_n <- tidy_outstandr(m_naive$absolute)
        if (!is.null(abs_n)) naive_list$abs <- mutate(abs_n, Type = "Absolute Estimates", Model = "Naive (Unadjusted)")
      }
    }
    # Check if naive is structured from a custom helper function (flat list)
    else if ("Estimate" %in% names(m_naive)) {
      if (type %in% c("both", "contrasts")) {
        naive_list$cont <- data.frame(
          Treatments = m_naive$Treatments,
          Estimate   = m_naive$Estimate,
          lower.0.95 = if (!is.null(m_naive$lower.0.95)) m_naive$lower.0.95 else NA,
          upper.0.95 = if (!is.null(m_naive$upper.0.95)) m_naive$upper.0.95 else NA,
          Type       = "Relative Contrasts",
          Model      = if (!is.null(m_naive$Method)) m_naive$Method else "Naive (Unadjusted)",
          stringsAsFactors = FALSE
        )
      }
    }
    
    plot_df <- bind_rows(plot_df, bind_rows(naive_list))
  }
  
  plot_df <- plot_df %>% dplyr::filter(!is.na(.data$Estimate))
  
  # Order models as factors to ensure the "Naive" method drops to the bottom of the legend
  model_levels <- unique(plot_df$Model)
  naive_names <- grep("Naive", model_levels, value = TRUE, ignore.case = TRUE)
  if (length(naive_names) > 0) {
    model_levels <- c(setdiff(model_levels, naive_names), naive_names)
    plot_df$Model <- factor(plot_df$Model, levels = model_levels)
  }
  
  # Define custom colors so that any naive model is plotted in black
  n_models <- length(model_levels)
  is_naive <- grepl("Naive", model_levels, ignore.case = TRUE)
  n_non_naive <- sum(!is_naive)
  
  if (n_non_naive > 0) {
    non_naive_colors <- grDevices::hcl(
      h = seq(15, 375, length.out = n_non_naive + 1)[1:n_non_naive], 
      c = 100, 
      l = 65
    )
  } else {
    non_naive_colors <- character(0)
  }
  
  model_colors <- stats::setNames(rep("black", n_models), model_levels)
  model_colors[!is_naive] <- non_naive_colors

  scale_lbl <- if (!is.null(x$scale)) {
    scale_map <- list(
      log_odds = "Log-Odds Ratio",
      mean_difference = "Mean Difference",
      delta_z = "Delta Z",
      log_relative_risk_rare_events = "Log Relative Risk (Rare)",
      log_relative_risk = "Log Relative Risk",
      risk_difference = "Risk Difference"
    )
    if (x$scale %in% names(scale_map)) scale_map[[x$scale]] else x$scale
  } else {
    "Unknown"
  }

  # combined forest plot
  ggplot(
    plot_df, aes(x = .data$Estimate, y = .data$Treatments, color = .data$Model)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbar(aes(xmin = .data$lower.0.95, xmax = .data$upper.0.95), 
                  position = position_dodge(width = 0.5), width = 0.2, na.rm = TRUE) +
    scale_color_manual(values = model_colors) +
    facet_wrap(~.data$Type, scales = "free") +
    geom_vline(data = dplyr::filter(plot_df, .data$Type == "Relative Contrasts"), 
               aes(xintercept = 0), linetype = "dashed", color = "gray50") +
    labs(title = "Population-Adjusted Indirect Comparison Results",
         x = paste0("Estimate (95% CI) [Scale: ", scale_lbl, "]"), y = NULL) +
    theme_bw() +
    theme(legend.position = "bottom")
}