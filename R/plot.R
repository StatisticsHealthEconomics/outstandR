#' Default Plot Method for outstandR Objects
#' 
#' @param x An object of class 'outstandR' or a list of 'outstandR' objects.
#' @param ... Additional 'outstandR' objects for comparison.
#' @param type Character, one of "both" (default), "contrasts", or "absolute".
#' @param labels Optional character vector of names for the models.
#' @import dplyr ggplot2 purrr
#' @export
#' 
plot.outstandR <- function(x, ..., 
                           type = c("both", "contrasts", "absolute"),
                           labels = NULL) {
  type <- match.arg(type)
  
  # Combine x with any additional objects passed via ...
  models <- list(x, ...)
  
  # Ensure all objects are of class 'outstandR'
  if (!all(sapply(models, function(m) inherits(m, "outstandR")))) {
    stop("All objects passed to plot() must be of class 'outstandR'.")
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
  
  plot_df <- plot_df %>% filter(!is.na(.data$Estimate))
  
  # combined forest plot
  ggplot(
    plot_df, aes(x = .data$Estimate, y = .data$Treatments, color = .data$Model)) +
    geom_point(position = position_dodge(width = 0.5), size = 3) +
    geom_errorbarh(aes(xmin = .data$lower.0.95, xmax = .data$upper.0.95), 
                   position = position_dodge(width = 0.5), height = 0.2) +
    facet_wrap(~.data$Type, scales = "free") +
    geom_vline(data = filter(plot_df, Type == "Relative Contrasts"), 
               aes(xintercept = 0), linetype = "dashed", color = "gray50") +
    labs(title = "Population-Adjusted Indirect Comparison Results",
         x = "Estimate (95% CI)", y = NULL) +
    theme_minimal() +
    theme(legend.position = "bottom")
}