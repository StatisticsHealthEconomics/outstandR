
#' @export
#' @method print strategy
print.strategy <- function(x, ...) {
  # Map internal names to clean titles
  method_names <- c(
    maic = "MAIC",
    stc = "STC",
    gcomp_ml = "G-Computation (ML)",
    gcomp_bayes = "G-Computation (Bayesian)",
    mim = "MIM"
  )
  
  method_title <- method_names[class(x)[1]]
  
  cat(pillar::style_bold(paste0("outstandR Strategy Specification: ", method_title)), "\n")
  cat("  Formula: ", deparse(x$formula), "\n")
  cat("  Family:  ", x$family$family, paste0("(", x$family$link, " link)"), "\n")
  
  invisible(x)
}

#' @export
#' @method print maic
print.maic <- function(x, ...) {
  NextMethod()
  cat(pillar::style_subtle("  Parameters:"), "\n")
  cat("    Bootstrap Resamples (n_boot):", x$n_boot, "\n")
  cat("    Treatment Variable:     ", x$trt_var, "\n")
  invisible(x)
}

#' @export
#' @method print stc
print.stc <- function(x, ...) {
  NextMethod()
  cat(pillar::style_subtle("  Parameters:"), "\n")
  cat("    Treatment Variable:     ", x$trt_var, "\n")
  invisible(x)
}

#' @export
#' @method print gcomp_ml
print.gcomp_ml <- function(x, ...) {
  NextMethod()
  cat(pillar::style_subtle("  Parameters:"), "\n")
  cat("    Synthetic Sample Size (N):", x$N, "\n")
  cat("    Bootstrap Resamples (n_boot):  ", x$n_boot, "\n")
  if (!all(is.na(x$marginal_distns))) {
    cat("    Marginal Distns:          ", paste(x$marginal_distns, collapse = ", "), "\n")
  }
  invisible(x)
}

#' @export
#' @method print gcomp_bayes
print.gcomp_bayes <- function(x, ...) {
  NextMethod()
  cat(pillar::style_subtle("  Parameters:"), "\n")
  cat("    Synthetic Sample Size (N):", x$N, "\n")
  if (!all(is.na(x$rho))) {
    cat("    Correlation (rho):         Provided\n")
  }
  invisible(x)
}

#' @export
#' @method print mim
print.mim <- function(x, ...) {
  NextMethod()
  cat(pillar::style_subtle("  Parameters:"), "\n")
  cat("    Synthetic Sample Size (N):", x$N, "\n")
  cat("    Treatment Variable:       ", x$trt_var, "\n")
  invisible(x)
}
