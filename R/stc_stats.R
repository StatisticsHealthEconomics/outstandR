
#'
stc_stats <- function(formula =
                        as.formula("y ~ X3 + X4 +
                                   trt*I(X1-BC.ALD$mean.X1) +
                                   trt*I(X2-BC.ALD$mean.X2)"),
                      data = AC.IPD) {
  fit_stc <-
    glm(formula, data = data,
        family = binomial)
  
  # fitted treatment coefficient is relative A vs C conditional effect
  list(mean = coef(fit_stc)["trt"],
       var = vcov(fit_stc)["trt", "trt"])
}
