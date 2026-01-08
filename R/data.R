
#' Individual-level patient data for binary outcome, continuous covariates
#'
#' This data set contains simulated patient covariate and outcome values.
#'
#' @name AC_IPD_binY_contX
#' @docType data
#' @format ## `y ~ PF_cont_1 + PF_cont_2 + trt + trt:(EM_cont_1 + EM_cont_2)`
#' \describe{
#'    \item{id}{Numeric unique identifier}
#'    \item{PF_cont_1}{Numeric prognostic factor continuous covariate}
#'    \item{PF_cont_2}{Numeric prognostic factor continuous covariate}
#'    \item{EM_cont_1}{Numeric effect modifier continuous covariate}
#'    \item{EM_cont_2}{Numeric effect modifier continuous covariate}
#'    \item{trt}{Factor treatment identifier. Levels A, C}
#'    \item{y}{Integer binary outcome}
#'    \item{true_eta}{Numeric linear predictor}
#' }
#' @references Remiro‐Azocar A, Heath A, Baio G (2022)
#' @source Simulated data
#' @keywords datasets
#' @usage data(AC_IPD_binY_contX)
"AC_IPD_binY_contX"

#' Aggregate level patient data for binary outcome, continuous covariates
#'
#' This data set contains summaries of simulated patient covariate and outcome values.
#'
#' @name BC_ALD_binY_contX
#' @docType data
#' @format ## `y ~ PF_cont_1 + PF_cont_2 + trt + trt:(EM_cont_1 + EM_cont_2)`
#' \describe{
#' \item{variable}{String covariate or outcome name. From EM_cont_1, EM_cont_2, PF_cont_1, PF_cont_2, y.}
#' \item{statistic}{String summary statistic name. From mean, sd, sum, N}
#' \item{value}{Numeric value}
#' \item{trt}{Treatment (arm) name. From B, C}
#' }
#' @references Remiro‐Azocar A, Heath A, Baio G (2022)
#' @source Simulated data
#' @keywords datasets
#' @usage data(BC_ALD_binY_contX)
"BC_ALD_binY_contX"

#' Individual-level patient data for continuous outcome, mixed covariates
#'
#' This data set contains simulated patient covariate and outcome values.
#' Corresponds to ALD data set.
#'
#' @name AC_IPD_contY_mixedX
#' @docType data
#' @format ## `y ~ X1 + X3 + X4 + trt + trt:(X2 + X3 + X4)`
#' \describe{
#'    \item{id}{Numeric unique identifier}
#'    \item{X1}{Numeric prognostic factor continuous covariate}
#'    \item{X2}{Numeric prognostic factor and effect modifier binary covariate}
#'    \item{X3}{Numeric prognostic factor and effect modifier continuous covariate}
#'    \item{X4}{Numeric effect modifier binary covariate}
#'    \item{trt}{Factor treatment identifier. Levels A, C}
#'    \item{y}{Integer binary outcome}
#'    \item{true_eta}{Numeric linear predictor}
#' }
#' @references Remiro‐Azocar A, Heath A, Baio G (2022)
#' @source Simulated data
#' @keywords datasets
#' @usage data(AC_IPD_contY_mixedX)
"AC_IPD_contY_mixedX"

#' Aggregate level patient data for continuous outcome, mixed covariates
#'
#' This data set contains summaries of simulated patient covariate and outcome values.
#' Corresponds to IPD data set.
#'
#' @name BC_ALD_contY_mixedX
#' @docType data
#' @format ## `y ~ X1 + X3 + X4 + trt + trt:(X2 + X3 + X4)`
#' \describe{
#' \item{variable}{String covariate or outcome name. From X1, X2, X3, X4, y.}
#' \item{statistic}{String summary statistic name. From mean, sd, prob, sum, N}
#' \item{value}{Numeric value}
#' \item{trt}{Treatment (arm) name. From B, C}
#' }
#' @references Remiro‐Azocar A, Heath A, Baio G (2022)
#' @source Simulated data
#' @keywords datasets
#' @usage data(BC_ALD_contY_mixedX)
"BC_ALD_contY_mixedX"

#' Individual-level patient data for count outcome, continuous covariates
#'
#' This data set contains simulated patient covariate and outcome values.
#' Corresponds to ALD data set.
#'
#' @name AC_IPD_countY_contX
#' @docType data
#' @format ## `y ~ PF_cont_1 + PF_cont_2 + trt + trt:(EM_cont_1 + EM_cont_2)`
#' \describe{
#'    \item{id}{Numeric unique identifier}
#'    \item{PF_cont_1}{Numeric prognostic factor continuous covariate}
#'    \item{PF_cont_2}{Numeric prognostic factor continuous covariate}
#'    \item{EM_cont_1}{Numeric effect modifier continuous covariate}
#'    \item{EM_cont_2}{Numeric effect modifier continuous covariate}
#'    \item{trt}{Factor treatment identifier. Levels A, C}
#'    \item{y}{Integer non-negative count outcome}
#'    \item{true_eta}{Numeric linear predictor}
#' }
#' @references Remiro‐Azocar A, Heath A, Baio G (2022)
#' @source Simulated data
#' @keywords datasets
#' @usage data(AC_IPD_countY_contX)
"AC_IPD_countY_contX"

#' Aggregate level patient data for count outcome, continuous covariates
#'
#' This data set contains summaries of simulated patient covariate and outcome values.
#' Corresponds to IPD data set.
#'
#' @name BC_ALD_countY_contX
#' @docType data
#' @format ## `y ~ PF_cont_1 + PF_cont_2 + trt + trt:(EM_cont_1 + EM_cont_2)`
#' \describe{
#' \item{variable}{String covariate or outcome name. From EM_cont_1, EM_cont_2, PF_cont_1, PF_cont_2, y.}
#' \item{statistic}{String summary statistic name. From mean, sd, sum, N}
#' \item{value}{Numeric value}
#' \item{trt}{Treatment (arm) name. From B, C}
#' }
#' @references Remiro‐Azocar A, Heath A, Baio G (2022)
#' @source Simulated data
#' @keywords datasets
#' @usage data(BC_ALD_countY_contX)
"BC_ALD_countY_contX"
