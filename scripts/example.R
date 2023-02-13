# Example R code implementing for
# MAIC, STC, G-computation methods
# on a simulated example


library(boot)      # non-parametric bootstrap in MAIC and ML G-computation
library(copula)    # simulating BC covariates from Gaussian copula
library(rstanarm)  # fit outcome regression, draw outcomes in Bayesian G-computation

set.seed(555)

AC.IPD <- read.csv(here::here("data", "AC_IPD.csv"))  # AC patient-level data
BC.ALD <- read.csv(here::here("data", "BC_ALD.csv"))  # BC aggregate-level data

# marginal effect variance using the delta method
marginal_variance <- function(x)
  1/x$y.C.sum + 1/(x$N.C - x$y.C.sum) + 1/x$y.B.sum + 1/(x$N.B - x$y.B.sum)

# B vs C marginal treatment effect from reported event counts
marginal_treatment_effect <- function(x)
  log(x$y.B.sum*(x$N.C - x$y.C.sum)/(x$y.C.sum*(x$N.B - x$y.B.sum)))

### MAIC ###

# non-parametric bootstrap with 1000 resamples
boot.object <- boot::boot(data=AC.IPD, statistic=maic.boot, R=1000)

# bootstrap mean of marginal A vs C treatment effect estimate
hat.Delta.AC <- mean(boot.object$t)

# bootstrap variance of A vs C treatment effect estimate
hat.var.Delta.AC <- var(boot.object$t)

hat.Delta.BC <- marginal_treatment_effect(BC.ALD)

hat.var.Delta.BC <- marginal_variance(BC.ALD)

hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC
hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC

# construct Wald-type normal distribution-based confidence interval
ci.Delta.AB <- hat.Delta.AB + qnorm(c(0.025, 0.975))*sqrt(hat.var.Delta.AB)

### STC (conventional outcome regression) ###

# fit regression model of outcome on treatment and covariates
# IPD effect modifiers centred at the mean BC values
# purely prognostic variables are included but not centred
outcome.model <-
  glm(y ~ X3 + X4 + trt*I(X1-BC.ALD$mean.X1) + trt*I(X2-BC.ALD$mean.X2),
      data = AC.IPD,
      family = binomial)

# fitted treatment coefficient is relative A vs C conditional effect
hat.Delta.AC <- coef(outcome.model)["trt"]
# estimated variance for A vs C from model fit
hat.var.Delta.AC <- vcov(outcome.model)["trt", "trt"]
# B vs C marginal treatment effect estimated from reported event counts
hat.Delta.BC <- marginal_treatment_effect(BC.ALD)

hat.var.Delta.BC <- marginal_variance(BC.ALD)

hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC
hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC

# construct Wald-type normal distribution-based confidence interval
ci.Delta.AB <- hat.Delta.AB + qnorm(c(0.025, 0.975))*sqrt(hat.var.Delta.AB)

### Parametric G-computation with maximum-likelihood estimation ###

# matrix of pairwise correlations between IPD covariates
rho <- cor(AC.IPD[,c("X1","X2","X3","X4")])
# covariate simulation for BC trial using copula package
cop <-
  normalCopula(param = c(rho[1,2],rho[1,3],rho[1,4],rho[2,3],rho[2,4],rho[3,4]),
               dim = 4,
               dispstr = "un") # AC IPD pairwise correlations

# sample covariates from approximate joint distribution using copula
mvd <- mvdc(copula=cop, margins=c("norm", "norm", # Gaussian marginals
                                  "norm", "norm"),
            # BC covariate means and standard deviations
            paramMargins = list(list(mean=BC.ALD$mean.X1, sd=BC.ALD$sd.X1),
                                list(mean=BC.ALD$mean.X2, sd=BC.ALD$sd.X2),
                                list(mean=BC.ALD$mean.X3, sd=BC.ALD$sd.X3),
                                list(mean=BC.ALD$mean.X4, sd=BC.ALD$sd.X4)))
# simulated BC pseudo-population of size 1000
x_star <- as.data.frame(rMvdc(1000, mvd))
colnames(x_star) <- c("X1", "X2", "X3", "X4")

# non-parametric bootstrap with 1000 resamples
boot.object <- boot::boot(data=AC.IPD, statistic=gcomp_ml, R=1000)

# bootstrap mean of marginal A vs C treatment effect estimate
hat.Delta.AC <- mean(boot.object$t)
# bootstrap variance of A vs C treatment effect estimate
hat.var.Delta.AC <- var(boot.object$t)

hat.Delta.BC <- marginal_treatment_effect(BC.ALD)

hat.var.Delta.BC <- marginal_variance(BC.ALD)

# marginal treatment effect for A vs B
hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC
# variance for A vs B
hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC

# construct Wald-type normal distribution-based confidence interval
ci.Delta.AB <- hat.Delta.AB + qnorm(c(0.025, 0.975))*sqrt(hat.var.Delta.AB)

### Bayesian G-computation with MCMC ###

# matrix of pairwise correlations between IPD covariates
rho <- cor(AC.IPD[,c("X1","X2","X3","X4")])
#  covariate simulation for BC trial using copula package
cop <-
  normalCopula(param = c(rho[1,2],rho[1,3],rho[1,4],rho[2,3],rho[2,4],rho[3,4]),
               dim = 4,
               dispstr = "un") # AC IPD pairwise correlations

# sample covariates from approximate joint distribution using copula
mvd <- mvdc(copula=cop, margins=c("norm", "norm", # Gaussian marginals
                                  "norm", "norm"),
            # BC covariate means and standard deviations
            paramMargins = list(list(mean=BC.ALD$mean.X1, sd=BC.ALD$sd.X1),
                                list(mean=BC.ALD$mean.X2, sd=BC.ALD$sd.X2),
                                list(mean=BC.ALD$mean.X3, sd=BC.ALD$sd.X3),
                                list(mean=BC.ALD$mean.X4, sd=BC.ALD$sd.X4)))
# simulated BC pseudo-population of size 1000
x_star <- as.data.frame(rMvdc(1000, mvd))
colnames(x_star) <- c("X1", "X2", "X3", "X4")

# outcome logistic regression fitted to IPD using MCMC (Stan)
outcome.model <- stan_glm(y ~ X3+X4+trt*X1+trt*X2,
                          data = AC.IPD,
                          family = binomial,
                          algorithm = "sampling",
                          iter = 4000, warmup = 2000, chains = 2)
# counterfactual datasets
data.trtA <- data.trtC <- x_star

# intervene on treatment while keeping set covariates fixed
data.trtA$trt <- 1 # everyone receives treatment A
data.trtC$trt <- 0 # all observations receive C

# draw binary responses from posterior predictive distribution
y.star.A <- posterior_predict(outcome.model, newdata=data.trtA)
y.star.C <- posterior_predict(outcome.model, newdata=data.trtC)

# compute marginal log-odds ratio for A vs C for each MCMC sample
# by transforming from probability to linear predictor scale
hat.delta.AC <- qlogis(rowMeans(y.star.A)) - qlogis(rowMeans(y.star.C))
hat.Delta.AC <- mean(hat.delta.AC)     # average over samples
hat.var.Delta.AC <- var(hat.delta.AC)  # sample variance

# B vs C from reported aggregate event counts in contingency table
hat.Delta.BC <-
  with(BC.ALD, log(y.B.sum*(N.C - y.C.sum)/(y.C.sum*(N.B - y.B.sum))))

hat.var.Delta.BC <- marginal_variance(BC.ALD)

# marginal treatment effect for A vs B
hat.Delta.AB <- hat.Delta.AC - hat.Delta.BC
# A vs B variance
hat.var.Delta.AB <- hat.var.Delta.AC + hat.var.Delta.BC

# construct Wald-type normal distribution-based confidence interval
ci.Delta.AB <- hat.Delta.AB + qnorm(c(0.025, 0.975))*sqrt(hat.var.Delta.AB)
