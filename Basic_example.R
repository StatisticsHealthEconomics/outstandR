## ----include = FALSE------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)


## ----setup, warning=FALSE, message=FALSE----------------------------------
library(boot)      # non-parametric bootstrap in MAIC and ML G-computation
library(copula)    # simulating BC covariates from Gaussian copula
library(rstanarm)  # fit outcome regression, draw outcomes in Bayesian G-computation
library(mimR)


## ----load-data------------------------------------------------------------
set.seed(555)

AC.IPD <- read.csv(here::here("data", "AC_IPD.csv"))  # AC patient-level data
BC.ALD <- read.csv(here::here("data", "BC_ALD.csv"))  # BC aggregate-level data


## -------------------------------------------------------------------------
head(AC.IPD)


## -------------------------------------------------------------------------
BC.ALD


## ----mimR_maic-------------------------------------------------
mimR_maic <- mimR(AC.IPD, BC.ALD, strategy = strategy_maic())


## ----mimR_maic-print-------------------------------------------
mimR_maic


## ----mimR_stc--------------------------------------------------
mimR_stc <- mimR(AC.IPD, BC.ALD, strategy = strategy_stc())
mimR_stc


## ----mimR_gcomp_ml---------------------------------------------
# mimR_gcomp_ml <- mimR(AC.IPD, BC.ALD, strategy = strategy_gcomp_ml())


## ----mimR_gcomp_stan-------------------------------------------
mimR_gcomp_stan <- mimR(AC.IPD, BC.ALD, strategy = strategy_gcomp_stan())
mimR_gcomp_stan

