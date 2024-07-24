## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup, warning=FALSE, message=FALSE--------------------------------------
library(boot)      # non-parametric bootstrap in MAIC and ML G-computation
library(copula)    # simulating BC covariates from Gaussian copula
library(rstanarm)  # fit outcome regression, draw outcomes in Bayesian G-computation
library(ModStanR)

## ----load-data----------------------------------------------------------------
set.seed(555)

AC.IPD <- read.csv(here::here("raw-data", "AC_IPD.csv"))  # AC patient-level data
BC.ALD <- read.csv(here::here("raw-data", "BC_ALD.csv"))  # BC aggregate-level data

## -----------------------------------------------------------------------------
head(AC.IPD)

## -----------------------------------------------------------------------------
BC.ALD

## -----------------------------------------------------------------------------
lin_form <- as.formula("y ~ X3 + X4 + trt*X1 + trt*X2")

## ----ModStanR_maic------------------------------------------------------------
ModStanR_maic <- ModStanR(AC.IPD, BC.ALD, strategy = strategy_maic(formula = lin_form))

## ----ModStanR_maic-print------------------------------------------------------
ModStanR_maic

## ----ModStanR_stc-------------------------------------------------------------
ModStanR_stc <- ModStanR(AC.IPD, BC.ALD, strategy = strategy_stc(formula = lin_form))
ModStanR_stc

## ----ModStanR_gcomp_ml--------------------------------------------------------
# ModStanR_gcomp_ml <- ModStanR(AC.IPD, BC.ALD, strategy = strategy_gcomp_ml(formula = lin_form))

## ----ModStanR_gcomp_stan------------------------------------------------------
ModStanR_gcomp_stan <- ModStanR(AC.IPD, BC.ALD, strategy = strategy_gcomp_stan(formula = lin_form))
ModStanR_gcomp_stan

