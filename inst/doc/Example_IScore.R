## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----eval = FALSE-------------------------------------------------------------
# install.packages("Iscores")

## ----eval = FALSE-------------------------------------------------------------
# if (!requireNamespace("devtools", quietly = TRUE)) {
#   install.packages("devtools")
# }
# devtools::install_github("missValTeam/Iscores")

## ----message=FALSE, warning=FALSE---------------------------------------------
library(Iscores)

## -----------------------------------------------------------------------------
set.seed(10)

X <- random_mcar_data(100, 4)

head(X)

## -----------------------------------------------------------------------------
impute_zero <- function(X) { 
  
  X[is.na(X)] <- 0
  
  return(X) 
}

## -----------------------------------------------------------------------------
sc <- energy_IScore(X = X, imputation_func = impute_zero)

sc

## -----------------------------------------------------------------------------
attr(sc, "dat")

## -----------------------------------------------------------------------------
sum(attr(sc, "dat")[["score"]] * attr(sc, "dat")[["weight"]]) / sum(attr(sc, "dat")[["weight"]])

## ----warning=FALSE, message=FALSE---------------------------------------------
energy_IScore(X = X, imputation_func = impute_zero, N = 5)

## ----warning=FALSE, message=FALSE---------------------------------------------
energy_IScore(X = X, imputation_func = impute_zero, max_length = 2)

## -----------------------------------------------------------------------------
set.seed(10)

X_cat <- random_mcar_mixed_data(100, 4)

head(X_cat)

## -----------------------------------------------------------------------------
impute_mean_mode <- median_mode_imputation

## -----------------------------------------------------------------------------
energy_IScore(X = X_cat, imputation_func = impute_mean_mode)

## -----------------------------------------------------------------------------
sc_dr <- DR_IScore(X = X,
                   imputation_func = impute_zero,
                   m = 3,
                   n_proj = 10,
                   n_trees_per_proj = 2,
                   n_cores = 1)

sc_dr

## ----warning=FALSE, message=FALSE---------------------------------------------
library(mice)


impute_mice_norm <- function(X) {
  imp <- mice(X, m = 1, method = "norm", maxit = 5, printFlag = FALSE)
  
  complete(imp)
}

impute_mice_rf <- function(X) {
  imp <- mice(X, m = 1, method = "rf", maxit = 5, printFlag = FALSE)
  
  complete(imp)
}

## -----------------------------------------------------------------------------
methods_list <- list(zero = impute_zero,
                     mice_norm = impute_mice_norm,
                     mice_rf = impute_mice_rf)

## -----------------------------------------------------------------------------
sc_comparison <- compare_Iscores(X = X,
                                 methods_list = methods_list,
                                 score = "energy_IScore",
                                 N = 10,
                                 silent = TRUE)

sc_comparison

## -----------------------------------------------------------------------------
comparison_all <- compare_Iscores(X = X,
                                  methods_list = methods_list,
                                  score = c("energy_IScore", "DR_IScore"),
                                  N = 10,
                                  m = 3,
                                  n_proj = 10,
                                  n_trees_per_proj = 2,
                                  silent = TRUE)

comparison_all

## -----------------------------------------------------------------------------

X_observed <- matrix(rnorm(2000), ncol = 4)  

X_miss <- X_observed
X_miss[runif(nrow(X_miss) * ncol(X_miss)) < 0.2] <- NA

edistance(X_observed, impute_zero(X_miss))

edistance(X_observed, impute_mice_norm(X_miss))


