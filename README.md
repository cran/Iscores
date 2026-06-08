# Iscores

<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/Iscores)](https://CRAN.R-project.org/package=Iscores)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![R-CMD-check](https://github.com/KrystynaGrzesiak/Iscores/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/KrystynaGrzesiak/Iscores/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->


`Iscores` provides scoring rules for evaluating and comparing imputation methods.

The package implements the methodology introduced in Näf et al. (2022) and Näf, Grzesiak, and Scornet (2025).
The package supports:

- numerical datasets,
- mixed numerical and categorical data,
- deterministic imputations,
- stochastic and multiple imputations,
- comparison of several imputation methods.

For more details about the energy-I-Score check our vignettes:

- [Energy-I-Score: Implementation Details](https://krystynagrzesiak.github.io/Iscores/articles/About_IScore.html)
- [Energy-I-Score: First Steps](https://krystynagrzesiak.github.io/Iscores/articles/Example_IScore.html)


---

## Installation

### CRAN

```r
install.packages("Iscores")
```

### Development version

```r
install.packages("devtools")

devtools::install_github("missValTeam/Iscores")
```

---

## Basic workflow

The package evaluates user-defined imputation methods.

An imputation function must:

- accept a dataset with missing values,
- return a completed dataset with the same dimensions.

Below we define a simple zero-imputation method.

```r
library(Iscores)

impute_zero <- function(X) {

  X[is.na(X)] <- 0

  X
}
```

We now generate example data with missing values.

```r
set.seed(10)

X <- Iscores:::random_mcar_data(100, 4)

head(X)
```

---

## Energy-I-Score

The `energy_IScore()` function evaluates the quality of an imputation method.

```r
sc <- energy_IScore(
  X = X,
  imputation_func = impute_zero,
  N = 10,
  silent = TRUE
)

sc
```

Detailed variable-level results are stored as an attribute:

```r
attr(sc, "dat")
```

---

## DR-I-Score

The package also provides the density-ratio based DR-I-Score.

```r
sc_dr <- DR_IScore(
  X = X,
  imputation_func = impute_zero,
  m = 3,
  n_proj = 10,
  n_trees_per_proj = 2,
  n_cores = 1
)

sc_dr
```

---

## Comparing imputation methods

Several methods can be compared simultaneously using `compare_Iscores()`.

```r
library(mice)

impute_mice_norm <- function(X) {

  imp <- mice(
    X,
    m = 1,
    method = "norm",
    maxit = 5,
    printFlag = FALSE
  )

  complete(imp)
}

impute_mice_rf <- function(X) {

  imp <- mice(
    X,
    m = 1,
    method = "rf",
    maxit = 5,
    printFlag = FALSE
  )

  complete(imp)
}

methods_list <- list(
  zero = impute_zero,
  norm = impute_mice_norm,
  rf = impute_mice_rf
)

compare_Iscores(
  X = X,
  methods_list = methods_list,
  score = c("energy_IScore", "DR_IScore"),
  N = 10,
  m = 3,
  silent = TRUE
)
```

---

## Documentation

See the vignette for a complete introduction:

```r
vignette("Example_IScore")
```

---

## References

Näf, Jeffrey, Krystyna Grzesiak, and Erwan Scornet. 2025. “How to Rank Imputation Methods?” 
https://arxiv.org/abs/2507.11297.

Näf, Jeffrey, Meta-Lina Spohn, Loris Michel, and Nicolai Meinshausen. 2022. “Imputation Scores.” 
https://arxiv.org/abs/2106.03742.
