#' Generate random data with MCAR missing values
#'
#' Generates a numerical dataset consisting of independent standard normal
#' variables and introduces missing values according to a Missing Completely
#' at Random (MCAR) mechanism.
#'
#' @param n Number of observations.
#' @param p Number of numerical variables.
#' @param ratio Proportion of entries to replace with missing values.
#'
#' @return A data frame with \code{n} rows and \code{p} numerical variables
#' containing missing values.
#'
#' @examples
#' X <- random_mcar_data(100, 3, ratio = 0.2)
#' head(X)
#'
#' @export

random_mcar_data <- function(n, p, ratio = 0.2) {
  X <- matrix(stats::rnorm(n * p), nrow = n)
  X[stats::runif(n * p) <= ratio] <- NA
  data.frame(X)
}



#' Generate random mixed data with MCAR missing values
#'
#' Generates a mixed dataset containing independent standard normal variables
#' and categorical variables, then introduces missing values according to a
#' Missing Completely at Random (MCAR) mechanism.
#'
#' @param n Number of observations.
#' @param p Number of numerical variables.
#' @param n_fac Number of categorical variables.
#' @param ratio Proportion of entries to replace with missing values.
#'
#' @return A data frame containing \code{p} numerical variables and
#' \code{n_fac} factor variables with missing values.
#'
#' @examples
#' X <- random_mcar_mixed_data(100, 3, n_fac = 2, ratio = 0.2)
#' str(X)
#'
#' @export

random_mcar_mixed_data <- function(n, p, n_fac = 1, ratio = 0.2) {
  X <- matrix(stats::rnorm(n * p), nrow = n)

  factors <- apply(
    matrix(sample(1:4, n * n_fac, replace = TRUE), nrow = n),
    2,
    function(x) factor(x)
  )

  factors[stats::runif(n * n_fac) <= ratio] <- NA
  X[stats::runif(n * p) <= ratio] <- NA

  X <- data.frame(X)
  factors <- data.frame(factors)

  X <- cbind(X, factors)
  colnames(X) <- paste0("col", 1:ncol(X))

  for (i in (ncol(X) - n_fac + 1):ncol(X)) {
    X[[i]] <- factor(X[[i]])
  }

  X
}


#' Standard exponential imputation
#'
#' Imputes all missing values by independent draws from an exponential
#' distribution with rate 1.
#'
#' @importFrom stats rexp
#'
#' @param X_miss A data set containing missing values.
#'
#' @return A completed data set with all missing values replaced by draws
#' from an \code{Exp(1)} distribution.
#'
#' @examples
#' X <- random_mcar_data(100, 3)
#' X_imp <- exp_imputation(X)
#'
#' @export

exp_imputation <- function(X_miss) {
  X_miss[is.na(X_miss)] <- stats::rexp(sum(is.na(X_miss)))

  X_miss
}


#' Standard normal imputation
#'
#' Imputes all missing values by independent draws from a standard normal
#' distribution.
#'
#' @importFrom stats rnorm
#'
#' @param X_miss A data set containing missing values.
#'
#' @return A completed data set with all missing values replaced by draws
#' from a \eqn{N(0,1)} distribution.
#'
#' @examples
#' X <- random_mcar_data(100, 3)
#' X_imp <- norm_imputation(X)
#'
#' @export

norm_imputation <- function(X_miss) {
  X_miss[is.na(X_miss)] <- stats::rnorm(sum(is.na(X_miss)))

  X_miss
}

#' Median/mode imputation
#'
#' Imputes numerical variables using their median and categorical variables
#' using their most frequent observed category.
#'
#' @param X_miss A data set containing missing values.
#'
#' @return A completed data set.
#'
#' @examples
#' X <- random_mcar_mixed_data(100, 3, n_fac = 1)
#' X_imp <- median_mode_imputation(X)
#'
#' @export

median_mode_imputation <- function(X_miss) {
  for (col in names(X_miss)) {
    if (is.numeric(X_miss[[col]])) {
      med <- stats::median(X_miss[[col]], na.rm = TRUE)
      X_miss[[col]][is.na(X_miss[[col]])] <- med
    }
    else if (is.factor(X_miss[[col]])) {
      mode_val <- names(sort(table(X_miss[[col]]), decreasing = TRUE))[1]
      X_miss[[col]][is.na(X_miss[[col]])] <- mode_val
      X_miss[[col]] <- factor(X_miss[[col]])
    }
  }

  X_miss
}

