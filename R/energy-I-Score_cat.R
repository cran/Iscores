
#' Internal function for changing factors to numerical
#'
#' A supplementary function for data management
#'
#' @param factor_col a factor column
#'
#' @details This function converts factor variables to numeric variables.
#'
#' @keywords internal
#'
#'

factor_to_numeric <- function(factor_col) {
  as.numeric(levels(factor_col))[factor_col]
}

#' One hot encoding
#'
#' A supplementary function for one-hot encoding
#'
#' @param dat a data containing some factor but numeric columns.
#'
#' @details This function converts factor variables into one-hot encoding
#'
#' @keywords internal

factor_to_onehot <- function(dat) {

  dat <- data.frame(dat)
  factor_columns <- which(sapply(as.data.frame(dat), is.factor))

  if(length(factor_columns) == 0) return(dat)

  n_levels <- c()

  for(ith_factor in factor_columns){
    n_levels <- c(n_levels, length(levels(dat[, ith_factor])))
    dat <- cbind(dat, do_one_hot(dat[, ith_factor]))
  }

  dat <- dat[, -factor_columns]
  attr(dat, "mask") <- rep(c(0, factor_columns),
                           times = c(ncol(dat) - sum(n_levels),  n_levels))
  attr(dat, "column_names") <- names(factor_columns)

  dat
}



#' Convert a factor vector to one-hot encoding
#'
#' @description Converts a factor vector into a one-hot encoded matrix with one
#' column per factor level.
#'
#' @param vec A factor vector to be encoded.
#'
#' @return A numeric matrix with one row per element of `vec` and one column per
#' factor level. Column names are prefixed with `"level_"`.
#'
#' @details
#' Missing values in `vec` are preserved as rows containing `NA` values.
#'
#' @keywords internal
#'

do_one_hot <- function(vec) {

  NA_mat <- matrix(NA, nrow = length(vec), ncol = length(levels(vec)))

  if(ncol(NA_mat) == 1) {
    NA_mat[, 1] <- vec
  } else {
    mm <- cbind(stats::model.matrix(~vec), 0)
    mm[, 1] <-  mm[, 1] - rowSums(data.frame(mm[, -1]))
    NA_mat[as.numeric(rownames(mm)), ] <- mm[, - ncol(mm)]
  }
  colnames(NA_mat) <- paste0("level_", sort(levels(vec)))

  NA_mat
}




#' energy-I-Score for imputation of mixed data (categorical and numerical)
#'
#' @importFrom scoringRules crps_sample
#' @importFrom pbapply pblapply
#' @importFrom stats model.matrix
#' @inheritParams energy_Iscore_num
#'
#' @return a numerical value denoting weighted Imputation Score obtained for
#' provided imputation function and a table with scores and weights calculated
#' for particular columns.
#'
#' @details
#' The categorical variables should be stored as factors. If you need additional
#' conversion of the data (for example one-hot encoding) for imputation, please,
#' implement everything within \code{imputation_func} parameter. You can use
#' \code{miceDRF:::onehot_to_factor} and \code{miceDRF:::factor_to_onehot}
#' functions.
#'
#' @references
#' This method is described in detail in:
#'
#' Näf, J., Grzesiak, K., and Scornet, E. (2025).
#' How to rank imputation methods?
#' arXiv preprint.
#' \doi{10.48550/arXiv.2507.11297}.
#'
#' @keywords internal
#'

energy_Iscore_cat <- function(X,
                              imputation_func,
                              X_imp = imputation_func(X),
                              multiple = TRUE,
                              N = 50,
                              max_length = NULL,
                              skip_if_needed = TRUE,
                              scale = FALSE,
                              n_cores = 1,
                              silent = TRUE){

  warnings_vec <- c()

  N <- ifelse(multiple, N, 1)

  n <- nrow(X)
  missings_per_col <- colSums(is.na(X))

  ## Missings pattern
  M <- is.na(X)

  dim_with_NA <- missings_per_col > 0

  if (is.null(max_length)) max_length <- sum(dim_with_NA)

  if (sum(dim_with_NA) < max_length){
    if(!silent) {
      warning("max_length is larger than the total number of columns with missing values!")
    }
    max_length <- sum(dim_with_NA)
  }

  factor_columns <- which(sapply(X, is.factor))

  cols_to_iterate <- intersect(order(missings_per_col, decreasing = TRUE),
                               which(dim_with_NA))[1:max_length]

  scores_dat <- pbmcapply::pbmclapply(cols_to_iterate, function(j) {

    weight <- (missings_per_col[j] / n) * ((n - missings_per_col[j]) / n)

    if(missings_per_col[j] < 10) {
      if(!silent) {
        warning('Sample size of missing and nonmissing too small for nonparametric distributional regression, setting to NA')
      }
      return(data.frame(column_id = j,
                        weight = weight,
                        score = NA,
                        n_columns_used = NA)) # return score = NA
    }

    observed_j_for_train <- !M[, j]

    # Fully observed columns except j
    Oj <- colSums(is.na(X[observed_j_for_train, ][, -j])) == 0

    if(!any(Oj)) {

      if(skip_if_needed) {

        Oj_candidates <- M[, -j]
        max_obs_Ojs <- colSums(!Oj_candidates[observed_j_for_train, ])
        observed_j_for_train <- !Oj_candidates[, which.max(max_obs_Ojs)] & !M[, j]
        if(!silent) {
          message(paste0("No complete variables for training column ", j,
                         ". Skipping some observations."))
        }

        Oj <- colSums(is.na(X[observed_j_for_train, ][, -j])) == 0

      } else {
        if(!silent) {
          warning("Oj was empty. There was no complete column for training.")
        }
        return(data.frame(column_id = j,
                          weight = weight,
                          score = NA,
                          n_columns_used = sum(Oj))) # return score = NA
      }
    }

    # Only take those that are fully observed H for all observed values of X_j
    X_imp_0 <- X_imp[observed_j_for_train, ]
    X_test <- X_imp_0[, -j][, Oj]
    Y_test <- X_imp_0[, j]

    # Only take those that are fully observed H for all missing values of X_j
    X_imp_1 <- X_imp[!observed_j_for_train, ]
    X_train <- X_imp_1[, -j][, Oj]
    Y_train <- X_imp_1[, j]

    if(sum(Oj) > 2) {
      names(X_test) <- paste0("1234", 1:ncol(X_test))
      names(X_train) <- names(X_test)
    }

    if(j %in% factor_columns) {

      X_artificial <- rbind(data.frame(y = as.factor(NA), X = X_test),
                            data.frame(y = Y_train, X = X_train))

      Y_test <- factor_to_onehot(Y_test)
    } else {

      X_artificial <- rbind(data.frame(y = NA_real_, X = X_test),
                            data.frame(y = as.numeric(Y_train), X = X_train))
    }

    imputation_list <- lapply(1:N, function(ith_imputation) {

      imputed <- try({imputation_func(X_artificial)})

      if(inherits(imputed, "try-error") | any(is.na(imputed)))
        return(NA)

      if(j %in% factor_columns) {
        res <- imputed[1:nrow(Y_test), 1]

        res <- factor_to_onehot(res)

      } else {
        res <- imputed[1:length(Y_test), 1]
      }

      res

    })

    if(length(imputation_list[!sapply(imputation_list,
                                      function(x) all(is.na(x)))]) < N) {
      if(!silent) {
        warning(sprintf("Unsuccessful imputation! Imputation function is unstable!
              Returning NA for column %i.", j))
      }
      return(data.frame(column_id = j, weight = weight, score = NA,
                        n_columns_used = sum(Oj))) # return score = NA
    }

    if(j %in% factor_columns) {
      Y_test <- factor_to_onehot(Y_test)

      if(!all(colnames(Y_test) %in% colnames(imputation_list[[1]]))) {

        imputation_list <- lapply(imputation_list, function(ith) {
          missing_cols <- setdiff(colnames(Y_test), colnames(ith))

          zeroes <- matrix(0, nrow = nrow(ith), ncol = length(missing_cols))
          colnames(zeroes) <- missing_cols

          tmp <- cbind(ith, zeroes)

          tmp[, colnames(Y_test)]
        })
      }

      Y_matrix <- do.call(cbind, imputation_list)

      score_j <- mean(sapply(1:nrow(Y_test), function(ith_obs) {
        scoringRules::es_sample(as.numeric(unlist(Y_test[ith_obs, ])),
                                t(matrix(as.numeric(Y_matrix[ith_obs, ]),
                                         ncol = ncol(Y_test), byrow = TRUE)))
      }))

    } else {
      Y_matrix <- do.call(cbind, imputation_list)

      if(scale) {
        Y_test <- (Y_test - mean(Y_test)) / sd(Y_test)
        Y_matrix <- (Y_matrix - mean(Y_test)) / sd(Y_test)
      }

      score_j <- mean(scoringRules::crps_sample(y = Y_test, dat = Y_matrix))
    }

    data.frame(column_id = j,
               weight = weight,
               score = score_j,
               n_columns_used = sum(Oj))

  }, mc.cores = n_cores)


  scores_dat <- do.call(rbind, scores_dat)

  weighted_score <- sum(scores_dat[["score"]] * scores_dat[["weight"]] /
                          (sum(scores_dat[["weight"]], na.rm = TRUE)),
                        na.rm = TRUE)

  weighted_score <- ifelse(all(is.na(scores_dat[["score"]])), NA, weighted_score)

  attr(weighted_score, "dat") <- scores_dat
  weighted_score
}





