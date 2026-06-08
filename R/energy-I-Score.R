
#' @title Calculates score for one imputation function
#'
#' @importFrom scoringRules crps_sample
#' @importFrom pbapply pblapply
#' @importFrom stats sd
#'
#' @param X data containing missing values denoted with NA's.
#' @param X_imp imputed dataset of the same size as \code{X}. It's \code{NULL}
#' by default meaning that it will be obtained by imputation of  \code{X} using
#' the \code{imputation_func}.
#' @param imputation_func a function that imputes data.
#' @param N a numeric value. Number of samples from imputation distribution H.
#' Default to 50.
#' @param max_length Maximum number of variables \eqn{X_j} to consider, can
#' speed up the code. Default to \code{NULL} meaning that all the columns will
#' be taken under consideration.
#' @param multiple a logical indicating whether provided imputation method is a
#' multiple imputation approach (i.e. it generates different values to impute
#' for each call). Default to TRUE. Note that if multiple equals to FALSE, N is
#' automatically set to 1.
#' @param skip_if_needed logical, indicating whether some observations should be
#' skipped to obtain complete columns for scoring. If FALSE, NA will be returned
#' for column with no observed variable for training.
#' @param scale a logical value. If TRUE, each variable is scaled in the score.
#' @param n_cores a number of cores for parallelization.
#' @param silent logical indicating whether warnings and messages should be
#' printed.
#'
#' @return a numerical value denoting weighted Imputation Score obtained for
#' provided imputation function and a table with scores and weights calculated
#' for particular columns.
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

energy_Iscore_num <- function(X,
                              imputation_func,
                              X_imp = imputation_func(X),
                              multiple = TRUE,
                              N = 50,
                              max_length = NULL,
                              skip_if_needed = TRUE,
                              scale = FALSE,
                              n_cores = 1,
                              silent = TRUE){

  N <- ifelse(multiple, N, 1)

  n <- nrow(X)

  missings_per_col <- colSums(is.na(X))

  M <- is.na(X)   # Missings pattern
  dim_with_NA <- missings_per_col > 0

  if (is.null(max_length))
    max_length <- sum(dim_with_NA)

  if (sum(dim_with_NA) < max_length){

    if(!silent) {
      warning(sprintf("max_length is larger than the total number of columns with missing values (%i)!
                    Setting max_length to %i", sum(dim_with_NA), sum(dim_with_NA)))
    }

    max_length <- sum(dim_with_NA)
  }

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
                        n_columns_used = NA))
    }

    observed_j_for_train <- !M[, j]

    # Fully observed columns except j
    Oj <- colSums(is.na(X[observed_j_for_train, ][, -j, drop = FALSE])) == 0

    if(!any(Oj)) {

      if(skip_if_needed) {
        if(!silent) {
          message(paste0("No complete variables for training column ", j,
                         ". Skipping some observations."))
        }

        Oj_candidates <- M[, -j]
        max_obs_Ojs <- colSums(!Oj_candidates[observed_j_for_train, ])
        observed_j_for_train <- !Oj_candidates[, which.max(max_obs_Ojs)] & !M[, j]

        Oj <- colSums(is.na(X[observed_j_for_train, ][, -j, drop = FALSE])) == 0

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
    X_test <- X_imp_0[, -j, drop = FALSE][, Oj]
    Y_test <- X_imp_0[, j]

    # Only take those that are fully observed H for all missing values of X_j
    X_imp_1 <- X_imp[!observed_j_for_train, ]
    X_train <- X_imp_1[, -j, drop = FALSE][, Oj]
    Y_train <- X_imp_1[, j]

    # Train DRF on imputed data
    X_artificial <- as.data.frame(rbind(cbind(y = NA, X_test),
                                        cbind(y = Y_train, X_train)))

    imputation_list <- lapply(1:N, function(ith_imputation) {

      imputed <- try({imputation_func(X_artificial)})

      if(inherits(imputed, "try-error") | any(is.na(imputed)))
        return(NULL)

      imputed[1:length(Y_test), 1]
    })

    if(sum(!sapply(imputation_list, is.null)) < N) {
      if(!silent) {
        warning(sprintf("Unsuccessful imputation! Imputation function is unstable!
              Returning NA for column %i.", j))
      }

      return(data.frame(column_id = j,
                        weight = weight,
                        score = NA,
                        n_columns_used = sum(Oj)))
    }

    Y_matrix <- do.call(cbind, imputation_list)

    if(scale) {
      Y_test <- (Y_test - mean(Y_test)) / sd(Y_test)
      Y_matrix <- (Y_matrix - mean(Y_test)) / sd(Y_test)
    }

    score_j <- mean(scoringRules::crps_sample(y = Y_test, dat = Y_matrix))

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





