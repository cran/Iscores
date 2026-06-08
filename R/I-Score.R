
#' @title Calculates Imputation Score for imputation function
#'
#' @inheritParams energy_Iscore_num
#'
#' @details This function relies on functions \link[Iscores]{energy_Iscore_num}
#' and \link[Iscores]{energy_Iscore_cat}. Depending on the presence of factor-type
#' data, these functions compute a score either for purely numerical data or for
#' mixed data types.
#'
#' If you want to compute the score for numerical data, make sure that the
#' dataset does not contain any factor-type variables.
#'
#' If you want to compute the score for categorical data, ensure that all
#' categorical variables are preserved as factors.
#'
#' If your imputation method does not support categorical variables represented
#' as factors, implement a wrapper function that handles the appropriate data
#' type conversions before and after imputation.
#'
#' @return a numerical value denoting weighted Imputation Score obtained for
#' provided imputation function and a table with scores and weights calculated
#' for particular columns.
#'
#' @examples
#' set.seed(111)
#' X <- random_mcar_data(100, 4)
#' imputation_func <- exp_imputation
#' energy_IScore(X, imputation_func)
#'
#' X <-  random_mcar_mixed_data(100, 4, 2)
#' imputation_func <- median_mode_imputation
#' energy_IScore(X, imputation_func)
#'
#' @references
#'
#' Näf, J., Grzesiak, K., and Scornet, E. (2025).
#' How to rank imputation methods?
#' arXiv preprint.
#' \doi{10.48550/arXiv.2507.11297}.
#'
#' @export
#'


energy_IScore <- function(X,
                          imputation_func,
                          X_imp = NULL,
                          multiple = TRUE,
                          N = 50,
                          max_length = NULL,
                          skip_if_needed = TRUE,
                          scale = FALSE,
                          n_cores = 1,
                          silent = TRUE) {

  X <- as.data.frame(X, check.names = FALSE)

  categoricals <- sapply(X, function(i) is.factor(i))

  if(any(categoricals)) {
    message( "Factor variables detected.
             Calculating the energy-I-Score for mixed data.")
    mixed <- TRUE
  } else {
    mixed <- FALSE
  }

  if(!is.function(imputation_func))
    stop("Imputation_func must be a function!")

  if(is.null(X_imp)) {
    X_imp <- try({imputation_func(X)})

    if(inherits(X_imp, "try-error") | any(is.na(X_imp)))
      stop("Errored imputing X using provided imputation_func!")

    X_imp <- as.data.frame(X_imp, check.names = FALSE)
  }



  if(mixed) {
    score <- energy_Iscore_cat(X, imputation_func, X_imp, multiple, N,
                               max_length, skip_if_needed, scale, n_cores,
                               silent)
  } else {
    score <- energy_Iscore_num(X, imputation_func, X_imp, multiple, N,
                               max_length, skip_if_needed, scale, n_cores,
                               silent)
  }

  score
}

#' @title Calculates IScores for multiple imputation functions
#'
#' @inheritParams energy_Iscore_num
#'
#' @param score a vector of names of scores to calculate. It can be
#' \code{"energy_IScore"} and \code{"DR_IScore"}.
#' @param methods_list a named list of imputing functions.
#' @param ... other arguments to be passed to  \link[Iscores]{energy_IScore} or
#'  \link[Iscores]{DR_IScore}
#'
#'
#' @return a vector of IScores for provided methods
#'
#' @examples
#' set.seed(111)
#' X <- random_mcar_data(100, 3, 0.2)
#' methods_list <- list(exp = exp_imputation,
#'                        norm = norm_imputation)
#' compare_Iscores(X, methods_list = methods_list, m = 2,
#'                 n_proj = 10, n_trees_per_proj = 2 )
#'
#' @export
#'

compare_Iscores <- function(X,
                            methods_list,
                            score = c("energy_IScore", "DR_IScore"),
                            ...) {

  score <- match.arg(score, c("energy_IScore", "DR_IScore"), several.ok = TRUE)

  do.call(rbind,
          lapply(score, function(ith_score) {

            score_fun <- get(ith_score)
            methods <- names(methods_list)

            do.call(rbind, lapply(seq_along(methods_list), function(ith_method) {

              message(sprintf("Calculating the %s for method %s ...", ith_score,
                              names(methods_list)[ith_method]))

              imputation_func <- methods_list[[ith_method]]

              add_args <- list(...)

              args <- c(list(X = X, imputation_func = imputation_func),
                        add_args[names(add_args) %in% names(formals(score_fun))])

              score <- do.call(score_fun, args)

              data.frame(score = as.numeric(score),
                         score_name = ith_score,
                         method = methods[ith_method])
            }))

          }))


}

