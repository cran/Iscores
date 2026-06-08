
#' Energy distance
#'
#' Calculating energy distance/statistic.
#'
#' @importFrom energy eqdist.e
#'
#' @param X a complete original dataset (without missing values).
#' @param X_imp an imputed dataset
#' @param rescale a logical, indicating whether the returned value should be
#' rescaled. Default to \code{FALSE}. See "details" section for more information.
#'
#' @details This function uses the \link[energy]{eqdist.e} function. According
#' to this implementation, by default, the function returns the energy statistic
#' which is given by
#' \deqn{E(X, Y) = \frac{nm}{n + m} \hat{\varepsilon}{(X, Y)},}
#' where \eqn{\hat{\varepsilon}{(X, Y)}} is the raw energy distance. To
#' obtain raw energy distance use \code{rescale = TRUE}.
#'
#' @return A numeric value giving the energy distance between the original
#' dataset and the imputed dataset.
#'
#' @examples
#' X <- matrix(rnorm(100), nrow = 25)
#' X_imp <- matrix(rnorm(100), nrow = 25)
#' edistance(X, X_imp)
#'
#' @export

edistance <- function(X, X_imp, rescale = FALSE){

  rescale <- ((nrow(X) + nrow(X_imp)) / (nrow(X) * nrow(X_imp))) ^ rescale

  energy::eqdist.e(rbind(X, X_imp), c(nrow(X), nrow(X_imp))) * rescale
}



