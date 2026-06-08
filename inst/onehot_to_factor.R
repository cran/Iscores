#' One hot encoding
#'
#' A supplementarty function for one-hot encoding
#'
#' @param onehot_dat a data coded with \code{factor_to_onehot} function.
#'
#' @details This function converts one hot variables into factor variables
#'
#' @keywords internal
#'

onehot_to_factor <- function(onehot_dat) {

  mask <- attr(onehot_dat, "mask")
  column_names <- attr(onehot_dat, "column_names")

  factor_dat <- onehot_dat
  factor_dat[, mask != 0] <- NULL

  for(ith_var in setdiff(unique(mask), 0)) {

    col_id <- which(setdiff(unique(mask), 0) == ith_var)

    onehot_part <- onehot_dat[, mask == ith_var]

    categories <- colnames(onehot_dat)[mask == ith_var]
    categories <- substr(categories, start = 7, stop = nchar(categories))
    categories <- sub("\\..*", "", categories)

    cat_column <- factor(apply(onehot_part, 1, function(ith_row) {
      category <- categories[which(as.logical(ith_row))]
      ifelse(length(category) == 0, NA, category)
    }), levels = as.numeric(categories))

    if(ith_var > ncol(factor_dat)) {
      factor_dat <- cbind(factor_dat, dummy_col_123_unique = cat_column)
    } else {
      factor_dat <- cbind(
        factor_dat[, 1:(ith_var - 1)],
        dummy_col_123_unique = cat_column,
        factor_dat[, (ith_var + 1):(ncol(factor_dat) + length(unique(mask)) - 1)]
      )
    }

    colnames(factor_dat)[ith_var] <- column_names[col_id]
  }
}
