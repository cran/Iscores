set.seed(111)
X <- Iscores:::random_mcar_data(100, 3, 0.2)
imputation_func <- norm_imputation

test_that("DR I Score works with imputation function", {

  set.seed(123)

  res <- DR_IScore(X, imputation_func)

  expect_type(res, "double")
  expect_length(res, 1)
  expect_true(is.finite(res))
})

test_that("DR I Score works with precomputed imputations", {

  set.seed(1)

  X_imp <- lapply(1:5, function(i) {
    exp_imputation(X)
  })

  set.seed(123)

  res <- DR_IScore(X, X_imp = X_imp)

  expect_type(res, "double")
  expect_length(res, 1)
  expect_true(is.finite(res))
})

test_that("DR I Score works on small datasets", {

  set.seed(111)

  X_small <- random_mcar_data(20, 4, 0.2)

  res <- DR_IScore(X_small, norm_imputation)

  expect_type(res, "double")
  expect_length(res, 1)
  expect_true(is.finite(res))
})

test_that("DR I Score works with low missingness", {

  set.seed(111)

  X_low_miss <-random_mcar_data(20, 6, 0.1)
  X_low_miss[2, 1] <- 0.1

  res <- DR_IScore(X_low_miss, norm_imputation)

  expect_type(res, "double")
  expect_length(res, 1)
  expect_true(is.finite(res))
})

test_that("DR I Score throws an error when there is no imputed data and no imputation function", {
  expect_error(
    DR_IScore(X),
    "You must provide one of imputation_func or X_imp!"
  )
})

test_that("DR I Score throws an error when imputation function fails", {
  imp_fun <- function(...) stop("imputation failed")

  expect_error(
    DR_IScore(X, imp_fun),
    "Errored imputing X using provided imputation_func!"
  )
})
