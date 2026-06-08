

test_that("Energy-I-Score works", {

  set.seed(111)
  X <- random_mcar_data(100, 4)
  imputation_func <- norm_imputation

  res <- energy_IScore(X, imputation_func)

  expect_equal(as.vector(round(res, 7)), 0.5785145, tolerance = 0.05)

  tbl <- attr(res, "dat")

  expect_identical(tbl,
                   structure(
                     list(
                       column_id = c(1L, 2L, 4L, 3L),
                       weight = c(0.1875, 0.1824, 0.16, 0.1476),
                       score = c(0.618580450631919, 0.587221424840952,
                                 0.566960886456749, 0.529382327748725),
                       n_columns_used = c(1L, 1L, 1L, 1L)),
                     row.names = c("X1", "X2", "X4", "X3"),
                     class = "data.frame"),
                   tolerance = 10^(-10))

  expect_warning(energy_IScore(X, imputation_func, max_length = 1000, silent = FALSE),
                 "max_length is larger than the total number of columns")
})



test_that("Energy-I-Score throws a proper error", {

  set.seed(111)
  X <- random_mcar_data(100, 4)

  expect_error(energy_IScore(X, imputation_func = 5),
               "Imputation_func must be a function!")

  imputation_func <- stop
  expect_error(energy_IScore(X, imputation_func = imputation_func),
               "Errored imputing X using provided imputation_func!")

  unstable_imp <- function(missdf) {
    if(sample(1:2, 1) == 1) {
      stop()
    } else {
      norm_imputation(missdf)
    }
  }

  expect_warning(energy_IScore(X,
                               imputation_func = unstable_imp,
                               X_imp = norm_imputation(X),
                               silent = FALSE),
               "Unsuccessful imputation! Imputation function is unstable!")

})


test_that("Energy-I-Score works for mixed data", {

  set.seed(123)
  X <-  random_mcar_mixed_data(100, 4, 2)
  imputation_func <- median_mode_imputation

  expect_message(energy_IScore(X, imputation_func, silent = FALSE),
                 "Factor variables detected.")


  res <- energy_IScore(X, imputation_func)

  expect_equal(as.vector(round(res, 7)), 0.822, tolerance = 0.05)

  tbl <- attr(res, "dat")

  expect_identical(tbl,
                   structure(list(
                     column_id = c(2L, 1L, 4L, 5L, 3L, 6L),
                     weight = c(0.1824, 0.16, 0.1539, 0.1539, 0.1476, 0.1476),
                     score = c(0.68651110663701, 0.708839930029851, 0.868062359850176,
                               0.96974644277, 0.684496594623652, 1.04528828523229),
                     n_columns_used = c(1L, 1L, 1L, 1L, 1L, 1L)),
                     row.names = c("col2", "col1", "col4", "col5", "col3", "col6"),
                     class = "data.frame"),
                   tolerance = 0.001)

})





test_that("Energy-I-Score cat throws a warning", {

  set.seed(123)
  X <-  random_mcar_mixed_data(100, 4, 2)
  imputation_func <- median_mode_imputation


  expect_warning(res <- energy_IScore(X, imputation_func, skip_if_needed = FALSE,
                                      silent = FALSE),
                 "Oj was empty. There was no complete column for training.")

  expect_true(all(is.na(attr(res, "dat")[["score"]])))



})





