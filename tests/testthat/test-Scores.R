set.seed(111)

X <- random_mcar_data(50, 3, 0.2)

methods_list <- list(
  exp = exp_imputation,
  norm = norm_imputation
)

test_that("We can compare IScores", {

  set.seed(123)

  res <- compare_Iscores(X, methods_list = methods_list)

  expect_s3_class(res, "data.frame")

  expect_equal(nrow(res), 4)
  expect_equal(ncol(res), 3)

  expect_true(all(is.finite(res[, 1])))

  expect_equal(as.numeric(res[, 1]),
               c(0.56466, 0.45178, 1.29138, 2.14579),
               tolerance = 0.25)

  expect_identical(as.character(res[, 2]), c("energy_IScore",
                                             "energy_IScore",
                                             "DR_IScore",
                                             "DR_IScore"))

  expect_identical(as.character(res[, 3]), c("exp", "norm", "exp", "norm"))

  expect_true(res[1, 1] > res[2, 1])
  expect_true(res[3, 1] > res[4, 1])
})
