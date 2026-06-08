
library(testthat)

set.seed(10)

X <- matrix(rnorm(100), nrow = 25)
X_imp <- matrix(rnorm(100), nrow = 25)

test_that("energy distance works", {

  expect_identical(round(edistance(X, X_imp), 4),
                   round(c(`E-statistic` = 3.4158), 4))

  expect_identical(round(edistance(X, X_imp, rescale = TRUE), 4),
                   round(c(`E-statistic` = 0.2733), 4))

})






