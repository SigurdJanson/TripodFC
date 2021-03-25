library(testthat)
#source("./Three-Point-Estimating.core.R")




test_that("Weighted Multiplicative Pooling", {
  # Tests without weights
  # Simple case
  expect_identical(1, weighted.prod(c(1, 1, 1), c(1, 1, 1))) # == 1
  expect_identical(1, weighted.prod(c(1, 1, 1))) # == 1
  expect_identical(prod(2, 18, -2^-20), weighted.prod(c(2, 18, -2^-20), c(1, 1, 1))) #
  # If there is only one 0 in x, then we get zero
  expect_identical(0, weighted.prod(c(rep(1, 25), 0), runif(26, 1, 999)))
  #
  expect_equal(prod((2:11)^(1:10)), weighted.prod(2:11, 1:10))
  
  # Bad input
  # without na.rm this is NA
  expect_identical(NA_real_, weighted.prod(c(NA, 18, 2), c(1, 1, 1)))
  # with na.rm it works again
  expect_identical(prod(1:5), weighted.prod(c(NA, 1:5), rep(1, 6), na.rm = TRUE)) # 
  # Meaningless when all weights are zero
  expect_identical(NaN, weighted.prod(c(1, 18, 2), c(0, 0, 0), na.rm = TRUE))
  expect_identical(NaN, weighted.prod(runif(25, 1, 999), rep(0, 25)))
})
