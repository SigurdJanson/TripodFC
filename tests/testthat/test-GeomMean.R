library(testthat)
#source("../R/threepoint_core.R")


test_that("All ones",{
  expect_identical(1, weighted.gmean(c(1, 1, 1), c(1, 1, 1))) # == 1
})


test_that("Single value with weight = 1 returns identical value",{
  #' Use single value
  expect_equal(12, weighted.gmean(12, 1))

  #' Random value
  for (Tests in 1:10) {
    value   <- runif(1, 0, 100)
    weights <- rep(1, 1)
    obs <- expect_silent(weighted.gmean(value, weights))
    expect_equal(obs, value)
  }
})




test_that("Zero values return 0",{
  #' Vector with increasing length
  #' But with all values = 0 (weights > 1)
  #' 
  for (Tests in 1:10) {
    value   <- rep(0L, Tests)
    weights <- runif(Tests, 0.0001, 1)
    obs <- expect_silent(weighted.gmean(value, weights))
    expect_equal(obs, 0)
  }
})



#' Vector with increasing length
#' But with all weights = 0, the result is NaN
test_that("Zero weights return NaN",{
  for (Tests in 1:10) {
    value   <- runif(Tests, 0, 100)
    weights <- rep(0L, Tests)
    obs <- expect_silent(weighted.gmean(value, weights))
    expect_true(is.nan(obs))
  }
})



test_that("Basic calculations with weight = 1",{
  #' Result: sqrt(1000)
  value   <- c(1, 1000)
  weights <- rep(1L, 2)
  obs <- expect_silent(weighted.gmean(value, weights))
  expect_equal(obs, sqrt(1000))
  #' 
  value   <- c(10, 100)
  weights <- rep(1L, 2)
  obs <- expect_silent(weighted.gmean(value, weights))
  expect_equal(obs, sqrt(1000))
  
  #' Result: 1
  for (Test in 1:10) {
    value <- runif(1L, 0.000001, 100)
    value <- c(value, 1/value)
    weights <- rep(1L, 2)
    obs <- expect_silent(weighted.gmean(value, weights))
    expect_equal(obs, 1)
  }
  
  #' Random values compared to `exp(mean(log(...)))`
  for (Test in 1L:10L) {
    value <- runif(Test, 0.000001, 100)
    weights <- rep(1L, Test)
    obs <- expect_silent(weighted.gmean(value, weights))
    expect_equal(obs, exp(mean(log(value))))
  }
  
  # 
  for(i in 2:100) {
    x <- runif(i, 2^-20, 2^10)
    expect_equal(prod(x)^(1/length(x)), weighted.gmean(x)) # no weights
  }
})


test_that("Basic calculations with all weights being equal (which is equivalent to all weigths = 1)",{
  #' Result: 1
  for (Test in 1:10) {
    value <- runif(1L, 0.000001, 100)
    value <- c(value, 1/value)
    weights <- rep(1L, 2)
    obs <- expect_silent(weighted.gmean(value, weights))
    expect_equal(obs, 1)
  }
  
  #' Random values compared to `exp(mean(log(...)))`
  for (Test in 1L:10L) {
    value <- runif(Test, 0.000001, 100)
    weights <- rep(1L, Test)
    obs <- expect_silent(weighted.gmean(value, weights))
    expect_equal(obs, exp(mean(log(value))))
  }
})



#' https://www.dummies.com/education/math/business-statistics/how-to-find-the-weighted-geometric-mean-of-a-data-set/
test_that("Specific Example 1",{
  #' 
  value <- 1:5
  weights <- c(2, 5, 6, 4, 3)
  obs <- expect_silent(weighted.gmean(value, weights))
  expect_equal(obs, 746496000^(1/20))
})



test_that("Specific Example 2",{
  # 2 and 18 = √(2 × 18) = 6; using all weights = 1
  expect_identical(6, weighted.gmean(c(2, 18), c(1, 1))) # 
  expect_identical(6, weighted.gmean(c(2, 18))) # no weights
})


test_that("Specific Example 3", {
  # Missing weights sets weights to 1

  # 3√(10 × 51,2 × 8) = 16
  expect_equal(16, weighted.gmean(c(10, 51.2, 8))) # no weights
  # 5√(1 × 3 × 9 × 27 × 81) = 9
  expect_equal(9, weighted.gmean(c(1, 3, 9, 27, 81))) # no weights
  # √(0.275 × 10^-9 × 8,8 × 10^3) ≈ 0,0016 m
  expect_equal(0.0015556349186104045, weighted.gmean(c(0.275e-9, 8.8e3))) # no weights
})


# Tests with weights (https://www.fxsolver.com/solve/)
test_that("Specific Example 4", {
  expect_equal(7.36806299728, weighted.gmean(c(10, 4), c(2, 1))) # integer weights
  expect_equal(8.72406186132, weighted.gmean(c(4, 8, 32, 64), c(4, 2, 1, 1))) # 
  
  expect_equal(9.617592736, weighted.gmean(c(10, 4), c(2.25, 0.1))) # float weights
  expect_equal(20.813830185, weighted.gmean(c(10, 4), c(2.25, -1))) # negative weights
})


  
test_that("Bad input", {
  # No negative numbers, no matter how small they may be
  expect_identical(NaN, weighted.gmean(c(2, 18, -2^-20), c(1, 1, 1)))
  # without na.rm this is NA
  expect_identical(NA_real_, weighted.gmean(c(NA, 18, 2), c(1, 1, 1)))
  # This works unlike the one above
  expect_identical(6, weighted.gmean(c(NA, 18, 2), c(1, 1, 1), na.rm = TRUE))
  # Meaningless when all weights are zero
  expect_identical(NaN, weighted.gmean(c(1, 18, 2), c(0, 0, 0), na.rm = TRUE))
})
