library(testthat)

test_that("Normal Distribution mean = 0, sd = 1; variable length of X", {
  for (Delta in c(0.25, 0.1, 0.05, 0.01, 0.001)) {
    x <- seq(-5, 5, Delta)
    obs <- expect_silent( AvgOfDensities(x, dnorm(x)) )
    expect_equal(obs, 0.00)
  }
})


test_that("Normal Distribution sd = 1 and variable distribution mean", {
  Delta <- 0.001
  for (DistrMean in c(0.0, 1.0, 5.0, 25.0, 100.0)) {
    x <- seq(-10.0, 10.0, Delta) + DistrMean
    obs <- expect_silent( AvgOfDensities(x, dnorm(x, DistrMean), TRUE) )
    expect_equal(obs, DistrMean)
  }
})


test_that("Normal Distribution with variable sd = 1", {
  Delta <- 0.001
  DistrMean <- 10
  for (DistrSd in c(0.1, 0.5, 1.0, 5.0, 25.0, 100.0)) {
    x <- seq(-10.0 * DistrSd, 10.0 * DistrSd, Delta) + DistrMean
    obs <- expect_silent( AvgOfDensities(x, dnorm(x, DistrMean, DistrSd), TRUE) )
    expect_equal(obs, DistrMean)
  }
})


#' Lower tolerance required to pass
test_that("Chi² - Asymmetric", {
  Delta <- 0.0001
  for (DoF in c(0.1, 1, 5.0, 25.0, 100.0)) {
    DistrMean <- DoF
    DistrVar  <- 2 * DoF
    x <- seq(.Machine$double.eps, 25.0 * max(1, DoF) + .Machine$double.eps, Delta)
    
    obs <- expect_silent( AvgOfDensities(x, dchisq(x, DoF), Truncated = TRUE) )
    expect_equal(obs, DistrMean, tolerance = 5E-5)
  }
})


#' Truncated is FALSE because the range of PERT is limited by min/max
test_that("Symmetric Pert with Mode 0", {
  Delta <- 0.001
  DistrMean <- 9.0
  for (Min in c(-1, -10, -100)) {
    Max <- -Min
    x <- seq(Min, Max, Delta) + DistrMean
    
    obs <- expect_silent( 
      AvgOfDensities(x, dpert(x, Min + DistrMean, DistrMean, Max + DistrMean), Truncated = FALSE) 
    )
    expect_equal(obs, DistrMean)
  }
})
