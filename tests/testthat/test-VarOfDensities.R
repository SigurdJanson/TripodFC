library(testthat)
library(mc2d)

test_that("Normal Distribution mean = 0, sd = 1; variable length of X", {
  DistrSd   <- 1.0
  DistrMean <- 0.0
  for (Delta in c(0.25, 0.1, 0.05, 0.01, 0.001)) {
    x <- seq(-10, 10, Delta)
    obs <- expect_silent( VarOfDensities(x, dnorm(x, DistrMean, DistrSd)Truncated = TRUE) )
    expect_equal(obs, DistrSd)
  }
})


test_that("Normal Distribution sd = 1 and variable distribution mean", {
  DistrSd   <- 1.0
  Delta <- 0.001
  for (DistrMean in c(0.0, 1.0, 5.0, 25.0, 100.0)) {
    x <- seq(-10.0, 10.0, Delta) + DistrMean
    obs <- expect_silent( VarOfDensities(x, dnorm(x, DistrMean, DistrSd), Truncated = TRUE) )
    expect_equal(obs, DistrSd)
  }
})


test_that("Normal Distribution with variable sd", {
  Delta <- 0.001
  DistrMean <- 10
  for (DistrSd in c(0.1, 0.5, 1.0, 5.0, 25.0, 100.0)) {
    x <- seq(-10.0 * DistrSd, 10.0 * DistrSd, Delta) + DistrMean
    obs <- expect_silent( VarOfDensities(x, dnorm(x, DistrMean, DistrSd), Truncated = TRUE) )
    expect_equal(obs, DistrSd^2)
  }
})



test_that("Chi² - Asymmetrisch", {
  Delta <- 0.001
  for (DoF in c(0.1, 1, 5.0, 25.0, 100.0)) {
    DistrMean <- DoF
    DistrVar  <- 2 * DoF
    x <- seq(Delta, 20.0 * max(1, DoF), Delta)
    
    obs <- expect_silent( VarOfDensities(x, dchisq(x, DoF), Truncated = TRUE) )
    expect_equal(obs, DistrVar)
  }
})



#' Truncated is FALSE because the range of PERT is limited by min/max
test_that("Symmetric Pert with Mode 0", {
  Delta <- 0.001
  DistrMean <- 0.00
  for (Min in c(-1, -10, -100)) {
    Max <- -Min
    x <- seq(Min, Max, Delta) + DistrMean
    DistrSd <- (Max-Min) / 6
    
    obs <- expect_silent( VarOfDensities(x, dpert(x, Min + DistrMean, DistrMean, Max + DistrMean), Truncated = FALSE) )
    expect_equal(obs, DistrSd^2)
  }
})


