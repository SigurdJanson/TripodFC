library(testthat)

# Pert2BetaParams(min=-1, mode=0, max=1, shape = 4, method = c("classic", "vose", "davis"))

# GENERAL TESTS ------
test_that("Do all methods return proper types", {
  for (mt in c("classic", "vose", "davis")) {
    min  <- runif(1, -100, 100)
    max  <- min + runif(1, 0, 100)
    mode <- runif(1, min, max)
    obs <- expect_silent(
      Pert2BetaParams(min, mode, max, shape = 4, mt)
    )
    expect_type(obs, "list")
    expect_s3_class(obs, "betaPERT")
    expect_named(obs, c("alpha", "beta", "min", "mode", "max", "shape", "method"))
    
    expect_identical(obs$min,  min)
    expect_identical(obs$mode, mode)
    expect_identical(obs$max,  max)
    expect_identical(obs$shape, shape)
    expect_identical(obs$method, mt)
  }
})


# CLASSIC


# GOLENKO-GINZBURG  -------

# example taken from 
# https://www.real-statistics.com/binomial-and-related-distributions/pert-distribution/
test_that("Specific examples", {
  min <- 1
  max <- 9
  mode <- 6
  obs <- expect_silent(
    Pert2BetaParams(min, mode, max, method = "golgin")
  )
  expect_equal(obs$alpha, 3.5)
  expect_equal(obs$beta,  2.5)
})



# VOSE -------






# DAVIS -------

#' In the symmetric case, alpha and beta must be 4
test_that("Symmetric Case", {
  for (scale in c(1, 10, 100, 1000)) {
    min <- ((scale-1) - 1) * scale # c(-1, 80, 9800, 998000)
    max <- min + scale
    mode <- (max+min) / 2
    obs <- expect_silent(
      Pert2BetaParams(min, mode, max, method = "davis")
    )
    expect_equal(obs$alpha, 4)
    expect_equal(obs$beta,  4)
  }
})


#' Relationship between alpha and beta
test_that("mode = min: then alpha + beta = 4", {
  for (min in c(-1000, -500, -50, 0, 50, 500, 1000))
    for(width in c(1, 10, 100, 1000)) {
      max <- min + width
      mode <- min
      obs <- expect_silent(
        Pert2BetaParams(min, mode, max, method = "davis")
      )
      expect_equal(obs$alpha + obs$beta, 4)
    }
})


#' Relationship between alpha and beta
test_that("mode = max: then alpha + beta = 4", {
  for (min in c(-1000, -500, -50, 0, 50, 500, 1000))
    for(width in c(1, 10, 100, 1000)) {
      max <- min + width
      mode <- max
      obs <- expect_silent(
        Pert2BetaParams(min, mode, max, method = "davis")
      )
      expect_equal(obs$alpha + obs$beta, 4)
    }
})


#' Relationship between alpha and beta
test_that("alpha + beta <= 8", {
  for (min in c(-1000, -500, -50, 0, 50, 500, 1000))
    for(width in c(1, 10, 100, 1000))
      for (i in c(0.01, 0.1, 0.25, 0.49, 0.51, 0.75, 0.9, 0.99)) {
        max <- min + width
        mode <- min + width * i
        obs <- expect_silent(
          Pert2BetaParams(min, mode, max, method = "davis")
        )
        expect_gte(obs$alpha + obs$beta, 4)
        expect_lte(obs$alpha + obs$beta, 8)
      }
})

# example taken from Davis (2008)
test_that("Specific examples", {
  min   <- c(3, 2, 2, 1, 0,  1,  3,  1,  5, 1)
  mode  <- c(6, 5, 3, 3, 7,  2,  4,  2, 10, 3)
  max   <- c(9, 6, 7, 3, 8, 10, 12, 15, 30, 4)
  alpha <- c(4, 14/3, 1.968, 10/3, 4.3125, 1.3434, 1.3434, 1.0845, 1.968, 4.6173)
  beta  <- c(4,  7/3, 4.592,  2/3, 1.4375, 4.2369, 4.2369, 3.9767, 4.592, 2.9383)
  # use the maximum possible tolerance
  tol   <- c(1E-5, 1E-15, 1E-5, 1E-15, 1E-5, 5E-5, 5E-5, 5E-5, 1E-5, 1E-5)
  
  for (i in 1:length(min)) {
    obs <- expect_silent(
      Pert2BetaParams(min[i], mode[i], max[i], method = "davis")
    )
    expect_equal(obs$alpha, alpha[i], tolerance = tol[i])
    expect_equal(obs$beta,  beta[i],  tolerance = tol[i])
  }
})

