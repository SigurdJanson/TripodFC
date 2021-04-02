library(testthat)



test_that("Area50p", {
  #Area50P(X, Values)
  
  # Symmetry tests; no interpolation required
  for(i in 1:50) {
    X <- 1:(i*2)
    Data <- sample(9999, i)
    Data <- c(Data, Data[i:1])
    O <- Area50P(X, Data)
    E <- i
    expect_equal(E, O)
  }
  # Symmetry tests; with interpolation
  for(i in seq(3, 101, by = 2)) {
    X <- 1:i
    Data <- sample(9999, i %/% 2)
    Data <- c(Data, sample(99, 1), Data[length(Data):1])
    O <- Area50P(X, Data)
    E <- i / 2
    expect_equal(E, O)
  }
  
  expect_error(Area50P(1:10, 1:11), "X and Values must have the same length")
  expect_warning(Area50P(1, 1), "A distribution of length 1")
  suppressWarnings(expect_equal(1, Area50P(1, 1)))
})




test_that("Area50P: Example by Vose 1", {
  # Using example by Vose
  Precision <- 100000
  Values = matrix(c(11,12,8,9,13,13,10,10,17,16,13,15), 4, 3,
                 dimnames = list(c("Peter", "Jane", "Paul", "Susan"),
                                 c("Optimal", "Typical", "Pessimistic")))
  X <- seq(min(Values), max(Values), length.out = Precision)
  Weight <- c(0.3, 0.2, 0.4, 0.1)
  Data <- NULL
  for(Line in 1:nrow(Values)) {
    NewData <- dpert(X, min  = Values[Line, "Optimal"], 
                     mode = Values[Line, "Typical"], 
                     max  = Values[Line, "Pessimistic"])
    Data <- rbind(Data, NewData)
  }
  rownames(Data) <- rownames(Values)
  
  # 
  PooledData <- PoolOpinions(Data, Weights = Weight, Method = "Arithmetic")
  # 
  Threshold50 <- Area50P(X, PooledData)
  expect_equal(11.8, Threshold50)
})



test_that("Area50P: Example by Vose 1", {
  # Using example by Vose
  Precision <- 100000
  Values = matrix(c(120, 140, 190, 150, 180, 220), 2, 3,
                  dimnames = list(c("Peter", "Jane"),
                                  c("Optimal", "Typical", "Pessimistic")))
  X <- seq(min(Values), max(Values), length.out = Precision)
  Weight <- c(0.5, 0.5)
  Data <- NULL
  for(Line in 1:nrow(Values)) {
    NewData <- dpert(X, min  = Values[Line, "Optimal"], 
                        mode = Values[Line, "Typical"], 
                        max  = Values[Line, "Pessimistic"])
    Data <- rbind(Data, NewData)
  }
  rownames(Data) <- rownames(Values)
  
  # 
  PooledData <- PoolOpinions(Data, Weights = Weight, Method = "Arithmetic")
  # 
  Threshold50 <- Area50P(X, PooledData)
  expect_equal(169.2857, Threshold50)
  
})


test_that("Compare with qpert", {
  MaxWidth <- 100
  StepSize <- 0.0001
  
  for (i in 1:10) {
    a <- runif(1, 0, MaxWidth)
    c <- a + runif(1, 0, MaxWidth)
    b <- runif(1, a, c)
    expect <- qpert(0.5, a, b, c)
    X <- seq(a, c, StepSize)
    Data <- dpert(X, a, b, c)
    obs <- expect_silent( Area50P(X, Data) )
    expect_equal(obs, expect, tolerance = StepSize)
  }
})
