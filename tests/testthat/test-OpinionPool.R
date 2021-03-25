library(testthat)
#source("./Three-Point-Estimating.core.R")



test_that("Opinion Pooling", {
  #PoolOpinions <- function( Values, Weights = 1/nrow(Values), 
  #                          Method = c("Arithmetic", "Geometric", "Multiplicative") ) {
  # Equal weights
  for(i in 1:20) {
    MaxCol <- 200
    MaxRow <- 20
    Size <- sample(MaxCol, 1)
    Data <- runif(Size, 1, MaxCol)
    for(R in 2:MaxRow) Data <- rbind(Data, runif(Size, 1, MaxCol))
    O <- PoolOpinions(Data, Method = "Arithmetic")
    E <- apply( Data, 2, mean )
    expect_equal(E, O)
    
    O <- PoolOpinions(Data, Method = "Geometric")
    E <- apply( Data, 2, function(vec) prod(vec^(1/length(vec))) )
    expect_equal(E, O)
    
    O <- PoolOpinions(Data, Method = "Multiplicative")
    E <- apply( Data, 2, function(vec) prod(vec^(1/length(vec))) )
    expect_equal(E, O)
  }
})
