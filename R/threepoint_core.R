#'
#' Basic app functions for "3-point Estimator" App
#' 
#'
library(distr)
library(mc2d)


#' @details The approach by Davis avoids assumptions regarding 
Pert.Params <- function(min=-1, mode=0, max=1, shape = 4, method = c("classic", "vose", "davis")) {
  ## check input
  if (!exists("min")) stop("'min' is missing")
  if (!exists("mode")) stop("'mode' is missing")
  if (!exists("max")) stop("'max' is missing")
  if (!exists("shape")) stop("'shape' is missing")
  if (!is.numeric(min)) stop("'min' must be a numeric value")
  if (!is.numeric(mode)) stop("'mode' must be a numeric value")
  if (!is.numeric(max)) stop("'max' must be a numeric value")
  if (!is.numeric(shape)) stop("'shape' must be a numeric value")
  
  if (!exists("method")) stop("'method' is missing")
  method <- match.arg(method)
  
  if (method == "classic") {
    mu <- (min + shape * mode + max) / (shape + 2)
    sdev <- (max - min) / (shape + 2)
    alpha <- ((mu - min) / (max - min)) * ( ((mu - min) * (max - mu) / (sdev^ 2 )) - 1 )
    beta <- alpha * (max - mu) / (mu - min)
  }
  
  if (method == "vose") {
    mu <- (min + shape * mode + max) / (shape + 2)
    alpha <- ifelse(mu == mode,
                    1 + shape / 2, # avoid div/0
                    ((mu - min) * (2 * mode - min - max)) / ((mode - mu) * (max - min))) #ok
    beta <- alpha * (max - mu) / (mu - min)
  }
  if (method == "davis") {
    alpha <- (2 * (max + 4*mode - 5*min) / (3*(max-min))) *
      (1 + 4 * ( (mode-min)*(max-mode)/(max-min)^2 ))
    beta <- alpha * (5*max-4*mode-min) / (max+4*mode-5*min)
  }
  if(is.null(alpha)) stop("Method not available")
    
  out <- list(alpha = alpha, beta = beta,
              min = min, mode = mode, max = max,
              method = method)
  class(out) <- "betaPERT"
  
  return(out)
}


#' betaPERT
#' Parametrize a generalized Beta distribution
#' @param min Pessimistic estimate (Minimum value)
#' @param mode Most likely estimate (Mode)
#' @param max Optimistic estimate (Maximum value)
#' @param shape Scale parameter
#' @param method "classic" or "vose"; see details below
#' @author Brecht Devleesschauwer <\email{brechtdv@@gmail.com}> 
betaPERT <- function(min=-1, mode=0, max=1, shape = 4, method = c("classic", "vose")) {
    ## check input
    if (!exists("min")) stop("'min' is missing")
    if (!exists("mode")) stop("'mode' is missing")
    if (!exists("max")) stop("'max' is missing")
    if (!exists("shape")) stop("'shape' is missing")
    if (!is.numeric(min)) stop("'min' must be a numeric value")
    if (!is.numeric(mode)) stop("'mode' must be a numeric value")
    if (!is.numeric(max)) stop("'max' must be a numeric value")
    if (!is.numeric(shape)) stop("'shape' must be a numeric value")
  
    if (!exists("method")) stop("'method' is missing")
    method <- match.arg(method)
    
    if (method == "classic") {
      mu <- (min + shape * mode + max) / (shape + 2)
      sdev <- (max - min) / (shape + 2)
      alpha <- ((mu - min) / (max - min)) * ( ((mu - min) * (max - mu) / (sdev^ 2 )) - 1 )
      beta <- alpha * (max - mu) / (mu - min)
    }
    
    if (method == "vose") {
      mu <- (min + shape * mode + max) / (shape + 2)
      alpha <- ifelse(mu == mode, 
                      1 + shape / 2, # avoid div/0
                      ((mu - min) * (2 * mode - min - max)) / ((mode - mu) * (max - min)))
      beta <- alpha * (max - mu) / (mu - min)
    }
    
    out <- list(alpha = alpha, beta = beta,
                min = min, mode = mode, max = max,
                method = method)
    class(out) <- "betaPERT"
    
    return(out)
  }

#' #' dpert
#' #' @description Density function for the PERT (aka Beta PERT) distribution with minimum 
#' #' equals to min, mode equals to mode and maximum equals to max.
#' #' @usage dpert(x, min=-1, mode=0, max=1, shape=4, log=FALSE)
#' #' @param x Vector of quantiles.
#' #' @param min Vector of minima.
#' #' @param mode Vector of modes.
#' #' @param max Vector of maxima.
#' #' @param shape	Vector of scaling parameters. Default value: 4.
#' #' @param log Logical; if TRUE, probabilities p are given as log(p).
#' #' @details Taken from the 'mc2d' package
#' #' @author Regis Pouillot
#' #' @author Matthew Wiener
#' #' @references R. Pouillot, M.-L. Delignette-Muller (2010), Evaluating variability and uncertainty in 
#' #' microbial quantitative risk assessment using two R packages. International Journal of Food Microbiology. 142(3):330-40 
#' dpert <- function (x, min = -1, mode = 0, max = 1, shape = 4, log = FALSE) 
#' {
#'   if (length(x) == 0) 
#'     return(numeric(0))
#'   
#'   min <- as.vector(min)
#'   mode <- as.vector(mode)
#'   max <- as.vector(max)
#'   shape <- as.vector(shape)
#' 
#'   a1 <- 1 + shape * (mode - min)/(max - min)
#'   a2 <- 1 + shape * (max - mode)/(max - min)
#'   
#'   oldw <- options(warn = -1)
#'   d <- (x - min)^(a1 - 1) * (max - x)^(a2 - 1) / 
#'     beta(a = a1, b = a2) / 
#'     (max - min)^(a1 + a2 - 1)
#'   options(warn = oldw$warn)
#'   
#'   d[x < min | x > max] <- 0
#'   d[mode < min | max < mode] <- NaN
#'   d[shape <= 0] <- NaN
#'   
#'   if (log) d <- log(d)
#'   if (any(is.na(d))) 
#'     warning("NaN in dpert")
#'   
#'   return(d)
#' }
#' 
#' 
#' #' qpert
#' #' @description Quantile function for the PERT (aka Beta PERT) distribution with minimum 
#' #' equals to min, mode equals to mode and maximum equals to max.
#' #' @usage qpert(p, min=-1, mode=0, max=1, shape=4, lower.tail=TRUE, log.p=FALSE)
#' #' @param p Vector of probabilities.
#' #' @param min Vector of minima.
#' #' @param mode Vector of modes.
#' #' @param max Vector of maxima.
#' #' @param shape	Vector of scaling parameters. Default value: 4.
#' #' @param log.p Logical; if TRUE, probabilities p are given as log(p).
#' #' @param lower.tail Logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' #' @details Taken from the 'mc2d' package
#' #' @author Regis Pouillot
#' #' @author Matthew Wiener
#' #' @references R. Pouillot, M.-L. Delignette-Muller (2010), Evaluating variability and uncertainty in 
#' #' microbial quantitative risk assessment using two R packages. International Journal of Food Microbiology. 142(3):330-40 
#' qpert <- function (p, min = -1, mode = 0, max = 1, shape = 4, lower.tail = TRUE, 
#'           log.p = FALSE) 
#' {
#'   if (length(p) == 0)
#'     return(numeric(0))
#'   min <- as.vector(min)
#'   mode <- as.vector(mode)
#'   max <- as.vector(max)
#'   shape <- as.vector(shape)
#' 
#'   lout <- max(length(p), length(min), length(mode), length(max), length(shape))
#'   min <- rep(min, length.out = lout)
#'   mode <- rep(mode, length.out = lout)
#'   max <- rep(max, length.out = lout)
#'   shape <- rep(shape, length.out = lout)
#' 
#'   if (log.p) 
#'     p <- exp(p)
#'   if (!lower.tail) 
#'     p <- 1 - p
#' 
#'   a1 <- 1 + shape * (mode - min)/(max - min)
#'   a2 <- 1 + shape * (max - mode)/(max - min)
#'   
#'   oldw <- options(warn = -1)
#'   q <- qbeta(p, shape1 = a1, shape2 = a2)
#'   options(warn = oldw$warn)
#'   
#'   q <- q * (max - min) + min
#'   minmodemax <- (abs(min - max) < (.Machine$double.eps^0.5))
#'   q <- ifelse(minmodemax, min, q)
#'   q[p < 0 | p > 1] <- NaN
#'   q[mode < min | max < mode] <- NaN
#'   q[shape <= 0] <- NaN
#'   
#'   if (any(is.na(q))) 
#'     warning("NaN in qpert")
#'   return(q)
#' }



#' weighted.gmean
#' @param x an object containing the values whose weighted mean is to be computed.
#' @param w	a numerical vector of weights the same length as x giving the weights 
#' to use for elements of x.
#' @param na.rm	a logical value indicating whether NA values in x should be 
#' stripped before the computation proceeds.
#' @details Returns NaN when all weights are 0.
#' @value A length-one numeric vector.
weighted.gmean <- function(x, w, scl = 1, na.rm = FALSE) {
  # PRECONDITIONS
  if(any(x < 0, na.rm = TRUE)) return(NaN)
  if(missing(w)) w <- rep(1, length(x))
  if(all(w == 0, na.rm = TRUE)) return(NaN)
  if(length(x) != length(w)) stop("Weights w must have the same length as x")
  if(na.rm) {
    w <- w[!is.na(x)]
    if(length(w) == 0) return(NaN)
    x <- x[!is.na(x)]
  }
  
  # CODE
  num <- sum(w * log(x)) # numerator
  den <- sum(w) # denominator
  return(scl * exp(num/den))
}



#' weighted.prod
#' @param x an object containing the values whose weighted mean is to be computed.
#' @param w	a numerical vector of weights the same length as x giving the weights 
#' to use for elements of x.
#' @param na.rm	a logical value indicating whether NA values in x should be 
#' stripped before the computation proceeds.
#' @note This function is not finished/correct.
#' @details Returns NaN when all weights are 0.
#' @value A length-one numeric vector.
weighted.prod <- function(x, w, scl = 1, na.rm = FALSE) {
  # PRECONDITIONS
  if(missing(w)) w <- rep(1, length(x))
  if(all(w == 0, na.rm = TRUE)) return(NaN)
  if(length(x) != length(w)) stop("Weights w must have the same length as x")

  # CODE
  return(scl * prod(x^w, na.rm))
}


#' PoolOpinions
#' Pools opinions according to the rationale laid out by Dietrich & List (2016).
#' @param Values A matrix with Pert distributions as rows
#' @param Weights Default: unweighted pooling
#' @param Method Arithmetic, geometric or multiplicative pooling. Default is 'Arithmetic'.
#' @details Instead of more complex pooling methods, this approach sticks to
#' the simple solutions that perform well, according to Clemen & Winkler (1999).
#' Weights must be non-negative (Clemen & Winkler, 1999). They will normalized be 
#' \code{PoolOpinions} so that they sum up to 1.
#' @references Dietrich, F. & List, C. (2016). Probabilistic Opinion Pooling. 
#' In A. Hajek & C. Hitchcock (eds.), Oxford Handbook of Philosophy and Probability. 
#' Oxford: Oxford University Press
#' Clemen, R.T. & Winkler R.L. (1999). Combining Probability Distributions 
#' From Experts in Risk Analysis. Risk Analysis 19, 187–203. doi:10.1023/A:1006917509560
PoolOpinions <- function( Values, Weights = rep(1/nrow(Values), nrow(Values)), 
                          Method = c("Arithmetic", "Geometric", "Multiplicative") ) {
  if(!is.matrix(Values)) { #try to coerce
    Values <- as.matrix(Values)
  }
  if(nrow(Values) < 2) {
    warning("Parameter 'Values' has only 1 row. A single vector cannot be pooled with itself.")
    return(Values)
  }
  if(any(Weights < 0)) stop("Weights must be > 1")
  if(sum(Weights != 1)) { #normalize weights if required
    Weights <- Weights / sum(Weights)
  }
  
  Method <- match.arg( Method )
  switch(Method,
         Arithmetic     = apply(Values, 2, weighted.mean,  Weights),
         Geometric      = apply(Values, 2, weighted.gmean, Weights), 
         Multiplicative = apply(Values, 2, weighted.prod,  Weights),
         apply(Values, 2, weighted.mean,  Weights)
         )
}


#' Area50P
#' Computes the average of the integral function based on empirical data.
#' @param X A numeric vector
#' @param dX 'histogram' data, i.e. observations based on some probability distribution
#' (\code{dX = f(X)} with \code{f} being an arbitrary density function).
#' @return The value of \code{X} at which the integral of the probability distribution 
#' exceeds 50% of it's area.
Area50P <- function(X, dX) {
  if(length(X) != length(dX)) stop("X and Values must have the same length")
  if(length(X) == 1) {
    warning("A distribution of length 1 cannot be divided in half")
    return(1)
  }
  Percentages <- dX / sum(dX)
  Accumulated <- cumsum(Percentages)
  
  # if there is an exact match for 0.5 take this exact index ...
  Index50 <- which(Accumulated == 0.5)
  # ... otherwise interpolate
  if(identical(Index50, integer(0))) {
    Index50 <- min(which(Accumulated > 0.5))  # get index where value goes over .5
    
    # Get an interpolated value: BaseX(<50) + (DeltaOfX(<50, >50) * Weight)
    Threshold50 <- X[Index50-1] + (X[Index50] - X[Index50-1]) * 
      ((0.5 - Accumulated[Index50-1]) / (Accumulated[Index50] - Accumulated[Index50-1]))
  } else {
    Threshold50 <- X[Index50] # convert index to x-value
  }
  
  return(Threshold50)
}


#' .dX2pX 
#'
#' @param X 
#' @param dX Probability densities of X
#' @param Truncated Does the range cover the whole distribution or only a part of it?
#' @return Gives the distribution mass derived from the densities.
.dX2pX <- function(X, dX, Truncated = FALSE) {
  Deltas <- range(diff(X)) # Intervals of X
  #if (diff(Deltas) > .Machine$double.eps/10) stop("Equal steps among X values is required")
  Deltas <- Deltas[1]
  if (Truncated)
    Deltas <- rep(Deltas, length(X))
  else
    Deltas <- c(0.5 * Deltas, rep(Deltas, length(X)-2), 0.5 * Deltas)
  
  dX * Deltas
}



#' Estimate distribution moments from a representation of the densities.
#' @description Computes the mean of a probability distribution from a series 
#' of density values.
#' @inheritParams .dX2pX
#' @note The correctness and precision depends on the quality of the data and
#' the shape of the probability distribution. Large quantile steps cannot lead 
#' to a precise estimate. Also problematic is the behaviour of the probability
#' distribution at it's outer edges. Optimal are probability distributions like
#' the beta Pert distribution which ends on both sides of its range with a density 
#' of zero. More difficult are distributions that converge to zero with the quantile
#' approaching infinity (like the normal distribution). Most problematic are 
#' functions that converge to infinity.
#' @export
AvgOfDensities <- function(X, dX, Truncated = FALSE) {
  pX <- .dX2pX(X, dX, Truncated)
  sum(X * pX)
}



#' @describeIn AvgOfDensities Variance of a probability distribution 
#' @export
VarOfDensities <- function(X, dX, Truncated = FALSE) {
  pX <- .dX2pX(X, dX, Truncated)
  Avg <- sum(X * pX)
  SumOfSquares <- sum(mapply(function(x, y) (x^2*y), X, pX))
  SumOfSquares - Avg^2
}

