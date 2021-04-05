#'
#' Basic app functions for "3-point Estimator" App
#' 
library(distr)
library(mc2d)


#' Pert2BetaParams
#' 
#' Determine the parameters of the beta distribution based on parameters of the 
#' Pert distribution. 
#' @param min The optimal estimate and minimum of the Pert distribution (double)
#' @param mode The most likely estimate and the mode of the Pert distribution (double)
#' @param max The pessimistic estimate and maximum of the Pert distribution (double)
#' @param shape The shape argument (double, default = 4). Ignored for 
#' the `davis` method.
#' @param method String indicating the estimation technique. Either 
#' "classic", "golgin", "vose" or "davis" (which is the default). See details.
#' @return A list (with class attribute "betaPERT") containing the elements
#' alpha, beta, min, mode, max, shape, and method (all doubles except method which 
#' is a string; see the arguments of this function).
#' @details 
#' The Beta-PERT methodology was developed in the context of Program Evaluation 
#' and Review Technique (PERT). Based on a pessimistic estimate (minimum value), 
#' a most likely estimate (mode), and an optimistic estimate (maximum value), 
#' typically derived through expert elicitation, the parameters of a 
#' Beta distribution can be calculated. The Beta-PERT distribution is used 
#' in stochastic modeling and risk assessment studies to reflect uncertainty 
#' regarding specific parameters.
#' 
#' Different methods exist in literature for defining the parameters of a 
#' Beta distribution based on PERT. The two most common methods are included 
#' in the `Pert2BetaParams` function:
#' 
#' Classic: The standard formulas for mean, standard deviation, α and β, are 
#' as follows:
#'    mean = (a + k*b + c) / (k + 2)
#'    sd = (c - a) / (k + 2)
#'    α = { (mean - a) / (c - a) } * { (mean - a) * (c - mean) / sd^{2} - 1 }
#'    β = α * (c - mean) / (mean - a)
#'
#' The resulting distribution is a 4-parameter Beta distribution: Beta(α, β, a, b).
#' 
#' Golenko-Ginzburg (1988; cited by Pleguezuelo, 2003) proposed this parametrization:
#' α = 1 + k * (b − a) / (c − a)
#' β = 1 + k * (c − b) / (c − a)
#' 
#' Vose: Vose (2008) describes a different formula for α:
#'    α = (mean - a) * (2 * m - a - b) / { (m - mean) * (b - a) }
#' 
#' Mean and β are calculated using the standard formulas; as for the 
#' classical PERT, the resulting distribution is a 4-parameter Beta 
#' distribution: Beta(α, β, a, b).
#'
#' Note: If m = mean, α is calculated as 1 + k/2, in accordance with the 
#' mc2d package (see 'Note').
#' 
#' Davis: the approach by Davis avoids unnecessary assumptions regarding mean
#' and standard deviation but actually makes sure that these assume the 
#' defined values (see the classic definition above).
#' 
#'   α = (2 * (b + 4*mode - 5*a) / (3*(b-a))) * 
#'       (1 + 4 * ( (mode-a)*(b-mode)/(b-a)^2 ))
#'   β = α * (5*b-4*mode-a) / (b+4*mode-5*a)
#' @references 
#' Davis, R. (2008). Teaching Note - Teaching Project Simulation in Excel 
#'  Using PERT-BetaDistributions. INFORMS Transactions on Education, 8(3), 
#'  139–148. https://doi.org/10.1287/ited.1080.0013
#' 
#' Malcolm, D.G., Roseboom, .J.H, Clark, C.E. & Fazar, W. (1959). 
#'  Application of a technique for research and development program evaluation. 
#'  Oper Res 7(5):646-669.
#'  
#' Vose, D. (2008). Risk analysis, a quantitative guide, 2nd edition. 
#'  Wiley and Sons.
#' @source This is an extension of the \link[prevalence]{betaPERT} function 
#' from the prevalence package.
Pert2BetaParams <- function(min=-1, mode=0, max=1, shape = 4, 
                            method = c("classic", "golgin", "vose", "davis")) {
  # PRECONDITIONS
  if (!is.numeric(min)) stop("Argument 'min' must be a numeric value")
  if (!is.numeric(mode)) stop("Argument 'mode' must be a numeric value")
  if (!is.numeric(max)) stop("Argument 'max' must be a numeric value")
  if (!is.numeric(shape)) stop("Argument 'shape' must be a numeric value")
  if (min > mode || mode > max) 
    stop("The order of argument values must be 'min' < 'mode' < 'max'")
  
  if (missing(method)) method <- "davis"
  method <- match.arg(method)
  
  if (method == "classic") {
    mu <- (min + shape * mode + max) / (shape + 2)
    sdev <- (max - min) / (shape + 2)
    alpha <- ((mu - min) / (max - min)) * ( ((mu - min) * (max - mu) / (sdev^2 )) - 1 )
    beta <- alpha * (max - mu) / (mu - min)
  }

  if (method == "golgin") {
    alpha <- 1 + shape * (mode - min) / (max - min)
    beta  <- 1 + shape * (max - mode) / (max - min)
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
    #beta <- alpha * (5*max-4*mode-min) / (max+4*mode-5*min)
    beta <- (2*(5*max-4*mode-min)/3/(max-min)) * 
      (1 + 4 * ((mode-min)*(max-mode)/(max-min)^2))
  }

  if(is.null(alpha)) stop("Method not available")
    
  out <- list(alpha = alpha, beta = beta,
              min = min, mode = mode, max = max, shape = shape,
              method = method)
  class(out) <- "betaPERT"
  
  return(out)
}





#' dBetaPert
#' @description Density function for the PERT (aka Beta PERT) distribution with minimum
#' equals to min, mode equals to mode and maximum equals to max.
#' @usage dBetaPert(x, min=-1, mode=0, max=1, shape=4, log=FALSE, method="davis")
#' @param x Vector of quantiles.
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @inheritParams Pert2BetaParams
#' @note Adapted from the 'mc2d' package
#' @author Regis Pouillot
#' @author Matthew Wiener
#' @references R. Pouillot, M.-L. Delignette-Muller (2010), Evaluating variability and uncertainty in
#' microbial quantitative risk assessment using two R packages. International Journal of Food Microbiology. 142(3):330-40
dBetaPert <- function(x, min = 0, mode = 0.5, max = 1, shape = 4, log = FALSE,
                      method=c("classic", "golgin", "vose", "davis")) {
  if (length(x) == 0) return(numeric(0))
  
  min   <- as.vector(min)
  mode  <- as.vector(mode)
  max   <- as.vector(max)
  shape <- as.vector(shape)
  
  Params <- Pert2BetaParams(min, mode, max, shape, method)
  alpha  <- Params$alpha
  beta   <- Params$beta
  
  oldw <- options(warn = -1)
  d <- (x - min)^(alpha - 1) * (max - x)^(beta - 1) /
    beta(a = alpha, b = beta) /
    (max - min)^(alpha + beta - 1)
  options(warn = oldw$warn)
  
  d[x < min | x > max] <- 0
  d[mode < min | max < mode] <- NaN
  d[shape <= 0] <- NaN
  
  if (log) d <- log(d)
  if (any(is.na(d)))
    warning("NaN in dBetaPert")
  
  return(d)
}



#' qBetaPert
#' @description Quantile function for the PERT (aka Beta PERT) distribution with minimum
#' equals to min, mode equals to mode and maximum equals to max.
#' @usage qBetaPert(p, min=-1, mode=0, max=1, shape=4, lower.tail=TRUE, log.p=FALSE, method="davis")
#' @param p Vector of probabilities.
#' @param log.p Logical; if TRUE, probabilities p are given as log(p).
#' @param lower.tail Logical; if TRUE (default), probabilities are P[X <= x], otherwise, P[X > x].
#' @inheritParams Pert2BetaParams
#' @note Adapted from the 'mc2d' package
#' @author Regis Pouillot
#' @author Matthew Wiener
#' @references R. Pouillot, M.-L. Delignette-Muller (2010), Evaluating variability and uncertainty in
#' microbial quantitative risk assessment using two R packages. International Journal of Food Microbiology. 142(3):330-40
qBetaPert <- function (p, min = -1, mode = 0, max = 1, shape = 4, 
                       lower.tail = TRUE, log.p = FALSE,
                       method=c("classic", "golgin", "vose", "davis"))
{
  if (length(p) == 0)
    return(numeric(0))
  min   <- as.vector(min)
  mode  <- as.vector(mode)
  max   <- as.vector(max)
  shape <- as.vector(shape)

  lout  <- max(length(p), length(min), length(mode), length(max), length(shape))
  min   <- rep(min, length.out = lout)
  mode  <- rep(mode, length.out = lout)
  max   <- rep(max, length.out = lout)
  shape <- rep(shape, length.out = lout)

  if (log.p)
    p <- exp(p)
  if (!lower.tail)
    p <- 1 - p

  Params <- Pert2BetaParams(min, mode, max, shape, method)
  alpha  <- Params$alpha
  beta   <- Params$beta

  oldw <- options(warn = -1)
  q <- qbeta(p, shape1 = alpha, shape2 = beta)
  options(warn = oldw$warn)

  q <- q * (max - min) + min
  minmodemax <- (abs(min - max) < (.Machine$double.eps^0.5))
  q <- ifelse(minmodemax, min, q)
  q[p < 0 | p > 1] <- NaN
  q[mode < min | max < mode] <- NaN
  q[shape <= 0] <- NaN

  if (any(is.na(q)))
    warning("NaN in qBetaPert")
  return(q)
}



#' weighted.gmean
#' Weighted gemotric mean
#' @param x an object containing the values whose weighted mean is to be computed.
#' @param w	a numerical vector of weights the same length as x giving the weights 
#' to use for elements of x.
#' @param scl arbitrary scaling factor
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
#' Computes a weighted product
#' @param x an object containing the values whose weighted product is to be computed.
#' @param w	a numerical vector of weights the same length as x giving the weights 
#' to use for elements of x.
#' @param scl arbitrary scaling factor
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
#' @note Precision mostly depends on the step size (X[n+1]-X[n]).
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
#' @note Precision mostly depends on the step size (X[n+1]-X[n]).
AvgOfDensities <- function(X, dX, Truncated = FALSE) {
  pX <- .dX2pX(X, dX, Truncated)
  sum(X * pX)
}



#' @describeIn AvgOfDensities Variance of a probability distribution 
#' @export
#' @note Precision mostly depends on the step size (X[n+1]-X[n]).
VarOfDensities <- function(X, dX, Truncated = FALSE) {
  pX <- .dX2pX(X, dX, Truncated)
  Avg <- sum(X * pX)
  SumOfSquares <- sum(mapply(function(x, y) (x^2*y), X, pX))
  SumOfSquares - Avg^2
}

