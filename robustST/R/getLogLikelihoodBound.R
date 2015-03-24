##' Get Log Likelihood Bound
##' 
##' For the robust skew-t, we will occassionally need to understand which
##' values of the negative log-likelihood are "extreme", and this depends on
##' the current estimate of the skew-t parameters.  This function computes the
##' negative log-likelihood of the 1-alpha quantile for the skew-t distribution
##' with provided density parameters.
##' 
##' Caution: This is an approximation!  It is very accurate for most cases, but
##' it is not perfect.  In particular, when alpha is close to 1, the
##' approximation becomes fairly poor.  Fortunately, this is usually not a
##' scenario of much importance as we tend to be more interested in the case of
##' alpha close to 0 (i.e. the tails) rather than 1.
##' 
##' @param dp A named list of the density parameters of the skew-t.
##' @param alpha The quantile of the negative log-likelihood.
##' 
##' @return The value of the negative log-likelihood that corresponds to the
##' 1-alpha quantile of the skew-t with density parameters as provided in dp.
##' 

getLogLikelihoodBound = function(dp, alpha=.01){

  if(alpha == 0)
    return(-1e10)
  ### Data Quality Checks
  stopifnot(setequal(names(dp), c("xi", "Omega", "alpha", "nu")) |
            setequal(names(dp), c("beta", "Omega", "alpha", "nu")))
  stopifnot(alpha >= 0)
  stopifnot(alpha <= 1)
  stopifnot(length(alpha) == 1)
  if(names(dp)[1] == "beta")
    names(dp)[1] = "xi"
  dimension = length(dp$xi)
  if(dimension > 2)
      stop("Function currently only implemented for univariate & bivariate!")
  stopifnot(dim(dp$Omega) == c(dimension, dimension))
  stopifnot(length(dp$alpha) == dimension)
  stopifnot(length(dp$nu) == 1)
  stopifnot(sapply(dp, "class") %in% c("numeric", "matrix"))
  stopifnot(dp$nu <= Inf)
  
  ### Rename alpha, since Azzalini's code uses alpha
  pValue = alpha

  if(dimension == 1){
    return(sn::qst(p = alpha, dp = do.call("c", dp)))
  } else if(dimension == 2){
    ### Mimic Azzalini's code for cutoff levels.  See sn:::plot.SECdistrMv and
    ### sn:::plot.SECdistrBv
    nu <- dp$nu
    Omega <- dp[[2]]
    Omega.bar <- cov2cor(Omega)
    alpha <- dp[[3]]
    alpha.star <- sqrt(sum(alpha * as.vector(Omega.bar %*% alpha)))
    # Values of 0 cause problems, but setting to a small value fixes these
    # without causing other problems
    if(alpha.star == 0)
      alpha.star = .00001
    omega <- sqrt(diag(Omega))
    l.nu <- (-1.3/nu - 4.93)
    h <- 100 * log(exp(((1.005 * alpha.star - 0.045) * l.nu - 
      1.5)/alpha.star) + 1)
    K <- h * (1.005 * alpha.star - 0.1) * (1 + nu)/(alpha.star * 
      nu)
    qF <- qf(pValue, 2, nu)
    # log.levels gives the value of the (positive) log-likelihood corresponding
    # to the provided value of alpha (or pValue, as it was renamed).
    log.levels <- (lgamma(nu/2 + 1) - lgamma(nu/2) - log(pi * 
      nu) - 0.5 * log(1 - Omega.bar[1, 2]^2) - (nu/2 + 
      1) * log(2 * qF/nu + 1) + K - sum(log(omega)))
    return(-log.levels)
  }
}



### Testing this function:
# library(sn)
# library(ggplot2)
# 
# n = 1000
# dp = list(xi = c(2,1), Omega = diag(2), alpha = c(3, -5), nu = 6)
# y = rmst(n, dp = dp)
# logLikelihood = sapply(1:n, function(i){
#     sn:::mst.pdev(y = y[i, , drop = FALSE], param = sn:::dplist2optpar(dp),
#               x = matrix(1, nrow = 1)) / 2 # divide by two to go dev -> LL
# })
# 
# bounds = sapply(seq(0, 1, .05), getLogLikelihoodBound, dp = dp)
# bounds = data.frame(bounds)
# bounds$level = seq(0, 1, .05)
# qplot(logLikelihood) +
#     geom_vline(data = bounds, aes(xintercept = bounds), color = "red")
# # Determine which bucket of the bounds each likelihood value falls into
# likelihoodIndex = findInterval(logLikelihood, bounds$bounds)
# bounds$count = tapply(rep(1, n), INDEX = likelihoodIndex, FUN = sum)
# # The next plot should be roughly horizontal: same counts for each group
# ggplot(bounds, aes(x = level, y = count)) + geom_point()
# # The next plot should be a straight line: observed quantiles vs theoretical
# ggplot(bounds, aes(x = level, y = cumsum(count))) + geom_point()