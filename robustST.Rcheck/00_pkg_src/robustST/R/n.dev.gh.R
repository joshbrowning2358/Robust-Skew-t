##' Gradient of Normal Deviance
##' 
##' This function computes the derivative of the normal deviance (where the
##' deviance is computed at observations y and with parameters dp).
##' 
##' @param dp The density parameters, in this case a vector of the mean and sd.
##' @param y A vector of observations.
##' 
##' @return A vector of length 2 giving the derivative of the deviance with
##' respect to the mean and standard deviation.
##' 
##' @export
##' 

n.dev.gh = function(dp, y){
    mu = dp[1]
    sigma = dp[2]
    return( c( -sum(y-mu)/(sigma^2), -sum((y-mu)^2/sigma^3) + length(y)/sigma ) )
}
