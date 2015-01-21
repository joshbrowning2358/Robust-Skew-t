##' Normal deviance
##' 
##' @param dp The "density parameters".  More specifically, just a vector of
##' length two containing the mean and standard deviation parameters.
##' @param y A vector of observed values.
##' 
##' @return The sum of the negative log-likelihoods of a normal with parameters
##' dp evaluated at the values y.
##' 
##' @export
##' 

n.dev = function(dp, y){
    mu = dp[1]
    sigma = dp[2]
    return( sum( (y-mu)^2/(2*sigma^2) ) + length(y)*log(sigma) )
}