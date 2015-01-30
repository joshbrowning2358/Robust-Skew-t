##' Gradient of Normal Deviance, Adjusted for Robustness
##' 
##' This function is the analog of n.dev.gh, but is applied with a robust
##' adjustment.  It takes advantage of the fact that the gradient of the robust
##' deviance for all observations y is the sum of the gradients of the
##' deviances for each individual observation.
##' 
##' @param dp The density parameters, in this case a vector of the mean and sd.
##' @param y A vector of observations.
##' @param k The robust adjustment parameter.
##' 
##' @return A vector of length 2 giving the derivative of the (robust) deviance
##' with respect to the mean and standard deviation.
##' 
##' @export
##' 

n.dev.gh.robust = function(dp, y, k=2){
    nonRobustGH = lapply(y, n.dev.gh, dp=dp)
    nonRobustGH = do.call("rbind", nonRobustGH)
    nonRobust = lapply(y, n.dev, dp=dp)
    nonRobust = do.call("c", nonRobust)
    # If -LogLikelihood<k, we don't adjust at all.  Thus, gradient doesn't change for that case.
    # If -LogLikelihood>k, we change -LL to Psi(-LL).  Thus, gradient becomes Psi'(-LL)*(-LL)'=psi(-LL)*(-LL)'
    robustGH = nonRobustGH
    if(any(nonRobust>k))
        robustGH[nonRobust>k,] = robustGH[nonRobust>k,]*sapply(nonRobust[nonRobust>k], psi.grad, k=k)
    return( colSums(robustGH) )
}