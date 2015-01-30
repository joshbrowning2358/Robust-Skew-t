##' Robust Normal deviance
##' 
##' This function computes the deviance for a set of observations (y), a
##' provided mean and standard deviation
##' 
##' @param dp The "density parameters".  More specifically, just a vector of
##' length two containing the mean and standard deviation parameters.
##' @param y The observed value.
##' @param k The parameter controlling the robustness adjustment.  See ?Psi.
##' 
##' @return The sum of the negative log-likelihoods of a normal with parameters
##' dp evaluated at the values y.  However, the individual negative
##' log-likelihood values are adjusted by Psi.
##' 
##' @export
##' 

n.dev.robust = function(dp, y, k=2){
    nonRobust = lapply(y, n.dev, dp=dp)
    nonRobust = do.call("c", nonRobust)
    robust = ifelse(nonRobust>k, sapply(nonRobust, Psi, k=k), nonRobust)
    return( sum(robust) )
}
