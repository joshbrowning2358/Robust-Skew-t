##' Computes the robust-ified, penalized deviance for the multivariate skew-t.
##' 
##' @param param Optimization parameters derived from xi, Omega, alpha, nu.
##' For conversion, see the functions optpar2dplist and dplist2optpar.
##' @param x A matrix of the independent variables for the fit. Typically just
##' a matrix of ones.
##' @param y A matrix of dependent variables.
##' @param k A parameter controlling the robustness of the fit.  The largest
##' possible value for the negative log-likelihood is 2*k, and the negative
##' log-likelihood is adjusted down whenever it is larger than k.
##' 
##' @return The "robust" skew-t deviance evaluated at observations (x,y) and
##' with parameters dp.
##'  
##' @export
##' 

mst.pdev.robust = function(param, x, y, k=2, ...){
    nonRobust = sapply(1:NROW(y), function(i){
        sn:::mst.pdev( param, x=matrix(1), y=y[i,,drop=F], w=1 )
    })
    robust = ifelse(nonRobust>k, sapply(nonRobust, Psi, k=k), nonRobust)
    return( sum(robust) )
}