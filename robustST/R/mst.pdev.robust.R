##' Computes the robust-ified, penalized deviance for the multivariate skew-t.
##' 
##' @param param Optimization parameters derived from xi, Omega, alpha, nu.
##' For conversion, see the functions optpar2dplist and dplist2optpar.
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

mst.pdev.robust = function(param, y, k = 2, ...){

    ## Data Quality Checks
    d = NCOL(y)
    paramDim = dimensionFromParamSize(length(param))
    if(d != paramDim)
        stop("parameter suggests a different dimension of y than what is ",
             "observed.")
    
    dp = sn:::optpar2dplist(param = param, p = 1, d = d)$dp
    nonRobust = -2 * log(dmst(x = y, dp = dp))
    robust = ifelse(nonRobust > k, Psi(nonRobust, k = k), nonRobust)
    return(sum(robust))
}