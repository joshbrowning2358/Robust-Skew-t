##' Gradient of robust, penalized deviance for multivariate skew-t
##' 
##' This function computes the gradient of the robust-ified, penalized deviance
##' for the multivariate skew-t (with respect to the optimization parameters).
##' 
##' @param param Optimization parameters derived from xi, Omega, alpha, nu. For
##' conversion, see the functions optpar2dplist and dplist2optpar.
##' @param x A matrix of the independent variables for the fit.  Typically just
##' a matrix of ones.
##' @param y A matrix of dependent variables.
##' @param k Parameter controlling the robustness of the fit.  The largest
##' possible value for the negative log-likelihood is 2*k, and the negative
##' log-likelihood is adjusted down whenever it is larger than k.
##' 
##' @return The gradient of the robust deviance with respect to dp.
##' 
##' @export
##' 

mst.pdev.grad.robust = function(param, y, k = 2, ...){
    nonRobustGH = mst.pdev.grad.vec(param = param, y = y, ...)
    dp = sn:::optpar2dplist(param, p = 1, d = NCOL(y))$dp
    nonRobust = -2 * dmst(x = y, dp = dp)
    # If -LogLikelihood<k, we don't adjust at all.  Thus, gradient doesn't
    # change for that case.
    # If -LogLikelihood>k, we change -LL to Psi(-LL).  Thus, gradient becomes
    # Psi'(-LL)*(-LL)'=psi(-LL)*(-LL)'
    robustGH = nonRobustGH
    if(any(nonRobust > k))
        robustGH[nonRobust > k, ] = robustGH[nonRobust > k, ] *
            sapply(nonRobust[nonRobust > k], psi.grad, k = k)
    #Some NA's can occur from st.pdev.gh.  Set them to 0, since I have no better idea as to why they occur (presumably really small density values)
    robustGH[is.na(robustGH)] = 0
    return(colSums(robustGH))
}