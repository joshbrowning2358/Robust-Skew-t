##' Gradient of the robust-ified, penalized deviance for the univariate skew-t.
##' 
##' @param dp "Density parameters"?  Vector of xi, omega, alpha, nu.
##' @param x A matrix of the independent variables for the fit.  Typically just
##' a matrix of ones.
##' @param y A vector of dependent variables
##' @param k Parameter controlling the robustness of the fit.  The largest
##' possible value for the negative log-likelihood is 2*k, and the negative
##' log-likelihood is adjusted down whenever it is larger than k.
##' 
##' @return The gradient of the robust deviance with respect to dp.
##' 
##' @export
##' 

st.pdev.gh.robust = function(dp, x, y, k=2){
    nonRobustGH = lapply(y, st.pdev.gh, dp=dp, x=matrix(1,nrow=1))
    nonRobustGH = do.call("rbind", nonRobustGH)
    nonRobust = lapply(y, st.pdev, dp=dp, x=matrix(1,nrow=1))
    nonRobust = do.call("c", nonRobust)
    #If -LogLikelihood<k, we don't adjust at all.  Thus, gradient doesn't change for that case.
    #If -LogLikelihood>k, we change -LL to Psi(-LL).  Thus, gradient becomes Psi'(-LL)*(-LL)'=psi(-LL)*(-LL)'
    robustGH = nonRobustGH
    if(any(nonRobust>k))
        robustGH[nonRobust>k,] = robustGH[nonRobust>k,]*sapply(nonRobust[nonRobust>k], psi.grad, k=k)
    #Some NA's can occur from st.pdev.gh.  Set them to 0, since I have no better idea as to why they occur (presumably really small density values)
    robustGH[is.na(robustGH)] = 0
    return( colSums(robustGH) )
}