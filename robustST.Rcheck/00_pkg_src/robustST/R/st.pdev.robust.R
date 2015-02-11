##' Computes the robust-ified, penalized deviance for the univariate skew-t.
##' 
##' @param dp "Density parameters"?  Vector of xi, omega, alpha, nu.
##' @param x matrix of the independent variables for the fit.  Typically just a
##' matrix of ones, but can theoretically be a design matrix if the end goal is
##' fitting a regression with skew-t errors.
##' @param y A vector of dependent variables.
##' @param k Parameter controlling the robustness of the fit.  The largest
##' possible value for the negative log-likelihood is 2*k, and the negative
##' log-likelihood is adjusted down whenever it is larger than k.
##' 
##' @return The "robust" skew-t deviance evaluated at observations (x,y) and
##' with parameters dp.
##' 
##' @export
##' 

st.pdev.robust = function(dp, x, y, k = 2, ...){
    nonRobust = sapply(1:length(y), sn:::st.pdev, y = y[i],
                       x = x[i, , drop = FALSE], dp = dp, w = 1)
    robust = ifelse(nonRobust > k, sapply(nonRobust, Psi, k = k), nonRobust)
    return(sum(robust))
}