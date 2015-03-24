##' Multivariate skew-t deviance
##' 
##' @param param
##' @param y
##' @param x
##' @param w
##' @param ... This argument is only for internal purposes.  The arguments
##' passed here do nothing.
##' 
##' @return
##' 

mst.pdev = function (param, y, x = matrix(1, NROW(y)), w = rep(1, NROW(y)), ...)
{
    
    ## Data Quality Checks
    d = NCOL(y)
    paramDim = dimensionFromParamSize(length(param))
    if(d != paramDim)
        stop("parameter suggests a different dimension of y than what is ",
             "observed.")
    
    if (missing(w)) 
        w = rep(1, NROW(y))
    p = ncol(x)
    dp.list = sn:::optpar2dplist(param, d, p)
    dp = dp.list$dp
    logL = sum(w * dmst(x = y, xi = x %*% dp$beta, Omega = dp$Omega,
                         alpha = dp$alpha, nu = dp$nu, log = TRUE))
    pdev = -2 * logL
    pdev
}
