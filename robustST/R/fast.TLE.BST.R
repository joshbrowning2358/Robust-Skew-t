##' Wrapper to fastTLE assuming bivariate skew-t data.
##' 
##' @param data The data used in the analysis.
##' @param k The proportion of observations used.  If k=1, the (non-robust) MLE
##' will be computed.
##' @param initial A vector of length eight containing an initial guess for
##' xi1, xi2, omega11, omega12, omega22, alpha1, alpha2, and nu.  Defaults to a
##' zero mean, identity omega, no skew, and large nu distribution (i.e.
##' approximately N(0, I)).  Note: do not set the nu parameter to Inf, as this
##' causes issues with the optimization.
##' 
##' @return Same as fastTLE.
##' 
##' @export
##' 

fast.TLE.BST = function(data, k, initial=c(0,0,1,0,1,0,0,1000), ...){
    #Data quality checks
    if(!is(data,"matrix")){
        warning("data is not a matrix.  Attempting to coerce...")
        data = matrix(data)
    }
    if(ncol(data)!=2)
        stop("This function is for bivariate skew-t only!")
    if(k>1 | k<=0)
        stop("k must be in (0,1]")
    if(length(initial)!=8)
        stop("initial must be of length 8 (xi1, xi2, om11, om12, om22, al1, al2, nu)!")
    
    MLE = function(data, start){
        fit = robustST(x=matrix(1,nrow(data)), y=data, robust=F, start=start)
        return( c(fit$beta[1], fit$beta[2], fit$Omega[1,1], fit$Omega[1,2], fit$Omega[2,2]
                  ,fit$alpha[1], fit$alpha[2], fit$nu) )
    }
    negLogLik = function(params, data){
        #We have to convert parameters to the type needed by mst.pdev:
        dp = list( beta=params[1:2], Omega=matrix(params[c(3,4,4,5)],nrow=2)
                   ,alpha=params[6:7], nu=params[8])
        params = dplist2optpar(dp)
        
        sapply(1:nrow(data), function(i){
            mst.pdev(y=data[i,,drop=F], param=params, x=matrix(1))
        })
    }
    fit = fastTLE(initial, MLE, negLogLik, data, k)
    return( list(xi=fit$MLE[1:2], omega=matrix(fit$MLE[c(3,4,4,5)],nrow=2)
                 ,alpha=fit$MLE[6:7], nu=fit$MLE[8]) )
}
