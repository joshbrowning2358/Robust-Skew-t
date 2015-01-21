##' Wrapper to fastTLE assuming univariate skew-t data.
##' 
##' @param data The data used in the analysis.
##' @param k The proportion of observations used.  If k=1, the (non-robust) MLE
##' will be computed.
##' @param initial A vector of length four containing an initial guess for xi,
##' omega, alpha, and nu.
##' 
##' @return Same as fastTLE.
##' 
##' @export
##' 

fast.TLE.ST = function(data, k, initial=c(0,1,0,1000), ...){
    #Data quality checks
    if(!is(data,"matrix")){
        warning("data is not a matrix.  Attempting to coerce...")
        data = matrix(data)
    }
    if(ncol(data)!=1)
        stop("This function is for univariate skew-t only!")
    if(k>1 | k<=0)
        stop("k must be in (0,1]")
    if(length(initial)!=4)
        stop("initial must be of length 4!")
    
    MLE = function(data, start){ robustST(x=matrix(1,nrow(data)), y=data, robust=F, start=start) }
    negLogLik = function(params, data)
        sapply(data, st.pdev, dp=params, x=matrix(1))
    fastTLE(initial, MLE, negLogLik, data, k)
}