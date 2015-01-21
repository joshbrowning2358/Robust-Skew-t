##' Wrapper to fastTLE assuming normal data.
##' 
##' @param data The data used in the analysis.
##' @param k The proportion of observations used.  If k=1, the (non-robust) MLE
##' will be computed.
##' @param initial A vector of length two containing an initial guess for the
##' mean and standard deviation.
##' 
##' @return Same as fastTLE.
##' 
##' @export
##' 

fast.TLE.normal = function(data, k, initial=c(0,1)){
    #Data quality checks
    if(!is(data,"matrix")){
        warning("data is not a matrix.  Attempting to coerce...")
        data = matrix(data)
    }
    if(ncol(data)!=1)
        stop("This function is for univariate normal only!")
    if(k>1 | k<=0)
        stop("k must be in (0,1]")
    if(length(initial)!=2)
        stop("initial must be of length 2!")
    
    MLE = function(data, start){ c(mean(data), sd(data)*NROW(data)/(NROW(data)-1)) }
    negLogLik = function(params, data){
        mu = params[1]
        sigma = params[2]
        log(sigma) + (data-mu)^2/sigma^2
    }
    fastTLE(initial, MLE, negLogLik, data, k)
}
