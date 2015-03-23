##' Fit a robust Skew-t, estimate k once
##' 
##' Fits a robust version of the multivariate skew-t, done by bounding the
##' negative log-likelihood for each observation.  Compute the value of k once
##' based on the initial (non-robust) estimates of the parameters.
##' 
##' @param y A vector or matrix of observations to fit the skew-t to.
##' @param method: constrOptim uses a constrained algorithm, forcing nu and
##' omega>0.  However, the implementation for multivariate skew-t fitting
##' enforces this by default, so nlminb and constrOptim should be very similar
##' for multivariate data.  For univariate, constrOptim is recommended.  For
##' multivariate, constrOptim is also recommended as it appears to be faster.
##' @param w A vector of case weights, defaults to a vector of ones.
##' @param pValue A parameter controlling the "robustness" of the fit.  Given
##' current parameter estimates, a (1-pValue)% confidence region can be
##' constructed, and observations in this region will not be adjusted during
##' the optimization.  However, values outside this region will have their
##' likelihood adjusted down, and hence will have less influence on the
##' M-estimator.  As pValue->1 the estimator approaches the MLE.
##' @param start The starting values for the optimization.  If NULL, reasonable
##' values are automatically chosen.
##' @param cooptimize If TRUE, then the k-value is updated during optimization.
##' Essentially, the density function being optimized varies with k.
##' 
##' @return A named list containing the results of the fit.  xi/omega/alpha/nu
##' are the parameters of the skew-t.  A convergence flag is also returned,
##' indicating if the solution is a true optimum.
##' 
##' @export
##' 

robustSTOnceK = function(y, family = c("ST", "SN", "T", "N"),
                         method = c("nlminb", "constrOptim"),
                         w = rep(1, NROW(y)), pValue = 0.01, start = NULL){
    
#################################### TO DO ####################################
# - family is currently only implemented for ST
# - getLogLikelihoodBound is currently only implemented for d = 2.
# - implement cooptimization
###############################################################################

    ## Data quality checks
    if(any(is.na(y))){
        if(is.null(dim(y)))
            filt = !is.na(y)
        else
            filt = !apply(y, 1, function(x){any(is.na(x))})
        
        w = w[filt]
        if(is.null(dim(y)))
            y = y[filt]
        else
            y = y[filt, ]
    }
    if(!is.matrix(y) & !is.numeric(y))
        stop("y must be a matrix or numeric vector!")
    if(!is.numeric(w))
        stop("w must be numeric!")
    if(length(method) > 1)
        method = method[1]
    if(!method %in% c("nlminb", "constrOptim"))
        stop("method must be one of nlminb or constrOptim!")
    if(length(family) > 1)
        family = family[1]
    
<<<<<<< HEAD
    ## Get the density functions based on the family
    func = getDensityFunction(family = family, dimension = NCOL(y), robust = T)
    densityFunction = func[[1]]
    gradientFunction = func[[2]]
    
    ## If cooptimize, update densityFunction
    if(cooptimize){
        oldDensityFunction = densityFunction
        densityFunction = function(param, x, y, k, ...){
            ## Passed k value is ignored, instead the update is computed
            k = getLogLikelihoodBound(
                dp = sn:::optpar2dplist(param, d = d, p = p)$dp,
                alpha = pValue)
            oldDensityFunction(param = param, x = x, y = y, k = k, ...)
        }
    }
    
    ## Assign useful variables
    n = NROW(y)
    p = 1
    d = NCOL(y)
    nw = sum(w)
    
    ## Compute starting estimate if needed
    if(is.null(start)){
        dp = getStartingEstimate(y = y, family = family, w = w)
    } else {
        dp = start
    }
    param = dplist2optpar(dp)
    
    ## Fit the models
    k = getLogLikelihoodBound(dp = dp, alpha = pValue)
    if(method == "nlminb"){
        fit = try(nlminb(start = param,
                    function(param){
                        densityFunction(param = param, y = y, k = k)},
                    gradient = function(param){
                        gradientFunction(param = param, y = y, k = k)}))
    } else if(method == "constrOptim"){
        fit = try(constrOptim(theta = param ,f = function(dp){
                        densityFunction(param = param, y = y, k = k)},
                    grad = function(dp){
                        gradientFunction(param = param, y = y, k = k)}
            # No need for constraints as optpar2dplist ensures nu>0 and
            # Omega is pos. def. So, set u_i to all 0's, and force this
            # to always be greater than -1 (which it always will be).
            ,ui=matrix(0,ncol=length(param))
            ,ci=-1))
    } else {
        stop("Invalid values of method and/or d!  method must be ",
             "'constrOptim' or 'nlminb' and d must be 1 or 2")
    }
    if(is(fit, "try-error")){
        output = list(parameters = rep(NA, length(param)),
                      convergence = NA, kValues = NA)
    }
    output = list(parameters = fit$par,
                  convergence = fit$convergence,
                  kValues = k)
            
    ## Post-process the parameters and output them
    optpar = sn:::optpar2dplist(output$parameters, p = p, d = d)
    output$xi = optpar$beta
    output$Omega = optpar$Omega
    output$alpha = optpar$alpha
    output$nu = optpar$nu
    return(output)
}
