##' Fit a robust Skew-t, iteratively update k
##' 
##' Fits a robust version of the multivariate skew-t, done by bounding the
##' negative log-likelihood for each observation.  Re-compute the value of k
##' each time new parameter estimates are available.
##' 
##' @param y A vector or matrix of observations to fit the skew-t to.
##' @param x A matrix of ones, or matrix of independent variables for skew-t
##' regression (use caution, as this feature has not been tested!)
##' @param robust Should the robust estimator be used?
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
##' 
##' @return A named list containing the results of the fit.  beta vector is
##' equivalent to the mean estimate if x = matrix of 1's, and omega/alpha/nu
##' are the parameters of the skew-t.  A convergence flag is also returned,
##' indicating if the solution is a true optimum.
##' 
##' @export
##' 

robustSTChangingK = function(y, x = matrix(1, nrow = NROW(y)), robust = T,
                             method = c("nlminb", "constrOptim"),
                             w = rep(1, nrow(x)), pValue = 0.01,
                             start = NULL){
    #library(sn)
    
    #Data quality checks
    if(any(is.na(y))){
        if(is.null(dim(y)))
            filt = !is.na(y)
        else
            filt = !apply(y, 1, function(x){any(is.na(x))})
        
        x = x[filt, ]
        w = w[filt]
        if(is.null(dim(y)))
            y = y[filt]
        else
            y = y[filt, ]
    }
    if(!is(x, "matrix"))
        stop("x must be a matrix!")
    if(!is.matrix(y) & !is.numeric(y))
        stop("y must be a matrix or numeric vector!")
    if(!is(robust,"logical"))
        stop("robust must be a logical!")
    if(nrow(x) != NROW(y))
        stop("x and y must have the same number of observations!")
    if(!is.numeric(w))
        stop("w must be numeric!")
    if(length(w) != nrow(x))
        stop("w must have the same length as ncol(x)!")
    if(length(method) > 1)
        method = method[1]
    if(!method %in% c("nlminb", "constrOptim"))
        stop("method must be one of nlminb or constrOptim!")
    
    n = nrow(x)
    p = ncol(x)
    d = NCOL(y)
    nw = sum(w)
    
    if(is.null(start)){
        getStartingEstimate(y = y, x = x, w = w)
    } else {
        param = start
    }
    
    kValues = NULL
    kDelta = Inf
    while(kDelta > .01){
        kOld = k
        k = getLogLikelihoodBound(
            dp = sn:::optpar2dplist(param, d = d, p = p)$dp,
            alpha = pValue)
        ## Compute the change in k, stop if small enough
        kDelta = abs(k - kOld)
        if(method == "nlminb" & d == 1){
            fit = try(nlminb( start = dp,
                        function(dp){
                            st.pdev.robust(dp, x, y, k = k)},
                        gradient = function(dp){
                            st.pdev.gh.robust(dp, x, y, k = k)}))
        } else if(method == "constrOptim" & d == 1){
            fit = try(constrOptim(theta = dp, f = function(dp){
                st.pdev.robust(dp, x, y, k = k)},
                grad = function(dp){
                    st.pdev.gh.robust(dp, x, y, k = k)},
                ui = matrix(c(0, 0, 1, 0, 0, 0, 0, 1),
                            nrow = 2),
                ci = rep(0,2)))
        } else if(method == "nlminb" & d > 1){
            fit = try(nlminb(start = param, function(param){
                    mst.pdev.robust(param, x, y, k= k)}, 
                gradient = function(param){
                    mst.pdev.grad.robust(param, x, y, k = k)}))
        } else if(method == "constrOptim" & d > 1){
            fit = try(constrOptim(theta=param ,f=function(param){
                mst.pdev.robust(param, x, y, k=k)},
                grad=function(param) mst.pdev.grad.robust(param, x, y, k=k)
                # No need for constraints as optpar2dplist ensures nu>0 and
                # Omega is pos. def. So, set u_i to all 0's, and force this
                # to always be greater than -1 (which it always will).
                ,ui=matrix(0,ncol=length(param))
                ,ci=-1))
        } else {
            stop("Invalid values of method and/or d!  method must be ",
                 "'constrOptim' or 'nlminb' and d must be 1 or 2")
        }
        if(is(fit, "try-error")){
            output = list(parameters = rep(NA, length(param)),
                          convergence = NA, kValues = NA)
            break
        }
        kValues = c(kValues, k)
    }
    output = list(parameters = fit$par,
                  convergence = fit$convergence,
                  kValues = kValues)
            
    ## Post-process the parameters and output them
    optpar = optpar2dplist(output$parameters, p = p, d = d)
    output$beta = optpar$beta
    output$Omega = optpar$Omega
    output$alpha = optpar$alpha
    output$nu = optpar$nu
    return(output)
}
