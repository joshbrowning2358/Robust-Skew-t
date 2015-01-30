##' Fit a robust Skew-t
##' 
##' Fits a robust version of the multivariate skew-t, done by bounding the
##' negative log-likelihood for each observation.
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
##' @param alpha A parameter controlling the "robustness" of the fit.  Given
##' current parameter estimates, a (1-alpha)% confidence region can be
##' constructed, and observations in this region will not be adjusted during
##' the optimization.  However, values outside this region will have their
##' likelihood adjusted down, and hence will have less influence on the
##' M-estimator.  As alpha->1 the estimator approaches the MLE.
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

robustST = function(y, x = matrix(1, nrow = NROW(y)), robust = T,
                    method = c("nlminb", "constrOptim"), w = rep(1,nrow(x)),
                    alpha = 0.01, start = NULL){
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
    
    #Treat univariate separately
    if(d == 1){
        if(is.null(start)){
            #Determine starting estimate (via logic from sn::st.mple function)
            ls <- lm.wfit(x, y, w)
            res <- ls$residuals
            s <- sqrt(sum(w * res^2)/nw)
            gamma1 <- sum(w * res^3)/(nw * s^3)
            gamma2 <- sum(res^4)/(nw * s^4) - 3
            cp <- c(ls$coef, s, gamma1, gamma2)
            dp <- st.cp2dp(cp, silent = TRUE)
            if (is.null(dp)) 
                dp <- rep(NA, length(cp))
            if (any(is.na(dp))) 
                dp <- c(cp[1:(p + 1)], 0, 10)
            names(dp) = c("xi", "omega", "alpha", "nu")
        } else {
            dp = start
        }
        
        if(!robust & method == "nlminb"){
            nlmEst = try(nlminb(start = dp,
                                function(dp){st.pdev(dp, x, y, w = w)},
                                gradient = function(dp){st.pdev.gh(dp, x, y)}))
            if(is(nlmEst, "try-error")){
                fit = rep(NA, length(dp) + 1)
            } else {
                fit = c(nlmEst$par, convergence = nlmEst$convergence)
            }
        }
        
        if(!robust & method == "constrOptim"){
            conEst = try(constrOptim(theta = dp,
                                     f = function(dp){st.pdev(dp, x, y)},
                                     grad = function(dp){st.pdev.gh(dp, x, y)},
                                     #Constraints: force omega>0 and nu>0
                                     ui = matrix(c(0,0,1,0,0,0,0,1), nrow=2),
                                     ci = rep(0,2)))
            if(is(conEst, "try-error")){
                fit = rep(NA, length(dp) + 1)
            } else {
                fit = c(conEst$par, convergence = conEst$convergence)
            }
        }
        
        if(robust & method == "nlminb"){
            nlmRobEst = try(nlminb( start = dp,
                                    function(dp){
                                        st.pdev.robust(dp, x, y, k = k)},
                                    gradient = function(dp){
                                        st.pdev.gh.robust(dp, x, y, k = k)}))
            if(is(nlmRobEst, "try-error")){
                fit = rep(NA, length(dp) + 1)
            } else {
                fit = c(nlmRobEst$par, convergence = nlmRobEst$convergence)
            }
        }
        
        if(robust & method == "constrOptim"){
            conRobEst = try(constrOptim(theta = dp,
                                        f = function(dp){
                                            st.pdev.robust(dp, x, y, k = k)},
                                        grad = function(dp){
                                            st.pdev.gh.robust(dp, x, y, k = k)},
                                        ui = matrix(c(0, 0, 1, 0, 0, 0, 0, 1),
                                                    nrow = 2),
                                        ci = rep(0,2)))
            if(is(conRobEst, "try-error")){
                fit = rep(NA, length(dp) + 1)
            } else {
                fit = c(conRobEst$par, convergence = conRobEst$convergence)
            }
        }
        
        return(fit)
        
    } else { #Now multivariate case
        if(is.null(start)){
            #Determine starting estimate (via logic from mst.mple function in sn)
            ls <- lm.wfit(x, y, w, singular.ok = FALSE)
            beta <- coef(ls)
            Omega <-  var(resid(ls))
            omega <- sqrt(diag(Omega))
            alpha <- rep(0, d)
            nu <- 8
            param <- dplist2optpar(list(beta = beta, Omega = Omega,
                                        alpha = alpha))
            param <- c(param, log(nu))
        } else {
            param = start
        }
        
        if(!robust & method=="nlminb"){
            nlmEst = try(nlminb( start=param
                                 ,function(param){mst.pdev(param, x, y)}
                                 ,gradient=function(dp){mst.pdev.grad(param, x, y, w=rep(1,NROW(x)))}))
            if(is(nlmEst, "try-error")){
                fit = rep(NA, length(param)+1)
            } else {
                fit = c(nlmEst$par, nlmEst$convergence)
            }
        }
        
        if(!robust & method=="constrOptim"){
            conEst = try(constrOptim(theta=param
                                     ,f=function(param){mst.pdev(param, x, y, w=w)}
                                     ,grad=function(param){mst.pdev.grad(param, x, y, w=w)}
                                     #No need for constraints as optpar2dplist ensures nu>0 and Omega is pos. def.
                                     #So, set u_i to all 0's, and force this to always be greater than -1
                                     # (which it always will).
                                     ,ui=matrix(0,ncol=length(param))
                                     ,ci=-1))
            if(is(conEst, "try-error")){
                fit = rep(NA, length(param)+1)
            } else {
                fit = c(conEst$par, conEst$convergence)
            }
        }
        
        if(robust & method=="nlminb"){
            nlmRobEst = try(nlminb( start=param
                                    ,function(param){mst.pdev.robust(param, x, y, k=k)}
                                    ,gradient=function(param){mst.pdev.grad.robust(param, x, y, k=k)}))
            if(is(nlmRobEst, "try-error")){
                fit = rep(NA, length(param)+1)
            } else {
                fit = c( nlmRobEst$par, nlmRobEst$convergence )
            }
        }
        
        if(robust & method=="constrOptim"){
            conRobEst = try(constrOptim(theta=param
                                        ,f=function(param){mst.pdev.robust(param, x, y, k=k)}
                                        ,grad=function(param){mst.pdev.grad.robust(param, x, y, k=k)}
                                        #No need for constraints as optpar2dplist ensures nu>0 and Omega is pos. def.
                                        #So, set u_i to all 0's, and force this to always be greater than -1
                                        # (which it always will).
                                        ,ui=matrix(0,ncol=length(param))
                                        ,ci=-1))
            if(is(conRobEst, "try-error")){
                fit = rep(NA, length(param)+1)
            } else {
                fit = c(conRobEst$par, conRobEst$convergence)
            }
        }
        
        optpar = optpar2dplist(fit[-length(fit)], p = p, d = d)
        return(list(beta = optpar$beta, Omega = optpar$Omega,
                    alpha = optpar$alpha, nu = optpar$nu,
                    convergence = fit[length(fit)]))
    }
}