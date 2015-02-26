##' Get Starting Estimate
##' 
##' This function is used to provide a reasonable starting estimate for the
##' skew-t parameters.  The logic replicates what is implemented in the sn
##' package (in particular, the st.mple function).
##' 
##' @param y The data for which the skew-t is being fit.
##' @param family The distribution of the data, should be one of "N", "T",
##' "SN", or "ST" for normal, Student's-t, skew normal, or skew-t,
##' respectively.
##' @param x The design matrix, if y is assumed to be residuals from some fit.
##' This parameter has not been tested thoroughly, so use with caution!  It
##' defaults to a matrix of 1's with 1 column, and this is equivalent to simply
##' fitting a skew-t to y.
##' @param w The observation weights used in the skew-t fitting.
##' 
##' @return An initial estimate for the density parameters.
##' 

getStartingEstimate = function(y, family,
                               x = matrix(1, nrow = NROW(y), ncol = 1),
                               w = rep(1, NROW(y))){
    d = NCOL(y)
    #Determine starting estimate (via logic from sn::st.mple function)
    if(d == 1){
        ls <- lm.wfit(x, y, w)
        res <- ls$residuals
        s <- sqrt(sum(w * res^2)/nw)
        gamma1 <- sum(w * res^3)/(nw * s^3)
        gamma2 <- sum(res^4)/(nw * s^4) - 3
        cp <- c(ls$coef, s, gamma1, gamma2)
        dp <- sn:::st.cp2dp(cp, silent = TRUE)
        if (is.null(dp)) 
            dp <- rep(NA, length(cp))
        if (any(is.na(dp))) 
            dp <- c(cp[1:(p + 1)], 0, 10)
        names(dp) = c("xi", "omega", "alpha", "nu")
        if(family == "N")
            return(dp[1:2])
        if(family == "T")
            return(dp[c(1,2,4)])
        if(family == "SN")
            return(dp[c(1,2,3)])
        if(family == "ST")
            return(dp)
    } else {
        ls <- lm.wfit(x, y, w, singular.ok = FALSE)
        beta <- coef(ls)
        Omega <-  var(resid(ls))
        omega <- sqrt(diag(Omega))
        alpha <- rep(0, d)
        nu <- 8
        param <- sn:::dplist2optpar(list(beta = beta, Omega = Omega,
                                         alpha = alpha))
        param <- c(param, log(nu))
        if(family == "N")
            return(param[1:2])
        if(family == "T")
            return(param[c(1,2,4)])
        if(family == "SN")
            return(param[c(1,2,3)])
        if(family == "ST")
            return(param)
    }
}