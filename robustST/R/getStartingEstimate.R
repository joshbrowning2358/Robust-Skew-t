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

getStartingEstimate = function(y, family, w = rep(1, NROW(y))){
    x = matrix(1, nrow = NROW(y), ncol = 1)
    d = NCOL(y)
    #Determine starting estimate (via logic from sn::st.mple function)
    if(d == 1){
        nw = sum(w)
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
            dp <- c(cp[1:2], 0, 10)
        if(family == "N")
            return(list(xi = dp[1], Omega = matrix(dp[2])))
        if(family == "T")
            return(list(xi = dp[1], Omega = matrix(dp[2]), nu = dp[4]))
        if(family == "SN")
            return(list(xi = dp[1], Omega = matrix(dp[2]), alpha = dp[3]))
        if(family == "ST")
            return(list(xi = dp[1], Omega = matrix(dp[2]), alpha = dp[3],
                        nu = dp[4]))
    } else {
        ls <- lm.wfit(x, y, w, singular.ok = FALSE)
        beta <- as.numeric(coef(ls))
        Omega <-  var(resid(ls))
        omega <- sqrt(diag(Omega))
        alpha <- rep(0, d)
        nu <- 8
        if(family == "N")
            return(list(xi = beta, Omega = Omega))
        if(family == "T")
            return(list(xi = beta, Omega = Omega, nu = nu))
        if(family == "SN")
            return(list(xi = beta, Omega = Omega, alpha = alpha))
        if(family == "ST")
            return(list(xi = beta, Omega = Omega, alpha = alpha, nu = nu))
    }
}