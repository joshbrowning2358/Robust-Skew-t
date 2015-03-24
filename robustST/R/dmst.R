dmst <- function (x, xi = rep(0, length(alpha)), Omega, alpha, nu = Inf, 
    dp = NULL, log = FALSE) 
{
    if (!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
        stop("You cannot set both component parameters and dp")
    if (!is.null(dp)) {
        if (length(dp) != 4) 
            stop("wrong length of non-null 'dp'")
        xi <- drop(dp[[1]])
        Omega <- dp[[2]]
        alpha <- dp[[3]]
        nu <- dp[[4]]
    }
    if (any(abs(alpha) == Inf)) 
        stop("Inf's in alpha are not allowed")
    if (nu == Inf) 
        return(sn:::dmsn(x, xi, Omega, alpha, log = log))
    d <- length(alpha)
    Omega <- matrix(Omega, d, d)
    if (!all(Omega - t(Omega) == 0)) 
        return(NA)
    invOmega <- mnormt:::pd.solve(Omega, silent = TRUE, log.det = TRUE)
    if (is.null(invOmega)) 
        return(NA)
    logDet <- attr(invOmega, "log.det")
    x <- if(is.vector(x)){
        matrix(x, ncol = 1, nrow = length(x))
    } else {
        data.matrix(x)
    }
    if (is.vector(xi)) 
        xi <- outer(rep(1, nrow(x)), xi)
    X <- t(x - xi)
    Q <- colSums((invOmega %*% X) * X)
    L <- as.vector(t(X/sqrt(diag(Omega))) %*% as.matrix(alpha))
    if (nu < 10000) {
        log.const <- lgamma((nu + d)/2) - lgamma(nu/2) - 0.5 * 
            d * logb(nu)
        log1Q <- logb(1 + Q/nu)
    } else {
        log.const <- (-0.5 * d * logb(2) + log1p((d/2) * (d/2 - 
            1)/nu))
        log1Q <- log1p(Q/nu)
    }
    log.dmt <- log.const - 0.5 * (d * logb(pi) + logDet + (nu + 
        d) * log1Q)
    log.pt <- pt(L * sqrt((nu + d)/(Q + nu)), df = nu + d, log.p = TRUE)
    logPDF <- logb(2) + log.dmt + log.pt
    if (log){
        return(logPDF)
    } else {
        return(exp(logPDF))
    }
}

