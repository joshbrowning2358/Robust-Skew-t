dmsn <- function (x, xi = rep(0, length(alpha)), Omega, alpha, tau = 0, 
    dp = NULL, log = FALSE) 
{
    if (!(missing(alpha) & missing(Omega)) && !is.null(dp)) 
        stop("You cannot set both component parameters and dp")
    if (!is.null(dp)) {
        if (length(dp) < 3) 
            stop("wrong length of non-null 'dp'")
        xi <- drop(dp[[1]])
        Omega <- dp[[2]]
        alpha <- dp[[3]]
        tau <- if (length(dp) == 4){
            dp[[4]]
        } else 0
    }
    if (any(abs(alpha) == Inf)) 
        stop("Inf's in alpha are not allowed")
    d <- length(alpha)
    Omega <- matrix(Omega, d, d)
    invOmega <- mnormt::pd.solve(Omega, silent = TRUE, log.det = TRUE)
    if (is.null(invOmega)) 
        stop("Omega matrix is not positive definite")
    logDet <- attr(invOmega, "log.det")
    x <- if (is.vector(x)){
        as.matrix(x)
    } else {data.matrix(x)}
    if (is.vector(xi)) 
        xi <- outer(rep(1, nrow(x)), xi)
    if (tau == 0) {
        log.const <- logb(2)
        alpha0 <- 0
    } else {
        log.const <- -pnorm(tau, log.p = TRUE)
        O.alpha <- cov2cor(Omega) %*% alpha
        alpha0 <- tau * sqrt(1 + sum(alpha * O.alpha))
    }
    X <- t(x - xi)
    Q <- colSums((invOmega %*% X) * X)
    L <- alpha0 + as.vector(t(X/sqrt(diag(Omega))) %*% as.matrix(alpha))
    logPDF <- (log.const - 0.5 * Q + pnorm(L, log.p = TRUE) - 
        0.5 * (d * logb(2 * pi) + logDet))
    if (log) 
        logPDF
    else exp(logPDF)
}
