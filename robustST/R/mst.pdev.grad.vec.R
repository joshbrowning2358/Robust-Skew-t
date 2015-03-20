##' Gradient of the Multivariate Skew-t Deviance
##' 
##' This function computes the gradient of the multivariate skew-t deviance for
##' each observation.
##' 
##' The mst.pdev.grad function is implemented in the sn package.  However, this
##' gradient is computed from the total deviance (i.e. considering all the
##' data) with respect to the parameters of the skew-t.  For the robust
##' implementation, we need the gradient of the deviance for each individual
##' observation.  This could be computed by applying mst.pdev.grad to each
##' observation, but this re-written function should be more efficient.
##' 
##' @param param The optimization parameters.  See sn:::optpar2dplist.
##' @param y The observed data values.  This should be a matrix of nxp, where
##' n = number of observations and p = number of dimensions of the data.
##' @param fixed.nu If NULL, the gradient with respect to nu will be estimated.
##' If a value is specified, that value will be used for nu and no gradient 
##' will be calculated for nu.
##' @param symmetr Should the distribution be assumed to be symmetric?  If so,
##' alpha is fixed and no gradient is computed for those parameters.
##' 
##' @return A matrix of dimensions length(param) x n.  The (i,j)-th entry
##' corresponds to the derivative of the deviance for observation i with
##' respect to the param variable j.
##' 

mst.pdev.grad.vec <- function (param, y, fixed.nu = NULL, symmetr = FALSE) 
{
    d <- ncol(y)
    p <- 1
    beta <- matrix(param[1:(p * d)], p, d)
    D <- exp(-2 * param[(p * d + 1):(p * d + d)])
    A <- diag(d)
    i0 <- p * d + d * (d + 1)/2
    if (d > 1) 
        A[!lower.tri(A, diag = TRUE)] <- param[(p * d + d + 
            1):i0]
    eta <- if (symmetr) 
        rep(0, d)
    else param[(i0 + 1):(i0 + d)]
    nu <- if (is.null(fixed.nu)) 
        exp(param[length(param)])
    else fixed.nu
    Oinv <- t(A) %*% diag(D, d, d) %*% A
    u <- y - x %*% beta
    Q <- as.vector(rowSums((u %*% Oinv) * u))
    L <- as.vector(u %*% eta)
    sf <- if (nu < 10000) 
        sqrt((nu + d)/(nu + Q))
    else sqrt((1 + d/nu)/(1 + Q/nu))
    t. <- L * sf
    dlogft <- (-0.5) * sf^2
    dt.dL <- sf
    dt.dQ <- (-0.5) * L * sf/(Q + nu)
    logT. <- pt(t., nu + d, log.p = TRUE)
    dlogT. <- exp(dt(t., nu + d, log = TRUE) - logT.)
    Dbeta <- (-2 * (u * dlogft) %*% Oinv - outer(as.vector( 
        (dlogT. * dt.dL)), eta) - 2 * (dlogT. * 
        dt.dQ * u) %*% Oinv)
    Deta <- dlogT. * sf * u
    if (d > 1) {
        M <- 2 * (diag(D, d, d) %*% A %*% t(u * dlogft + u * 
            dlogT. * dt.dQ))
        DA = NULL
        for(i in 1:(nrow(M)-1))
            for(j in (i+1):nrow(M))
                DA = cbind(DA, M[i, ] * u[, j])
    }
    else DA <- NULL
    M <- (A %*% t(u * dlogft + u * dlogT. * dt.dQ))
    DD = NULL
    for(i in 1:nrow(M))
        DD = cbind(DD, M[i, ] * (u %*% t(A))[, i] * D[i] + .5 / D[i])
    grad <- (-2) * cbind(Dbeta, -2 * DD, DA, if (!symmetr) Deta)
    if (is.null(fixed.nu)) {
        df0 <- min(nu, 100000000)
        if (df0 < 10000) {
            diff.digamma <- digamma((df0 + d)/2) - digamma(df0/2)
            log1Q <- log(1 + Q/df0)
        }
        else {
            diff.digamma <- log1p(d/df0)
            log1Q <- log1p(Q/df0)
        }
        dlogft.ddf <- 0.5 * (diff.digamma - d/df0 + (1 + d/df0) * 
            Q/((1 + Q/df0) * df0) - log1Q)
        eps <- 0.0001
        df1 <- df0 + eps
        sf1 <- if (df0 < 10000) 
            sqrt((df1 + d)/(Q + df1))
        else sqrt((1 + d/df1)/(1 + Q/df1))
        logT.eps <- pt(L * sf1, df1 + d, log.p = TRUE)
        dlogT.ddf <- (logT.eps - logT.)/eps
        Ddf <- dlogft.ddf + dlogT.ddf
        grad <- cbind(grad, -2 * Ddf * df0)
    }
    return(grad)
}
