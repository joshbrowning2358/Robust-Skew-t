mst.pdev.grad <- function (param, y, x = matrix(1, NROW(y)),
                           w = rep(1, NROW(y)), fixed.nu = NULL,
                           symmetr = FALSE)
{
    d <- NCOL(y)
    p <- ncol(x)
    beta <- matrix(param[1:(p * d)], p, d)
    D <- exp(-2 * param[(p * d + 1):(p * d + d)])
    A <- diag(d)
    i0 <- p * d + d * (d + 1)/2
    if (d > 1) 
        A[!lower.tri(A, diag = TRUE)] <- param[(p * d + d + 1):i0]
    eta <- if (symmetr){
        rep(0, d)
    } else {param[(i0 + 1):(i0 + d)]}
    nu <- if (is.null(fixed.nu)){
        exp(param[length(param)])
    } else {fixed.nu}
    Oinv <- t(A) %*% diag(D, d, d) %*% A
    u <- y - x %*% beta
    u.w <- u * w
    Q <- as.vector(rowSums((u %*% Oinv) * u.w))
    L <- as.vector(u.w %*% eta)
    sf <- if (nu < 10000) {
        sqrt((nu + d)/(nu + Q))
    } else {sqrt((1 + d/nu)/(1 + Q/nu))}
    t. <- L * sf
    dlogft <- (-0.5) * sf^2
    dt.dL <- sf
    dt.dQ <- (-0.5) * L * sf/(Q + nu)
    logT. <- pt(t., nu + d, log.p = TRUE)
    dlogT. <- exp(dt(t., nu + d, log = TRUE) - logT.)
    Dbeta <- (-2 * t(x) %*% (u.w * dlogft) %*% Oinv - outer(as.vector(t(x) %*% 
        (dlogT. * dt.dL * w)), eta) - 2 * t(x) %*% (dlogT. * 
        dt.dQ * u.w) %*% Oinv)
    Deta <- colSums(dlogT. * sf * u.w)
    if (d > 1) {
        M <- 2 * (diag(D, d, d) %*% A %*% t(u * dlogft + u * 
            dlogT. * dt.dQ) %*% u.w)
        DA <- M[!lower.tri(M, diag = TRUE)]
    } else DA <- NULL
    M <- (A %*% t(u * dlogft + u * dlogT. * dt.dQ) %*% u.w %*% 
        t(A))
    if (d > 1){
        DD <- diag(M) + 0.5 * sum(w)/D
    } else DD <- as.vector(M + 0.5 * sum(w)/D)
    grad <- (-2) * c(Dbeta, DD * (-2 * D), DA, if (!symmetr) Deta)
    if (is.null(fixed.nu)) {
        df0 <- min(nu, 1e+08)
        if (df0 < 10000) {
            diff.digamma <- digamma((df0 + d)/2) - digamma(df0/2)
            log1Q <- log(1 + Q/df0)
        } else {
            diff.digamma <- log1p(d/df0)
            log1Q <- log1p(Q/df0)
        }
        dlogft.ddf <- 0.5 * (diff.digamma - d/df0 + (1 + d/df0) * 
            Q/((1 + Q/df0) * df0) - log1Q)
        eps <- 1e-04
        df1 <- df0 + eps
        sf1 <- if (df0 < 10000){
            sqrt((df1 + d)/(Q + df1))
        } else sqrt((1 + d/df1)/(1 + Q/df1))
        logT.eps <- pt(L * sf1, df1 + d, log.p = TRUE)
        dlogT.ddf <- (logT.eps - logT.)/eps
        Ddf <- sum((dlogft.ddf + dlogT.ddf) * w)
        grad <- c(grad, -2 * Ddf * df0)
    }
    return(grad)
}

