mst.pdev <- function (param, x, y, w, fixed.nu = NULL, symmetr = FALSE, penalty = NULL, 
    trace = FALSE) 
{
    if (missing(w)) 
        w <- rep(1, nrow(y))
    d <- NCOL(y)
    p <- ncol(x)
    npar0 <- (p * d + d * (d + 1)/2)
    param1 <- c(param[1:npar0], if (symmetr) rep(0, d) else param[npar0 + 
        (1:d)], if (is.null(fixed.nu)) param[length(param)])
    dp.list <- sn:::optpar2dplist(param1, d, p)
    dp <- dp.list$dp
    nu <- if (is.null(fixed.nu)){
        dp$nu
    } else {
        fixed.nu
    }
    logL <- sum(w * dmst(x = y, xi = x %*% dp$beta, Omega = dp$Omega,
                         alpha = dp$alpha, nu = nu, log = TRUE))
    Q <- if (is.null(penalty)){
        0
    } else {
        penalty(list(alpha = dp$alpha, Omega.bar = cov2cor(dp$Omega)), 
                nu, der = 0)
    }
    pdev <- (-2) * (logL - Q)
    if (trace) 
        cat("mst.pdev: ", pdev, "\nparam:", format(param), "\n")
    pdev
}
